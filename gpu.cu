#include "gpu.h"
#include "common.h"
#include <math.h> 
using namespace std;



// Step 1
// Generate constant variables 
__host__ void variables_initialization(int n){
    //initialize a bunck of stuff here
    GRID_SIZE = size;
    NUM_BINS_PER_DIM = int(ceil(sqrt(n)));
    BIN_SIZE = GRID_SIZE / NUM_BINS_PER_DIM;
    if (BIN_SIZE < cutoff){
        NUM_BINS_PER_DIM = int((GRID_SIZE /  cutoff));
        BIN_SIZE = GRID_SIZE / NUM_BINS_PER_DIM;
    }

    //debug
    if (DEBUG){
        printf("Grid Size is %f \n", GRID_SIZE);
        printf("NUM_BINS_PER_DIM is %d \n", NUM_BINS_PER_DIM);
        printf("BIN_SIZE is %f \n", BIN_SIZE);
        printf("Cutoff %f \n", cutoff);
        printf("dt %f \n", dt);
    }
    
}





// Step 2
// Generate the grid
__host__ Bin* generateGrid(particle_t* particles, int n){
    //initialize the grid with NUM_BINS_PER_DIM^2
    Bin *grid = (Bin *)malloc(NUM_BINS_PER_DIM * NUM_BINS_PER_DIM * sizeof(Bin)); 

    for (int i = 0; i < NUM_BINS_PER_DIM * NUM_BINS_PER_DIM; i++){
        grid[i] = Bin();
    }

    //store the point into the grid
    for (int i = 0; i < n; i++){
        int which_block_x = min((int)(particles[i].x / BIN_SIZE), NUM_BINS_PER_DIM - 1);
        int which_block_y = min((int)(particles[i].y / BIN_SIZE), NUM_BINS_PER_DIM - 1);
        int index = FIND_POS_HOST(which_block_y, which_block_x, NUM_BINS_PER_DIM);
        grid[index].addParticle(particles[i], i);
    }
    return grid;
}

// Step 3
// Push the grid to the gpu
__host__ Bin* push_data_to_device(Bin* grid){

    // First, allocate enough bins on GPUs
    Bin* bins;
    GPUERRCHK(cudaMalloc((void **) &bins, NUM_BINS_PER_DIM * NUM_BINS_PER_DIM * sizeof(Bin)));

    // Second, memcpy from grid to bins
    GPUERRCHK(cudaMemcpy(bins, grid, NUM_BINS_PER_DIM * NUM_BINS_PER_DIM * sizeof(Bin), cudaMemcpyHostToDevice));

    return bins;
}

// Step 4
// Generate a dummy buffer for the bins
__host__ Bin* generateRedundantBins(){
    Bin* redundantBins; // here is where we can change it to a 2d array? which then allows CUDA 2-d? 
    GPUERRCHK(cudaMalloc((void **) &redundantBins, NUM_BINS_PER_DIM * NUM_BINS_PER_DIM * sizeof(Bin)));
    return redundantBins;
}

// Step 0
// Clear the redundant bins
__global__ void clear_bins(Bin* redundantBins, int NUM_BINS_PER_DIM){
    int tid_x = threadIdx.x + blockIdx.x * blockDim.x;
    int tid_y = threadIdx.y + blockIdx.y * blockDim.y; 
    if (tid_x >= NUM_BINS_PER_DIM){ 
        return;
    }
    if (tid_y >= NUM_BINS_PER_DIM){ 
        return;
    }

    redundantBins[tid_y*NUM_BINS_PER_DIM + tid_x].currentSize = 0;
}


__device__ void compute_force_grid(Bin* bins, int NUM_BINS_PER_DIM, int i, int j, Bin* local_memory){



    //int currentIndex = FIND_POS_DEVICE(i, j, NUM_BINS_PER_DIM);
    //Bin& currentBin = bins[currentIndex];

    Bin& currentBin = local_memory[ LOCAL_BLOCK_INDEX(i , j) ];

    //set acceleration to zero
    for (int k = 0; k < currentBin.currentSize; k++){
        currentBin.particles[k].ax = currentBin.particles[k].ay = 0;
    }

    //check right
    if (j != NUM_BINS_PER_DIM - 1){
        compute_force_between_blocks(currentBin, local_memory[ LOCAL_BLOCK_INDEX(i , j+1) ]);
    }
    //check diagonal right bot
    if (j != NUM_BINS_PER_DIM - 1 && i != NUM_BINS_PER_DIM - 1){
        compute_force_between_blocks(currentBin, local_memory[ LOCAL_BLOCK_INDEX(i+1 , j+1) ] );
    }
    //check diagonal right top
    if (j != NUM_BINS_PER_DIM - 1 && i != 0){
        compute_force_between_blocks(currentBin, local_memory[ LOCAL_BLOCK_INDEX(i-1 , j+1) ] );
    }
    //check left
    if (j != 0){
        compute_force_between_blocks(currentBin, local_memory[ LOCAL_BLOCK_INDEX(i , j-1) ]);
    }
    //check diagonal left bot
    if (j != 0 && i != NUM_BINS_PER_DIM - 1){
        compute_force_between_blocks(currentBin, local_memory[ LOCAL_BLOCK_INDEX(i+1 , j-1) ] );
    }
    //check diagonal left top
    if (j != 0 && i != 0){
        compute_force_between_blocks(currentBin, local_memory[ LOCAL_BLOCK_INDEX(i-1 , j-1) ]);
    }
    //check top
    if (i != 0){
        compute_force_between_blocks(currentBin, local_memory[ LOCAL_BLOCK_INDEX(i-1 , j) ] );
    }
    //check bot
    if (i != NUM_BINS_PER_DIM - 1){
        compute_force_between_blocks(currentBin, local_memory[ LOCAL_BLOCK_INDEX(i+1, j) ]);
    }
    //compute within itself
    compute_force_within_block(currentBin);
}

__device__ void move_particles(Bin* bins, Bin* redundantBins, double BIN_SIZE, int NUM_BINS_PER_DIM, double GRID_SIZE, int i,int j, Bin* local_memory){
    Bin& bin = local_memory[ LOCAL_BLOCK_INDEX(i , j) ];
    for (int k = 0; k < bin.currentSize; k++){
        particle_t p = bin.particles[k];
        move_gpu(&p, GRID_SIZE);
        bin_change(redundantBins, p, bin.ids[k], BIN_SIZE, NUM_BINS_PER_DIM);
    }
}




__global__ void compute_move_particles(Bin* bins, Bin* redundantBins, double BIN_SIZE, int NUM_BINS_PER_DIM, double GRID_SIZE){
    int tid_x = threadIdx.x + blockIdx.x * blockDim.x;
    int tid_y = threadIdx.y + blockIdx.y * blockDim.y;
    if (tid_x >= NUM_BINS_PER_DIM){ 
        return;
    }
    if (tid_y >=NUM_BINS_PER_DIM ){
        return;
    }
    
    __shared__ Bin local_memory[18*18]; 

    int i = threadIdx.y;
    int j = threadIdx.x;

// this numbering might be wrong because i -> y, j -> x but this is wrong passing for tid_x, tid_y

    if( (1 <= i && i <= 14  && 1 <= j && j <= 14) ){
        local_memory[ LOCAL_BLOCK_INDEX(i, j) ] = *BINS_INDEX_RETURN(tid_x, tid_y, bins);
    }
    else if( (i == 0 || i == 15) && (j == 0  || j == 15)  ){
        local_memory[ LOCAL_BLOCK_INDEX(i,j) ] = *BINS_INDEX_RETURN(tid_x, tid_y, bins);

        int i_offset;
        int j_offset; 

        if(i ==0){
            i_offset = -1;
        }
        else{
            i_offset = 1;
        }
        if(j == 0){
            j_offset = -1;
        }
        else{
            j_offset = 1; 
        }
        local_memory[ LOCAL_BLOCK_INDEX(i+i_offset,j + j_offset) ] = *BINS_INDEX_RETURN(tid_x + i_offset, tid_y + j_offset, bins);
        local_memory[ LOCAL_BLOCK_INDEX(i+i_offset,j) ] = *BINS_INDEX_RETURN(tid_x+i_offset, tid_y, bins);
        local_memory[ LOCAL_BLOCK_INDEX(i,j+j_offset) ] = *BINS_INDEX_RETURN(tid_x, tid_y+j_offset, bins);
        
    }
    else if(i== 0 || i == 15){
        local_memory[ LOCAL_BLOCK_INDEX(i,j) ] = *BINS_INDEX_RETURN(tid_x, tid_y, bins);
        int i_offset;
        if(i ==0){
            i_offset = -1;
        }
        else{
            i_offset = 1;
        }
        local_memory[ LOCAL_BLOCK_INDEX(i + i_offset, j) ] =* BINS_INDEX_RETURN(tid_x+i_offset, tid_y, bins);

    } 
    else if (j == 0 || j == 15){
        local_memory[ LOCAL_BLOCK_INDEX(i,j) ] = *BINS_INDEX_RETURN(tid_x, tid_y, bins);
        int j_offset;
        if(j == 0){
            j_offset = -1;
        }
        else{
            j_offset = 1; 
        }
        local_memory[ LOCAL_BLOCK_INDEX(i , j+j_offset) ] = *BINS_INDEX_RETURN(tid_x, tid_y+j_offset, bins);
    }
    
    

    
    __syncthreads(); 

    compute_force_grid(bins, NUM_BINS_PER_DIM, i, j, local_memory);

    __syncthreads(); 
    move_particles(bins, redundantBins, BIN_SIZE, NUM_BINS_PER_DIM, GRID_SIZE, i, j, local_memory);

    __syncthreads(); 
}



// Simulation begins
__host__ void simulate_particles(FILE* fsave, particle_t* particles, Bin* grid, Bin* bins, Bin* redundantBins, int n){

    //int num_blocks = (NUM_BINS_PER_DIM * NUM_BINS_PER_DIM + NUM_THREADS - 1) / NUM_THREADS;



    // New square blocking

    dim3 threadsPerBlock(16, 16); // 256 threads 
    dim3 numBlocks(ceil(NUM_BINS_PER_DIM * 1.0/threadsPerBlock.x), ceil(NUM_BINS_PER_DIM*1.0/threadsPerBlock.y) );

    for(int step = 0; step < NSTEPS; step++ ) {

        compute_move_particles <<< numBlocks, threadsPerBlock >>> (bins, redundantBins, BIN_SIZE, NUM_BINS_PER_DIM, GRID_SIZE);

        swap(bins, redundantBins);

        clear_bins <<< numBlocks, threadsPerBlock >>> (redundantBins, NUM_BINS_PER_DIM);


        if( fsave && (step%SAVEFREQ) == 0 ) {
            // Copy the particles back to the CPU
            cudaMemcpy(grid, bins, NUM_BINS_PER_DIM * NUM_BINS_PER_DIM * sizeof(Bin), cudaMemcpyDeviceToHost);
            int count = 0;
            for (int p = 0; p < NUM_BINS_PER_DIM; p++){
                for (int q = 0; q < NUM_BINS_PER_DIM; q++){
                    Bin& thisBin = grid[FIND_POS_HOST(p,q, NUM_BINS_PER_DIM)];
                    for (int k = 0; k < thisBin.currentSize; k++){
                        count++;
                        particles[thisBin.ids[k]] = thisBin.particles[k];
                    }
                }
            }
            save( fsave, count, particles);
        }
    }
}


int main( int argc, char **argv )
{    
    // This takes a few seconds to initialize the runtime
    cudaThreadSynchronize(); 

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );


    //generate particles
    set_size( n );
    init_particles( n, particles );
    
    // Step 1 initialize constants and grids
    variables_initialization(n);

    // Step 2 generate grid
    Bin* grid = generateGrid(particles, n);
    
    // Get copy timestamp
    double copy_time = read_timer( );
    
    // Step 3 push the grid to the gpu
    Bin* device_bins = push_data_to_device(grid);
    
    // Step 4 generate redundant bins for devie_bins
    // We use it store data and later swap with device bins after each round
    Bin* device_redundant_bins = generateRedundantBins();
    
    // Barrier here
    cudaThreadSynchronize();
    copy_time = read_timer( ) - copy_time;
    

    //begin simulation
    cudaThreadSynchronize();
    double simulation_time = read_timer( );
    // Simulate Particles
    simulate_particles(fsave, particles, grid, device_bins, device_redundant_bins, n);
    simulation_time = read_timer( ) - simulation_time;
    cudaThreadSynchronize();
    
    //printf( "CPU-GPU copy time = %g seconds\n", copy_time);
    printf( "%d %g\n", n, simulation_time );
    
    free( particles );
    free(grid);
    cudaFree(device_bins);
    cudaFree(device_redundant_bins);

    if( fsave )
        fclose( fsave );
    
    return 0;
}