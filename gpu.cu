#include "gpu.h"
#include "common.h"

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
        grid[index].addParticle_host(i);
    }
    return grid;
}

// Step 3
// Push the grid to the gpu
__host__ Bin* push_data_to_device(Bin* grid){
    Bin* bins;
    GPUERRCHK(cudaMalloc((void **) &bins, NUM_BINS_PER_DIM * NUM_BINS_PER_DIM * sizeof(Bin)));
    GPUERRCHK(cudaMemcpy(bins, grid, NUM_BINS_PER_DIM * NUM_BINS_PER_DIM * sizeof(Bin), cudaMemcpyHostToDevice));
    return bins;
}

__host__ particle_t* push_particles_to_device(particle_t* particles, int n){
    particle_t* device_particles;
    GPUERRCHK(cudaMalloc((void **) &device_particles, n * sizeof(particle_t)));
    GPUERRCHK(cudaMemcpy(device_particles, particles, n * sizeof(particle_t), cudaMemcpyHostToDevice));
    return device_particles;
}

// Step 4
// Generate a dummy buffer for the bins
__host__ Bin* generateRedundantBins(){
    Bin* redundantBins;
    GPUERRCHK(cudaMalloc((void **) &redundantBins, NUM_BINS_PER_DIM * NUM_BINS_PER_DIM * sizeof(Bin)));
    return redundantBins;
}

// Step 0
// Clear the redundant bins
__global__ void clear_bins(Bin* redundantBins, int NUM_BINS_PER_DIM){
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= NUM_BINS_PER_DIM * NUM_BINS_PER_DIM){ 
        return;
    }
    redundantBins[tid].currentSize = 0;
}


__device__ void compute_force_grid(particle_t* device_particles, Bin* bins, int NUM_BINS_PER_DIM, int tid){

    int i = tid / NUM_BINS_PER_DIM;
    int j = tid % NUM_BINS_PER_DIM;

    int currentIndex = FIND_POS_DEVICE(i, j, NUM_BINS_PER_DIM);
    Bin& currentBin = bins[currentIndex];

    //set acceleration to zero
    for (int k = 0; k < currentBin.currentSize; k++){
        device_particles[currentBin.ids[k]].ax = device_particles[currentBin.ids[k]].ay = 0;
    }
    //check right
    if (j != NUM_BINS_PER_DIM - 1){
        compute_force_between_blocks(device_particles, currentBin, bins[FIND_POS_DEVICE(i, j+1, NUM_BINS_PER_DIM)]);
    }
    //check diagonal right bot
    if (j != NUM_BINS_PER_DIM - 1 && i != NUM_BINS_PER_DIM - 1){
        compute_force_between_blocks(device_particles, currentBin, bins[FIND_POS_DEVICE(i+1, j+1, NUM_BINS_PER_DIM)]);
    }
    //check diagonal right top
    if (j != NUM_BINS_PER_DIM - 1 && i != 0){
        compute_force_between_blocks(device_particles, currentBin, bins[FIND_POS_DEVICE(i-1, j+1, NUM_BINS_PER_DIM)]);
    }
    //check left
    if (j != 0){
        compute_force_between_blocks(device_particles, currentBin, bins[FIND_POS_DEVICE(i, j-1, NUM_BINS_PER_DIM)]);
    }
    //check diagonal left bot
    if (j != 0 && i != NUM_BINS_PER_DIM - 1){
        compute_force_between_blocks(device_particles, currentBin, bins[FIND_POS_DEVICE(i+1, j-1, NUM_BINS_PER_DIM)]);
    }
    //check diagonal left top
    if (j != 0 && i != 0){
        compute_force_between_blocks(device_particles,currentBin, bins[FIND_POS_DEVICE(i-1, j-1, NUM_BINS_PER_DIM)]);
    }
    //check top
    if (i != 0){
        compute_force_between_blocks(device_particles, currentBin, bins[FIND_POS_DEVICE(i-1, j, NUM_BINS_PER_DIM)]);
    }
    //check bot
    if (i != NUM_BINS_PER_DIM - 1){
        compute_force_between_blocks(device_particles, currentBin, bins[FIND_POS_DEVICE(i+1, j, NUM_BINS_PER_DIM)]);
    }
    //compute within itself
    compute_force_between_blocks(device_particles, currentBin, currentBin);
}

__device__ void move_particles(particle_t* device_particles, Bin* bins, Bin* redundantBins, double BIN_SIZE, int NUM_BINS_PER_DIM, double GRID_SIZE, int tid){
    Bin& bin = bins[tid];
    for (int k = 0; k < bin.currentSize; k++){
        move_gpu(&device_particles[bin.ids[k]], GRID_SIZE);
        bin_change(redundantBins, device_particles[bin.ids[k]], bin.ids[k], BIN_SIZE, NUM_BINS_PER_DIM);
    }
}

__global__ void compute_move_particles(particle_t* device_particles, Bin* bins, Bin* redundantBins, double BIN_SIZE, int NUM_BINS_PER_DIM, double GRID_SIZE){
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= NUM_BINS_PER_DIM * NUM_BINS_PER_DIM){ 
        return;
    }
    compute_force_grid(device_particles, bins, NUM_BINS_PER_DIM, tid);
    move_particles(device_particles, bins, redundantBins, BIN_SIZE, NUM_BINS_PER_DIM, GRID_SIZE, tid);
}



// Simulation begins
__host__ void simulate_particles(FILE* fsave, particle_t* particles, particle_t* device_particles, Bin* grid, Bin* bins, Bin* redundantBins, int n){

    int num_blocks = (NUM_BINS_PER_DIM * NUM_BINS_PER_DIM + NUM_THREADS - 1) / NUM_THREADS;

    for(int step = 0; step < NSTEPS; step++ ) {

        compute_move_particles <<< num_blocks, NUM_THREADS >>> (device_particles, bins, redundantBins, BIN_SIZE, NUM_BINS_PER_DIM, GRID_SIZE);

        swap(bins, redundantBins);

        clear_bins <<< num_blocks, NUM_THREADS >>> (redundantBins, NUM_BINS_PER_DIM);


        if( fsave && (step%SAVEFREQ) == 0 ) {
            // Copy the particles back to the CPU
            cudaMemcpy(particles, device_particles, n * sizeof(particle_t), cudaMemcpyDeviceToHost);
            save( fsave, n, particles);
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

    particle_t* device_particles = push_particles_to_device(particles, n);
    
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
    simulate_particles(fsave, particles, device_particles, grid, device_bins, device_redundant_bins, n);
    simulation_time = read_timer( ) - simulation_time;
    cudaThreadSynchronize();
    
    printf( "CPU-GPU copy time = %g seconds\n", copy_time);
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );
    
    free( particles );
    free(grid);
    cudaFree(device_bins);
    cudaFree(device_redundant_bins);

    if( fsave )
        fclose( fsave );
    
    return 0;
}
