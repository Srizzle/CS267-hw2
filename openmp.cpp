
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <vector>
#include <set>
#include <iostream>
#include <omp.h>

using namespace std;

#define FIND_POS(ROW_INDEX, COL_INDEX, NUM_BLOCKS_PER_DIM) (ROW_INDEX * NUM_BLOCKS_PER_DIM + COL_INDEX)

//the length of each blocks
double BLOCK_SIZE = 0;
//the maximum span of the particles in both dimension
double GRID_SIZE = 1;
//number of blocks per dim, so the overall number of number of blocks is its square
int NUM_BLOCKS_PER_DIM = -1;

struct Serial_output{
    int nabsavg;
    double absavg;
    double absmin;
    Serial_output(){
        nabsavg = 0;
        absavg = 0;
        absmin = 1;
    }
};

vector<vector<set<int> > > grid;
vector<omp_lock_t> block_locks; //for each block, there is a lock



//Called in compute_force_grid
void compute_force_within_block(set<int>& block, particle_t* particles,
          int& navg, double& davg, double& dmin){

    for (set<int>::iterator it_1 = block.begin(); it_1 != block.end(); it_1++){
      for (set<int>::iterator it_2 = block.begin(); it_2 != block.end(); it_2++){
        apply_force(particles[*it_1], particles[*it_2], &dmin, &davg, &navg);
      }
    }
}

//Called in compute_force_grid
void compute_force_between_blocks(set<int>& block_A, set<int>& block_B, particle_t* particles,
          int navg, double davg, double dmin){
  for (set<int>::iterator it_A = block_A.begin(); it_A != block_A.end(); it_A++){
    for (set<int>::iterator it_B = block_B.begin(); it_B != block_B.end(); it_B++){
      apply_force(particles[*it_A], particles[*it_B], &dmin, &davg, &navg);
    }
  }
}

//Called in move_particles to change the membership
void move_to_another_block(int i, double old_x, double old_y,
                particle_t* particles){

      int which_block_x_old = min((int)(old_x / BLOCK_SIZE), NUM_BLOCKS_PER_DIM - 1);
      int which_block_y_old = min((int)(old_y / BLOCK_SIZE), NUM_BLOCKS_PER_DIM - 1);

      int which_block_x = min((int)(particles[i].x / BLOCK_SIZE), NUM_BLOCKS_PER_DIM - 1);
      int which_block_y = min((int)(particles[i].y / BLOCK_SIZE), NUM_BLOCKS_PER_DIM - 1);


      if (which_block_x_old != which_block_x || which_block_y_old != which_block_y){
        int old_block_index = FIND_POS(which_block_y_old, which_block_x_old, NUM_BLOCKS_PER_DIM);
        int new_block_index = FIND_POS(which_block_y, which_block_x, NUM_BLOCKS_PER_DIM);
        //lock, update, and release the old block
        omp_set_lock(&block_locks[old_block_index]);
        grid[which_block_y_old][which_block_x_old].erase(i);
        omp_unset_lock(&block_locks[old_block_index]);
        //lock, update, and release the new block
        omp_set_lock(&block_locks[new_block_index]);
        grid[which_block_y][which_block_x].insert(i);
        omp_unset_lock(&block_locks[new_block_index]);
      }
}

//Called in simulate_particles
void compute_force_grid(particle_t* particles, int& navg, double& davg, double& dmin){

  #pragma omp parallel for reduction (+:navg) reduction(+:davg)
  for (int k = 0; k < NUM_BLOCKS_PER_DIM * NUM_BLOCKS_PER_DIM; k++){
      int i = k / NUM_BLOCKS_PER_DIM;
      int j = k % NUM_BLOCKS_PER_DIM;
      
      int navg_ = 0;
      double davg_ = 0.0;

      //set acceleration to zero
      for (set<int>::iterator it = grid[i][j].begin(); it != grid[i][j].end(); it++){
        particles[*it].ax = particles[*it].ay = 0;
      }

      //check right
      if (j != NUM_BLOCKS_PER_DIM - 1){
        compute_force_between_blocks(grid[i][j], grid[i][j+1], particles, navg_, davg_, dmin);
      }
      //check diagonal right bot
      if (j != NUM_BLOCKS_PER_DIM - 1 && i != NUM_BLOCKS_PER_DIM - 1){
        compute_force_between_blocks(grid[i][j], grid[i+1][j+1], particles, navg_, davg_, dmin);
      }
      //check diagonal right top
      if (j != NUM_BLOCKS_PER_DIM - 1 && i != 0){
        compute_force_between_blocks(grid[i][j], grid[i-1][j+1], particles, navg_, davg_, dmin);
      }
      //check left
      if (j != 0){
        compute_force_between_blocks(grid[i][j], grid[i][j-1], particles, navg_, davg_, dmin);
      }
      //check diagonal left bot
      if (j != 0 && i != NUM_BLOCKS_PER_DIM - 1){
        compute_force_between_blocks(grid[i][j], grid[i+1][j-1], particles, navg_, davg_, dmin);
      }
      //check diagonal left top
      if (j != 0 && i != 0){
        compute_force_between_blocks(grid[i][j], grid[i-1][j-1], particles, navg_, davg_, dmin);
      }
      //check top
      if (i != 0){
        compute_force_between_blocks(grid[i][j], grid[i-1][j], particles, navg_, davg_, dmin);
      }
      //check bot
      if (i != NUM_BLOCKS_PER_DIM - 1){
        compute_force_between_blocks(grid[i][j], grid[i+1][j], particles, navg_, davg_, dmin);
      }
      //compute within itself
      compute_force_within_block(grid[i][j], particles, navg_, davg_, dmin);
      navg += navg_;
      davg += davg_;
    }
}

//Call in simulate_particles, after finsih computing force
void move_particles(particle_t* particles, int n){
  #pragma omp parallel
  {
    #pragma omp for
    for (int i = 0; i < n; i++){
      double old_x = particles[i].x;
      double old_y = particles[i].y;
      move(particles[i]);
      //check if the particle might move to another block
      move_to_another_block(i, old_x, old_y, particles);
    }
  }
}

//Called once to generate the grid
void generateGrid(particle_t* particles, int n){
    //initialize the grid with NUM_BLOCKS_PER_DIM^2
    grid = vector<vector<set<int> > >(NUM_BLOCKS_PER_DIM, vector< set<int> >(NUM_BLOCKS_PER_DIM, set<int>()));
    
    //store the point into the grid
    #pragma omp parallel
    {
      #pragma omp for
      for (int i = 0; i < n; i++){
          int which_block_x = min((int)(particles[i].x / BLOCK_SIZE), NUM_BLOCKS_PER_DIM - 1);
          int which_block_y = min((int)(particles[i].y / BLOCK_SIZE), NUM_BLOCKS_PER_DIM - 1);

          //set the lock, update it and release
          int block_index = FIND_POS(which_block_y, which_block_x, NUM_BLOCKS_PER_DIM);
          omp_set_lock(&block_locks[block_index]);
          grid[which_block_y][which_block_x].insert(i);
          omp_unset_lock(&block_locks[block_index]);
      }
    }
}

//Called once to initialize block_lockss
void initializeBlockLocks(){
  //initialize the block_locks with the number of blocks
  block_locks = vector<omp_lock_t>(NUM_BLOCKS_PER_DIM * NUM_BLOCKS_PER_DIM);
  #pragma omp parallel
  {
    #pragma omp for
    for (int i = 0; i < block_locks.size(); i++){
      omp_init_lock(&block_locks[i]);
    }
  }
}

//Called once, no need to parallized
double findSize(particle_t* particles, int n){
  double min_x = 1 << 30;
  double min_y = 1 << 30;
  double max_x = -1;
  double max_y = -1;
  #pragma omp parallel for reduction(max : max_x) reduction(max : max_y) reduction(min : min_x) reduction(min : min_y)
  for (int i = 0; i < n; i++){
      min_x = min(particles[i].x, min_x);
      max_x = max(particles[i].x, max_x);
      min_y = min(particles[i].y, min_y);
      max_y = max(particles[i].y, max_y);
  }
  double size = max(max_x - min_x, max_y - min_y);
  return size;
}

void simulate_particles(particle_t* particles, int n, Serial_output& output, FILE* fsave, int argc, char** argv){

    int navg = 0;
    double davg = 0.0;
    double dmin = 1.0;
    double cutoff = 0.01;

    //initialize a bunck of stuff here

    GRID_SIZE = findSize(particles, n);
    NUM_BLOCKS_PER_DIM = int(sqrt(ceil(n/64.0)*64)) ;
    BLOCK_SIZE = GRID_SIZE / NUM_BLOCKS_PER_DIM;
    if (BLOCK_SIZE < 0.01){
      NUM_BLOCKS_PER_DIM = int((GRID_SIZE /  cutoff));
      BLOCK_SIZE = GRID_SIZE / NUM_BLOCKS_PER_DIM;
    }
    //initialize locks
    initializeBlockLocks();
    //gerenate grid
    generateGrid(particles, n);
   
    
    //debug
    for(int step = 0; step < NSTEPS; step++ ) {
        navg = 0;
        davg = 0.0;
        dmin = 1.0;

        //double computation_begin = read_timer();
        compute_force_grid(particles, navg, davg, dmin);

        //move the particles
        move_particles(particles, n);

        //book keeping
        if( find_option( argc, argv, "-no" ) == -1 ) {

          if (navg) {
            output.absavg +=  davg/navg;
            output.nabsavg++;
          }

          if (dmin < output.absmin){
            output.absmin = dmin;
          }

          if( fsave && (step%SAVEFREQ) == 0 ){
              save( fsave, n, particles );
          }
        }
    }
}

int main( int argc, char **argv )
{
    cout << "Optimized Openmp " << endl;
    int num_threads = -1;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    #pragma omp parallel
    {
      num_threads = omp_get_num_threads();
    }

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );
    //simulate a number of time steps
    double simulation_time = read_timer( );

    //call the simulator
    Serial_output output;
    simulate_particles(particles, n, output, fsave, argc, argv);

    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d,threads = %d, simulation time = %g seconds", n, num_threads, simulation_time);


    if( find_option( argc, argv, "-no" ) == -1 )
    {
        if (output.nabsavg){
            output.absavg /= output.nabsavg;
        }
        printf( ", absmin = %lf, absavg = %lf", output.absmin, output.absavg);
        if (output.absmin < 0.4)
            printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
        if (output.absavg < 0.8) printf
            ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");

    if( fsum)
        fprintf(fsum,"%d %d %g\n",n,num_threads,simulation_time);

    if( fsum )
        fclose( fsum );
    free( particles );
    if( fsave )
        fclose( fsave );

    return 0;
}
