
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <vector>
#include <set>
#include <iostream>

using namespace std;

#define FIND_POS(ROW_INDEX, COL_INDEX, NUM_BLOCKS_PER_DIM) (ROW_INDEX * NUM_BLOCKS_PER_DIM + COL_INDEX)

//the length of each blocks
double BLOCK_SIZE = 0;
//the maximum span of the particles in both dimension
double GRID_SIZE = 1;
//number of blocks per dim, so the overall number of number of blocks is its square
int NUM_BLOCKS_PER_DIM = -1;


vector<vector<set<int> > > grid;

//Called in compute_force_grid
void compute_force_within_block(set<int>& block, particle_t* particles){

    for (set<int>::iterator it_1 = block.begin(); it_1 != block.end(); it_1++){
      for (set<int>::iterator it_2 = block.begin(); it_2 != block.end(); it_2++){
        apply_force(particles[*it_1], particles[*it_2]);
      }
    }
}

//Called in compute_force_grid
void compute_force_between_blocks(set<int>& block_A, set<int>& block_B, particle_t* particles){
  for (set<int>::iterator it_A = block_A.begin(); it_A != block_A.end(); it_A++){
    for (set<int>::iterator it_B = block_B.begin(); it_B != block_B.end(); it_B++){
      apply_force(particles[*it_A], particles[*it_B]);
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
        grid[which_block_y_old][which_block_x_old].erase(i);
        grid[which_block_y][which_block_x].insert(i);
        // cout << "Membership change" << endl;
        // printf("From %d %d to %d %d \n", which_block_y_old, which_block_x_old, which_block_y, which_block_x);
      }
}

//Called in simulate_particles
void compute_force_grid(particle_t* particles){

  for (int k = 0; k < NUM_BLOCKS_PER_DIM * NUM_BLOCKS_PER_DIM; k++){
      int i = k / NUM_BLOCKS_PER_DIM;
      int j = k % NUM_BLOCKS_PER_DIM;

      //set acceleration to zero
      for (set<int>::iterator it = grid[i][j].begin(); it != grid[i][j].end(); it++){
        particles[*it].ax = particles[*it].ay = 0;
      }

      //check right
      if (j != NUM_BLOCKS_PER_DIM - 1){
        compute_force_between_blocks(grid[i][j], grid[i][j+1], particles);
      }
      //check diagonal right bot
      if (j != NUM_BLOCKS_PER_DIM - 1 && i != NUM_BLOCKS_PER_DIM - 1){
        compute_force_between_blocks(grid[i][j], grid[i+1][j+1], particles);
      }
      //check diagonal right top
      if (j != NUM_BLOCKS_PER_DIM - 1 && i != 0){
        compute_force_between_blocks(grid[i][j], grid[i-1][j+1], particles);
      }
      //check left
      if (j != 0){
        compute_force_between_blocks(grid[i][j], grid[i][j-1], particles);
      }
      //check diagonal left bot
      if (j != 0 && i != NUM_BLOCKS_PER_DIM - 1){
        compute_force_between_blocks(grid[i][j], grid[i+1][j-1], particles);
      }
      //check diagonal left top
      if (j != 0 && i != 0){
        compute_force_between_blocks(grid[i][j], grid[i-1][j-1], particles);
      }
      //check top
      if (i != 0){
        compute_force_between_blocks(grid[i][j], grid[i-1][j], particles);
      }
      //check bot
      if (i != NUM_BLOCKS_PER_DIM - 1){
        compute_force_between_blocks(grid[i][j], grid[i+1][j], particles);
      }
      //compute within itself
      compute_force_within_block(grid[i][j], particles);
    }
}

//Call in simulate_particles, after finsih computing force
void move_particles(particle_t* particles, int n){

  for (int i = 0; i < n; i++){
    double old_x = particles[i].x;
    double old_y = particles[i].y;
    move(particles[i]);
    //check if the particle might move to another block
    move_to_another_block(i, old_x, old_y, particles);
  }
}


//Called once to generate the grid
void generateGrid(particle_t* particles, int n){
  //initialize the grid with NUM_BLOCKS_PER_DIM^2
  grid = vector<vector<set<int> > >(NUM_BLOCKS_PER_DIM, vector< set<int> >(NUM_BLOCKS_PER_DIM, set<int>()));
  
  //store the point into the grid
  for (int i = 0; i < n; i++){
      int which_block_x = min((int)(particles[i].x / BLOCK_SIZE), NUM_BLOCKS_PER_DIM - 1);
      int which_block_y = min((int)(particles[i].y / BLOCK_SIZE), NUM_BLOCKS_PER_DIM - 1);
      grid[which_block_y][which_block_x].insert(i);
  }
}

//Called once, no need to parallized
double findSize(){
   extern double size;
   return size;
}

void simulate_particles(particle_t* particles, int n, FILE* fsave, int argc, char** argv){

    //initialize a bunck of stuff here
    GRID_SIZE = findSize();
    NUM_BLOCKS_PER_DIM = int(sqrt(ceil(n/64.0)*64)) ;
    BLOCK_SIZE = GRID_SIZE / NUM_BLOCKS_PER_DIM;
    if (BLOCK_SIZE < 0.01){
      NUM_BLOCKS_PER_DIM = int((GRID_SIZE /  cutoff));
      BLOCK_SIZE = GRID_SIZE / NUM_BLOCKS_PER_DIM;
    }
    //gerenate grid
    generateGrid(particles, n);

    for(int step = 0; step < 50; step++ ) {

        //double computation_begin = read_timer();
        compute_force_grid(particles);

        //move the particles
        move_particles(particles, n);

        //save if necessary
        if( fsave && (step%1) == 0 ) {
          save( fsave, n, particles);
        }
    }
}

int main( int argc, char **argv ) {
  
    if( find_option( argc, argv, "-h" ) >= 0 ) {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );
    //simulate a number of time steps
    double simulation_time = read_timer( );

    //call the simulator
    simulate_particles(particles, n, fsave, argc, argv);

    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d simulation time = %g seconds \n", n, simulation_time);

    if( fsum)
        fprintf(fsum,"%d %g\n",n,simulation_time);

    if( fsum )
        fclose( fsum );
    free( particles );
    if( fsave )
        fclose( fsave );

    return 0;
}