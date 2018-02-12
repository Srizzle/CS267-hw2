#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <vector>
#include <set>
#include <iostream>
#include <pthread.h>

using namespace std;

double GRID_SIZE = 0;
double TOTAL_SIZE = 1;

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

void compute_force_within_block(set<int>& block, particle_t* particles,
          int* navg, double* davg, double* dmin){

    for (set<int>::iterator it_1 = block.begin(); it_1 != block.end(); it_1++){
      for (set<int>::iterator it_2 = block.begin(); it_2 != block.end(); it_2++){
        apply_force(particles[*it_1], particles[*it_2], dmin, davg, navg);
      }
    }
}

void compute_force_between_blocks(set<int>& block_A, set<int>& block_B, particle_t* particles,
          int* navg, double* davg, double* dmin){
  for (set<int>::iterator it_A = block_A.begin(); it_A != block_A.end(); it_A++){
    for (set<int>::iterator it_B = block_B.begin(); it_B != block_B.end(); it_B++){
      apply_force(particles[*it_A], particles[*it_B], dmin, davg, navg);
    }
  }
}

void move_to_another_block(int i, double old_x, double old_y,
                vector<vector<set<int> > >& grid,
                particle_t* particles, double GRID_SIZE){

      int numRows = TOTAL_SIZE / GRID_SIZE;
      int numCols = TOTAL_SIZE / GRID_SIZE;

      int which_block_x_old = min((int)(old_x / GRID_SIZE), numRows - 1);
      int which_block_y_old = min((int)(old_y / GRID_SIZE), numCols - 1);

      int which_block_x = min((int)(particles[i].x / GRID_SIZE), numRows - 1);
      int which_block_y = min((int)(particles[i].y / GRID_SIZE), numCols - 1);

      // cout << which_block_x_old << "," << which_block_y_old << "=>" << which_block_x << "," << which_block_y << endl;

      if (which_block_x_old != which_block_x || which_block_y_old != which_block_y){
        grid[which_block_y_old][which_block_x_old].erase(i);
        grid[which_block_y][which_block_x].insert(i);
      }
}

void compute_force_grid(vector<vector<set<int> > >& grid,
            particle_t* particles, double GRID_SIZE, int* navg, double* davg, double* dmin){
    int numRows = TOTAL_SIZE / GRID_SIZE;
    int numCols = TOTAL_SIZE / GRID_SIZE;
    for (int i = 0; i < numRows; i++){
      for (int j = 0; j < numCols; j++){
        //printf("Block %d, %d has %d particles \n", i, j, grid[i][j].size());

        //set acceleration to zero
        for (set<int>::iterator it = grid[i][j].begin(); it != grid[i][j].end(); it++){
          particles[*it].ax = particles[*it].ay = 0;
        }

        //check right
        if (j != numCols - 1){
          compute_force_between_blocks(grid[i][j], grid[i][j+1], particles, navg, davg, dmin);
        }
        //check diagonal right bot
        if (j != numCols - 1 && i != numRows - 1){
          compute_force_between_blocks(grid[i][j], grid[i+1][j+1], particles, navg, davg, dmin);
        }
        //check diagonal right top
        if (j != numCols - 1 && i != 0){
          compute_force_between_blocks(grid[i][j], grid[i-1][j+1], particles, navg, davg, dmin);
        }
        //check left
        if (j != 0){
          compute_force_between_blocks(grid[i][j], grid[i][j-1], particles, navg, davg, dmin);
        }
        //check diagonal left bot
        if (j != 0 && i != numRows - 1){
          compute_force_between_blocks(grid[i][j], grid[i+1][j-1], particles, navg, davg, dmin);
        }
        //check diagonal left top
        if (j != 0 && i != 0){
          compute_force_between_blocks(grid[i][j], grid[i-1][j-1], particles, navg, davg, dmin);
        }
        //check top
        if (i != 0){
          compute_force_between_blocks(grid[i][j], grid[i-1][j], particles, navg, davg, dmin);
        }
        //check bot
        if (i != numRows - 1){
          compute_force_between_blocks(grid[i][j], grid[i+1][j], particles, navg, davg, dmin);
        }
        //compute within itself
        compute_force_within_block(grid[i][j], particles, navg, davg, dmin);
      }
    }
}

void move_particles(particle_t* particles, int n, vector<vector<set<int> > >& grid, double grid_size){
  for (int i = 0; i < n; i++){
    double old_x = particles[i].x;
    double old_y = particles[i].y;
    move(particles[i]);
    //check if the particle might move to another block
    move_to_another_block(i, old_x, old_y, grid, particles, grid_size);
  }
}


vector<vector<set<int> > > generateGrid(particle_t* particles, int n, double GRID_SIZE){
    //initialize the grid with 10 * 10;
    int numRols = TOTAL_SIZE / GRID_SIZE;
    int numCols = TOTAL_SIZE / GRID_SIZE;
    vector<vector<set<int> > > grid(numRols, vector< set<int> >(numCols, set<int>() ) );

    //store the point into the grid
    for (int i = 0; i < n; i++){
        int which_block_x = (int)(particles[i].x / GRID_SIZE);
        int which_block_y = (int)(particles[i].y / GRID_SIZE);
        grid[min(which_block_y, numRols - 1)][min(which_block_x, numCols - 1)].insert(i);
        //cout << which_block_x << " " << which_block_y << endl;
    }
    for (int i = 0; i < numRols; i++){
      for (int j = 0; j < numCols; j++){
        //cout << "Block " << i << "," << j << " contains " << grid[i][j].size() << " particles" << endl;
      }
    }

    return grid;
}

double findSize(particle_t* particles, int n){
  double min_x = 1 << 30;
  double min_y = 1 << 30;
  double max_x = -1;
  double max_y = -1;
  for (int i = 0; i < n; i++){
      min_x = min(particles[i].x, min_x);
      max_x = max(particles[i].x, max_x);
      min_y = min(particles[i].y, min_y);
      max_y = max(particles[i].y, max_y);
  }
  // cout << min_x << " " << max_x << endl;
  // cout << min_y << " " << max_y << endl;
  double size = max(max_x - min_x, max_y - min_y);
  return size;
}

void simulate_particles(particle_t* particles, int n, Serial_output& output, FILE* fsave, int argc, char** argv){

    int navg = 0;
    double davg = 0.0;
    double dmin = 1.0;
    double cutoff = 0.01;

    TOTAL_SIZE = findSize(particles, n);
    GRID_SIZE = TOTAL_SIZE / ((int)sqrt(n));
    if (GRID_SIZE < 0.01)
      GRID_SIZE = (TOTAL_SIZE / ((int)(ceil(TOTAL_SIZE / cutoff))));
    //cout << TOTAL_SIZE << " " << GRID_SIZE << endl;

    //gerenate grid
    vector<vector<set<int> > > grid = generateGrid(particles, n, GRID_SIZE);

    for( int step = 0; step < NSTEPS; step++ ) {
        navg = 0;
        davg = 0.0;
        dmin = 1.0;

        //compute force within the grid
        compute_force_grid(grid, particles, GRID_SIZE, &navg, &davg, &dmin);

        //move the particles
        move_particles(particles, n, grid, GRID_SIZE);

        //book keeping
        if( find_option( argc, argv, "-no" ) == -1 )
        {
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
        //cout << "-------" << endl;
    }
}


void simulate_particles_naives(particle_t* particles, int n, Serial_output& output, FILE* fsave, int argc, char** argv){

    int navg = 0;
    double davg = 0.0;
    double dmin = 1.0;

    for( int step = 0; step < NSTEPS; step++ ) {
      navg = 0;
        davg = 0.0;
      dmin = 1.0;

        //compute forces
        for( int i = 0; i < n; i++ ){
            particles[i].ax = particles[i].ay = 0;
            for (int j = 0; j < n; j++ )
        apply_force( particles[i], particles[j],&dmin,&davg,&navg);
        }

        //move particles
        for( int i = 0; i < n; i++ )
            move( particles[i] );

        //book keeping
        if( find_option( argc, argv, "-no" ) == -1 )
        {
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

    //int navg,nabsavg=0;
    //double absmin=1.0, absavg=0.0;

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

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );

    //  simulate a number of time steps
    double simulation_time = read_timer( );

    //call the simulator
    Serial_output output;
    //simulate_particles_naives(particles, n, output, fsave, argc, argv);
    simulate_particles(particles, n, output, fsave, argc, argv);

    simulation_time = read_timer( ) - simulation_time;

    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
        if (output.nabsavg){
            output.absavg /= output.nabsavg;
        }
        //
        //  -The minimum distance absmin between 2 particles during the run of the simulation
        //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
        //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
        //
        //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
        //
        printf( ", absmin = %lf, absavg = %lf", output.absmin, output.absavg);
        if (output.absmin < 0.4)
            printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
        if (output.absavg < 0.8) printf
            ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");

    //
    // Printing summary data
    //
    if( fsum)
        fprintf(fsum,"%d %g\n",n,simulation_time);

    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );
    free( particles );
    if( fsave )
        fclose( fsave );

    return 0;
}
