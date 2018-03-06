#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <vector>
#include <iostream>
#include "common.h"


using namespace std;

#define NUM_THREADS 256

extern double size;

// Based on StackOverflow
#define GPUERRCHK(ans) { gpuAssert((ans), __FILE__, __LINE__); }
#define MAX_BIN_SIZE 10
#define DEBUG 0

__host__ int FIND_POS_HOST(int i, int j, int stride){
  return i * stride + j;
}

__device__ int FIND_POS_DEVICE(int i, int j, int stride){
  return i * stride + j;
}

struct Bin{
  particle_t particles[MAX_BIN_SIZE];
  int ids[MAX_BIN_SIZE];
  int currentSize;
  __device__ Bin(){
    currentSize = 0;
  };
  __host__ Bin(int x){
    currentSize = 0;
  }
  void addParticle(particle_t particle, int id){
    if (currentSize >= MAX_BIN_SIZE){
      printf("HOST WARNING: Not enough room in the bin, current size: %d\n", currentSize);
      return;
    }
    particles[currentSize] = particle;
    ids[currentSize] = id;
    currentSize += 1;
  };
};

//the length of each blocks
double BIN_SIZE = 0;
//the maximum span of the particles in both dimension
double GRID_SIZE = 1;
//number of blocks per dim, so the overall number of number of blocks is its square
int NUM_BINS_PER_DIM = -1;

//How many gpu_blocks to generate for gpus
int NUM_GPU_BLOCKS = -1;


inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

__device__ void addParticle(Bin& bin, particle_t p, int id){
  int old_pos = atomicAdd(&(bin.currentSize), 1);
  if (old_pos >= MAX_BIN_SIZE){
    printf("DEVICE WARNING: Not enough room in the bin, current size: %d\n", old_pos);
    return;
  }
  bin.particles[old_pos] = p;
  bin.ids[old_pos] = id;
};

__device__ void apply_force_gpu(particle_t &particle, particle_t &neighbor) {
  double dx = neighbor.x - particle.x;
  double dy = neighbor.y - particle.y;
  double r2 = dx * dx + dy * dy;
  if( r2 > cutoff*cutoff )
      return;
  r2 = fmax( r2, min_r*min_r );
  //r2 = (r2 > min_r*min_r) ? r2 : min_r*min_r;
  double r = sqrt( r2 );
  double coef = ( 1 - cutoff / r ) / r2 / mass;
  particle.ax += coef * dx;
  particle.ay += coef * dy;

}

// Called in compute_force_grid
__device__ void compute_force_within_block(Bin& bin){
    for (int i = 0; i < bin.currentSize; i++){
      for (int j = 0; j < bin.currentSize; j++){
        apply_force_gpu(bin.particles[i], bin.particles[j]);
      }
    }
}

// Called in compute_force_grid
__device__ void compute_force_between_blocks(Bin& bin_A, Bin& bin_B){
  for (int i = 0; i < bin_A.currentSize; i++){
    for (int j = 0; j < bin_B.currentSize; j++){
      apply_force_gpu(bin_A.particles[i], bin_B.particles[j]);
    }
  }
}

// Called in move_particles
__device__ void move_gpu(particle_t* p, double GRID_SIZE){

  p->vx += p->ax * dt;
  p->vy += p->ay * dt;
  p->x  += p->vx * dt;
  p->y  += p->vy * dt;
  while( p->x < 0 || p->x > GRID_SIZE ) {
      p->x  = p->x < 0 ? -(p->x) : 2*GRID_SIZE-p->x;
      p->vx = -(p->vx);
  }
  while( p->y < 0 || p->y > GRID_SIZE ) {
      p->y  = p->y < 0 ? -(p->y) : 2*GRID_SIZE-p->y;
      p->vy = -(p->vy);
  }
}

// Called in move_particles
// Add the particle to redundantBins, not our current bins
__device__ void bin_change(Bin* redundantBins, particle_t p, int id, double BIN_SIZE, int NUM_BINS_PER_DIM){
  int which_block_x = min((int)(p.x / BIN_SIZE), NUM_BINS_PER_DIM - 1);
  int which_block_y = min((int)(p.y / BIN_SIZE), NUM_BINS_PER_DIM - 1);
  //atomic operation for addParticle
  addParticle(redundantBins[FIND_POS_DEVICE(which_block_y, which_block_x, NUM_BINS_PER_DIM)], p, id);
}
