#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <vector>
#include <set>
#include <iostream>
#include <iterator>

using namespace std;

#define BLOCKS_INITIALIZATION_TAG 1
#define CLUSTERINFO_INITIALIZATION_TAG 2
#define REQUEST_AND_FEED_EDGES 3
#define MASTER 0
#define MAX_RECV_BUFFER_SIZE 100 //1000 particle data type
#define DEBUG 0

//define some datatype
MPI_Datatype PARTICLE;
MPI_Datatype CLUSTERINFO;
MPI_Datatype METADATA;


struct Block{
  vector<particle_t> particles;
  Block(){};
  Block(vector<particle_t> buffer){
    this->particles = buffer;
  }
};

struct ClusterInfo{
  int start_row;
  int end_row;
  int start_col;
  int end_col;
  ClusterInfo(int start_row_, int end_row_, int start_col_, int end_col_){
    start_row = start_row_;
    start_col = start_col_;
    end_row = end_row_;
    end_col = end_col_;
  };
  ClusterInfo(){
    start_row = -1;
    end_row = -2;
    start_col = -1;
    end_col = -2;
  }
};

struct MetaData{
  double NUM_BLOCKS_PER_DIM;
  double NUM_PARTICLES;
  double BLOCK_SIZE;
  double GRID_SIZE;
  MetaData(double grid_size, double block_size, double num_blocks_per_dim, double n){
    BLOCK_SIZE = block_size;
    GRID_SIZE = grid_size;
    NUM_BLOCKS_PER_DIM = num_blocks_per_dim;
    NUM_PARTICLES = n;
  };
  MetaData(){
    NUM_BLOCKS_PER_DIM = -1;
    NUM_PARTICLES = -1;
    BLOCK_SIZE = -1;
    GRID_SIZE = -1;
  };
};
