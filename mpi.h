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
#include <cassert>

using namespace std;

#define BLOCKS_INITIALIZATION_TAG 1
#define CLUSTERINFO_INITIALIZATION_TAG 2
#define REQUEST_AND_FEED_EDGES_TAG 3
#define TRANSFER_PARTICLE_TAG 4
#define END_OF_TRANSMISSION_TAG 5

#define MASTER 0
#define MAX_RECV_BUFFER_SIZE 100 //50 particle data type
#define DEBUG 0

extern double size;

//define some datatype
MPI_Datatype PARTICLE;
MPI_Datatype CLUSTERINFO;
MPI_Datatype METADATA;

/*--------------------define some structure we need--------------------*/

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

/*--------------------define some global variables we need--------------------*/
double GRID_SIZE = -1;
double BLOCK_SIZE = -1;
int NUM_BLOCKS_PER_DIM = -1;
double CUT_OFF = 0.01;
int NUM_PARTICLES = -1;
int NUM_PROC = -1;
int RANK = -1;
int NUM_ACTIVE_PROC = -1;
ClusterInfo myClusterInfo;

particle_t TERMINATE_SYMBOL;

vector<vector<Block> > myBlocks;
vector<Block> topEdge;
vector<Block> botEdge;
vector<ClusterInfo> cluster_layout;


/*--------------------define some helpers we need--------------------*/
bool withinRange(int rank, int x, int y){
  ClusterInfo myInfo = cluster_layout[rank];
  if (y < myInfo.start_row || y > myInfo.end_row)
    return false;
  if (x < myInfo.start_col || x > myInfo.end_col)
    return false;
  return true;
}

int locateRecipient(int x, int y, int currentRank){
    int left = currentRank, right = currentRank;
    while (left >= 0 || right < NUM_PROC){
        if (left >= 0){
            if (withinRange(left, x, y))
                return left;
        }
        if (right < NUM_PROC){
            if (withinRange(right, x, y))
                return right;
        }
        left--;
        right++;
    }
    return -1;
}

bool isValidCluster(int rank){
    ClusterInfo thisCluster = cluster_layout[rank];
    if (thisCluster.start_row <= thisCluster.end_row && thisCluster.end_row >=0)
        return true;
}

bool isFirstCluster(int rank){
    return (cluster_layout[rank].start_row == 0) && isValidCluster(rank);
}

bool isLastCluster(int rank){
    return (cluster_layout[rank].end_row == NUM_BLOCKS_PER_DIM - 1) && isValidCluster(rank);
}


void printSummary(){
  printf("Processor %d: GRID_SIZE %f, BLOCK_SIZE %f, NUM_BLOCKS_PER_DIM %d, NUM_PARTICLES %d \n",
        RANK, GRID_SIZE, BLOCK_SIZE, NUM_BLOCKS_PER_DIM, NUM_PARTICLES);
  printf("Processor %d: ClusterInfo; start_row %d, end_row %d, start_col %d, end_col %d \n",
        RANK, cluster_layout[RANK].start_row, cluster_layout[RANK].end_row,
        cluster_layout[RANK].start_col, cluster_layout[RANK].end_col);
}

void printBlocks(){
  ClusterInfo myInfo = cluster_layout[RANK];
  int count = 0;
  for (int i = myInfo.start_row; i <= myInfo.end_row; i++){
    for (int j = myInfo.start_col; j <= myInfo.end_col; j++){
      if (DEBUG == 3)
        printf("%d ", myBlocks[i-myInfo.start_row][j].particles.size());
      count += myBlocks[i-myInfo.start_row][j].particles.size();
    }
    if (DEBUG == 3)
      printf("\n");
  }
  printf("Processor %d: Total %d particles------------------------------- \n", RANK, count);
}

void printParticles(particle_t* particles, int n){
  for (int i = 0; i < n; i++){
    cout << particles[i].x << " " << particles[i].y << endl;
  }
}

void printParticle(particle_t particle){
  cout << particle.x << " " << particle.y << endl;
}
