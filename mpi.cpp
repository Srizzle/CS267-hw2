#include "mpi.h"
using namespace std;

/*-----------------------------------------------Master Method--------------------------------------------*/
//figure out the overall span
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
  double size = max(max_x - min_x, max_y - min_y);
  return size;
}

//Called once to generate the grid at the really beginning
vector<vector<Block> > initializeGrid(particle_t* particles, int n){
    //initialize the grid with NUM_BLOCKS_PER_DIM^2
    vector<vector<Block> > grid = vector<vector<Block> >(NUM_BLOCKS_PER_DIM, vector<Block>(NUM_BLOCKS_PER_DIM, Block()));

    //store the point into the grid
    for (int i = 0; i < n; i++){
        int which_block_x = min((int)(particles[i].x / BLOCK_SIZE), NUM_BLOCKS_PER_DIM - 1);
        int which_block_y = min((int)(particles[i].y / BLOCK_SIZE), NUM_BLOCKS_PER_DIM - 1);
        grid[which_block_y][which_block_x].particles.push_back(particles[i]);
    }
    return grid;
}

//Called once to initialize which processor is reponsible for which blocks/rows of grid
vector<ClusterInfo> initializeClusterInfos(int NUM_PROC){
  vector<ClusterInfo> clusterInfo = vector<ClusterInfo>();
  int row_stride = ceil(NUM_BLOCKS_PER_DIM / double(NUM_PROC));
  for (int start_row = 0; start_row < NUM_BLOCKS_PER_DIM; start_row += row_stride){
    int end_row = min(start_row + row_stride - 1, NUM_BLOCKS_PER_DIM - 1);
    int start_col = 0;
    int end_col = NUM_BLOCKS_PER_DIM - 1;
    clusterInfo.push_back(ClusterInfo(start_row, end_row, start_col, end_col));
  }
  for (int i = clusterInfo.size(); i < NUM_PROC; i++){
    clusterInfo.push_back(ClusterInfo(-1,-2,-1,-2));
  }
  return clusterInfo;
}

//Called by the master, to dispense the blocks to the processors at really beginning
void dispense_blocks(vector<ClusterInfo> clusterInfo, vector<vector<Block> > grid){
  //MPI_Request req[NUM_BLOCKS_PER_D-IM *  NUM_BLOCKS_PER_DIM];
  int count = 0;
  //Directly put the data to master's own block
  for (int i = clusterInfo[MASTER].start_row; i <= clusterInfo[MASTER].end_row; i++){
    vector<Block> line;
    for (int j = clusterInfo[MASTER].start_col; j <= clusterInfo[MASTER].end_col; j++){
      line.push_back(grid[i][j]);
    }
    myBlocks.push_back(line);
  }

  //iterate through the ClusterInfo
  for (int i = 1; i < clusterInfo.size(); i++){
    //it is not initialized
    if (clusterInfo[i].start_row == -1)
      continue;
    //send blocks to the correpsonding processor
    for (int which_row = clusterInfo[i].start_row; which_row <= clusterInfo[i].end_row; which_row++){
      for (int which_col = clusterInfo[i].start_col; which_col <= clusterInfo[i].end_col; which_col++){
        MPI_Send(&grid[which_row][which_col].particles.front(), grid[which_row][which_col].particles.size(), PARTICLE,
          i, BLOCKS_INITIALIZATION_TAG, MPI_COMM_WORLD);
      }
    }
    // printf("Master send %d particles in Block [%d, %d] to %d process \n",
    //        grid[clusterInfo[i].end_row][clusterInfo[i].end_col].particles.size(), clusterInfo[i].end_row, clusterInfo[i].end_col, i);
  }
  //MPI_Waitall(count, req, MPI_STATUSES_IGNORE);
}

//Called by the master, to broadcast the metaData
void dispense_meta_data(double GRID_SIZE, double BLOCK_SIZE, int NUM_BLOCKS_PER_DIM, int n){
  MetaData metaData(GRID_SIZE, BLOCK_SIZE, NUM_BLOCKS_PER_DIM, n);
  MPI_Bcast(&metaData, 1, METADATA, 0, MPI_COMM_WORLD);
}

//Called by the master, to dispense the clusterInfo to each process
void dispense_clusterInfo(vector<ClusterInfo> clusterInfo){
  MPI_Bcast(&clusterInfo.front(), clusterInfo.size(), CLUSTERINFO, 0, MPI_COMM_WORLD);
}

//Master routine
void master_routine(particle_t* particles, int n){
  //initialize particles
  init_particles( n, particles );
  //find the size
  GRID_SIZE = findSize(particles, n);
  //generate the important constants
  NUM_BLOCKS_PER_DIM = int(sqrt(ceil(n/64.0)*64)) ;
  BLOCK_SIZE = GRID_SIZE / NUM_BLOCKS_PER_DIM;
  NUM_PARTICLES = n;
  if (BLOCK_SIZE < 0.01){
    NUM_BLOCKS_PER_DIM = int(GRID_SIZE /  CUT_OFF);
    BLOCK_SIZE = GRID_SIZE / NUM_BLOCKS_PER_DIM;
  }
  //broadcast metadata
  dispense_meta_data(GRID_SIZE, BLOCK_SIZE, NUM_BLOCKS_PER_DIM, n);
  //initialize the grid
  vector<vector<Block> > grid = initializeGrid(particles, n);
  //initialize the ClusterInfos info (which processor is reponsible for which rows)
  cluster_layout = initializeClusterInfos(NUM_PROC);
  //dispense clusterInfo to all processes
  dispense_clusterInfo(cluster_layout);
  //send the blocks to correpsonding the processors
  dispense_blocks(cluster_layout, grid);
}

/*-----------------------------------------------Worker Method---------------------------------------------------------*/

//Called by all the processors to receive clusterInfo from the master
void receive_clusterInfo_from_master(int source){
  cluster_layout = vector<ClusterInfo>(NUM_PROC, ClusterInfo());
  MPI_Bcast(&cluster_layout.front(), cluster_layout.size(), CLUSTERINFO, source, MPI_COMM_WORLD);

  //initialize the top and bot edge buffer
  ClusterInfo myInfo = cluster_layout[RANK];
  // topEdge = vector<Block>(myInfo.end_col - myInfo.start_col + 1);
  // botEdge = vector<Block>(myInfo.end_col - myInfo.start_col + 1);
}

//Called by all the processors to receive metadata broadcasted from the master
void receive_metaData_from_master(int source){
  MetaData metaData;
  MPI_Bcast(&metaData, 1, METADATA, source, MPI_COMM_WORLD);
  GRID_SIZE = metaData.GRID_SIZE;
  BLOCK_SIZE = metaData.BLOCK_SIZE;
  NUM_PARTICLES = metaData.NUM_PARTICLES;
  NUM_BLOCKS_PER_DIM = metaData.NUM_BLOCKS_PER_DIM;
}

//Called by all the processors to receive blocks from the master
void receive_blocks_from_master(int source, int tag){
  int real_num_particles;
  MPI_Status status;
  for (int i = cluster_layout[RANK].start_row; i <= cluster_layout[RANK].end_row; i++){
    vector<Block> line;
    for (int j = cluster_layout[RANK].start_col; j <= cluster_layout[RANK].end_col; j++){
      vector<particle_t> buffer(MAX_RECV_BUFFER_SIZE, particle_t());
      MPI_Recv(&buffer.front(), MAX_RECV_BUFFER_SIZE, PARTICLE, source, tag, MPI_COMM_WORLD, &status);
      MPI_Get_count(&status, PARTICLE, &real_num_particles);
      buffer.resize(real_num_particles);
      line.push_back(Block(buffer));
      // printf("Processor %d: Received %d particles in Block [%d, %d] tag = %d \n",
      //   RANK, real_num_particles, i, j, status.MPI_TAG);
    }
    myBlocks.push_back(line);
  }
  // printf("Process %d: Block [%d, %d] has %d particles \n",
  //   RANK, cluster_layout[RANK].end_row, cluster_layout[RANK].end_col, myBlocks.back().back().particles.size());
}

//request and edge blocks to other processors (currently the one above and below it)
void request_and_feed_edges(int tag){
  ClusterInfo myInfo = cluster_layout[RANK];
  /* ---------------------------------- Send out my edges ----------------------------------*/
  int recepient = -1;
  //send top edge
  if (RANK != 0){
    recepient = RANK - 1;
    for (int col = myInfo.start_col; col <= myInfo.end_col; col++){
      MPI_Request request;
      MPI_Isend(&myBlocks[0][col].particles.front(), myBlocks[0][col].particles.size(),
                PARTICLE, recepient, tag, MPI_COMM_WORLD, &request);
      if (DEBUG == 3)
        printf("Processor %d: Send out block [%d, %d] to Processor %d, with %d particles \n",
          RANK, myInfo.start_row, col, recepient, myBlocks[0][col].particles.size());
    }
  }
  //send bot edge
  if (RANK != NUM_PROC - 1){
    recepient = RANK + 1;
    for (int col = myInfo.start_col; col <= myInfo.end_col; col++){
      MPI_Request request;
      MPI_Isend(&myBlocks[myInfo.end_row - myInfo.start_row][col].particles.front(), myBlocks[myInfo.end_row - myInfo.start_row] [col].particles.size(),
                PARTICLE, recepient, tag, MPI_COMM_WORLD, &request);
      if (DEBUG == 3)
        printf("Processor %d: Send out block [%d, %d] to Processor %d, with %d particles \n",
          RANK, myInfo.end_row, col, recepient, myBlocks[myInfo.end_row - myInfo.start_row][col].particles.size());
    }
  }

  /* ---------------------------------- Receive edges ----------------------------------*/
  int sender = -1;
  int real_num_particles;
  MPI_Status status;
  //receive top edge
  if (RANK != 0){
    sender = RANK - 1;
    for (int col = myInfo.start_col; col <= myInfo.end_col; col++){
      vector<particle_t> buffer(MAX_RECV_BUFFER_SIZE, particle_t());
      MPI_Recv(&buffer.front(), MAX_RECV_BUFFER_SIZE, PARTICLE, sender, tag, MPI_COMM_WORLD, &status);
      MPI_Get_count(&status, PARTICLE, &real_num_particles);
      buffer.resize(real_num_particles);
      topEdge.push_back(Block(buffer));
      if (DEBUG == 3)
        printf("Processor %d: Receive block [%d, %d] from Processor %d, with %d particles \n",
          RANK, myInfo.start_row - 1, col, sender, topEdge[col].particles.size());
    }
  }
  //receive bot edge
  if (RANK != NUM_PROC-1){
    sender = RANK + 1;
    for (int col = myInfo.start_col; col <= myInfo.end_col; col++){
      vector<particle_t> buffer(MAX_RECV_BUFFER_SIZE, particle_t());
      MPI_Recv(&buffer.front(), MAX_RECV_BUFFER_SIZE, PARTICLE, sender, tag, MPI_COMM_WORLD, &status);
      MPI_Get_count(&status, PARTICLE, &real_num_particles);
      buffer.resize(real_num_particles);
      botEdge.push_back(Block(buffer));
      if (DEBUG == 3)
        printf("Processor %d: Receive block [%d, %d] from Processor %d, with %d particles \n",
          RANK, myInfo.end_row + 1, col, sender, botEdge[col].particles.size());
    }
  }
}

//send particle to another block on a separate processor
void transfer_particle(particle_t particle, int recepient, int tag){
  MPI_Request request;
  MPI_Isend(&particle, 1, PARTICLE, recepient, tag, MPI_COMM_WORLD, &request);
}

//decide memmbership;
/*
1. No membership change.
2. Membership changed to another block within the processor.
3. Membership changed to another block outside the processor
*/
void decide_membership(Block& currentBlock, double old_x, double old_y, particle_t particle){
  ClusterInfo myInfo = cluster_layout[RANK];
  int which_block_x_old = min((int)(old_x / BLOCK_SIZE), NUM_BLOCKS_PER_DIM - 1);
  int which_block_y_old = min((int)(old_y / BLOCK_SIZE), NUM_BLOCKS_PER_DIM - 1);
  int which_block_x = min((int)(particle.x / BLOCK_SIZE), NUM_BLOCKS_PER_DIM - 1);
  int which_block_y = min((int)(particle.y / BLOCK_SIZE), NUM_BLOCKS_PER_DIM - 1);

  if (which_block_x_old != which_block_x || which_block_y_old != which_block_y){
    //case 2
    if (withinRange(RANK, which_block_x, which_block_y)){
      myBlocks[which_block_y - myInfo.start_row][which_block_x].particles.push_back(particle);
    //case 3
    }else{
      transfer_particle(particle, locateRecipient(which_block_x, which_block_y), TRANSFER_PARTICLE_TAG);
    }
  }else{
    //case 1
    currentBlock.particles.push_back(particle);
  }
}


//receive particle from another processor, figure out in which block to put it.
void poll_particles(int tag){
  ClusterInfo myInfo = cluster_layout[RANK];
  MPI_Status status;
  particle_t placeholder;
  int finished_processes = 0; //indicates how many other processes havve finished sending their particles to this process
  //we wont stop until receive terminate symbols from all  processes (including itself)
  while (finished_processes != NUM_PROC){
    //polling any incoming particles
    MPI_Recv(&placeholder, 1, PARTICLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
    //check if it is a terminate symbol //x == y == -1 is our termiante symbols
    if (placeholder.x == -1 && placeholder.y == -1){
      finished_processes++;
      if (DEBUG == 2)
        printf("Processor: %d: Finish polling from Processor %d \n", status.MPI_SOURCE);
      continue;
    }
    //otherwise, first we need to compute which block it belongs to
    int which_block_x = min((int)(placeholder.x / BLOCK_SIZE), NUM_BLOCKS_PER_DIM - 1);
    int which_block_y = min((int)(placeholder.y / BLOCK_SIZE), NUM_BLOCKS_PER_DIM - 1);
    //assert here
    assert((myInfo.start_row <= which_block_y) && (myInfo.end_row >= which_block_y));
    assert((myInfo.start_col <= which_block_x) && (myInfo.end_col >= which_block_x));
    //insert
    myBlocks[which_block_y - myInfo.start_col][which_block_x].particles.push_back(placeholder);
  }

  if (DEBUG == 2){
    printf("Processor %d: Finish polling particles from all other processes \n", RANK);
  }
}

//move the particles to appropriate blocks (can be remote) after computing the force
void move_particles(){
  ClusterInfo myInfo = cluster_layout[RANK];
  vector<vector<Block> > oldBlocks = myBlocks;
  //clear everything
  for (int i = 0; i <= myInfo.end_row - myInfo.start_row; i++){
    for (int j = myInfo.start_col; j <= myInfo.end_col; j++){
      myBlocks[i][j].particles.clear();
    }
  }
  //loop through the old block
  for (int i = 0; i <= myInfo.end_row - myInfo.start_row; i++){
    for (int j = myInfo.start_col; j <= myInfo.end_col; j++){
      Block oldBlock = oldBlocks[i][j];
      for (int k = 0; k < oldBlock.particles.size(); k++){
        double old_x = oldBlock.particles[k].x;
        double old_y = oldBlock.particles[k].y;
        //update position
        move(oldBlock.particles[k]);
        //check if the particle might move to another block
        decide_membership(myBlocks[i][j], old_x, old_y, oldBlock.particles[k]);
      }
    }
  }
  //let all the processes excluding itself knows that it is done
  MPI_Request request;
  particle_t termiate_symbol;
  termiate_symbol.x = -1;
  termiate_symbol.y = -1;
  for (int proc = 0; proc < NUM_PROC; proc++){
    MPI_Isend(&termiate_symbol, 1, PARTICLE, proc, TRANSFER_PARTICLE_TAG, MPI_COMM_WORLD, &request);
  }
  //polling for the particles we need
  poll_particles(TRANSFER_PARTICLE_TAG);

  if (DEBUG == 2)
    printf("Process %d: Finish moving particles \n", RANK);
}





int main( int argc, char **argv )
{
    int navg, nabsavg=0;
    double dmin, absmin=1.0,davg,absavg=0.0;
    double rdavg,rdmin;
    int rnavg;

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

    //set up MPI
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &NUM_PROC );
    MPI_Comm_rank( MPI_COMM_WORLD, &RANK );
    //define MPI Particle
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );
    //define MPI ClusterInfo
    MPI_Type_contiguous( 4, MPI_INT, &CLUSTERINFO );
    MPI_Type_commit( &CLUSTERINFO );
    //define MPI MetaData
    MPI_Type_contiguous( 4, MPI_DOUBLE, &METADATA );
    MPI_Type_commit( &METADATA );

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );

    if(RANK == MASTER){
      master_routine(particles, n);
    }else{
      receive_metaData_from_master(MASTER);
      receive_clusterInfo_from_master(MASTER);
      receive_blocks_from_master(MASTER, BLOCKS_INITIALIZATION_TAG);
    }

    //debug
    if (DEBUG){
      printSummary();
      printBlocks();
    }

    request_and_feed_edges(REQUEST_AND_FEED_EDGES_TAG);
    move_particles();

    if (DEBUG){
      printBlocks();
    }

    free( particles );

    MPI_Finalize( );

    return 0;
}
