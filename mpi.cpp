#include "mpi.h"
using namespace std;

/*-----------------------------------------------Master Method--------------------------------------------*/
//figure out the overall span
double findSize(particle_t* particles, int n){
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

  int count = 0;
  int start_row = 0;
  int row_stride = (NUM_BLOCKS_PER_DIM / double(NUM_PROC));
  while (start_row < NUM_BLOCKS_PER_DIM && count < NUM_PROC){
    int end_row = start_row + row_stride - 1;
    if ((NUM_BLOCKS_PER_DIM - start_row) / (double)(NUM_PROC - count) > row_stride){
        end_row++;
    }
    end_row = min(end_row, NUM_BLOCKS_PER_DIM - 1);
    int start_col = 0;
    int end_col = NUM_BLOCKS_PER_DIM - 1;
    clusterInfo.push_back(ClusterInfo(start_row, end_row, start_col, end_col));
    count++;
    start_row = end_row + 1;
  }


  for (int i = clusterInfo.size(); i < NUM_PROC; i++){
    clusterInfo.push_back(ClusterInfo(-1,-2,-1,-2));
  }
  return clusterInfo;
}

//Called by the master, to dispense the blocks to the processors at really beginning
void dispense_blocks(vector<ClusterInfo> clusterInfo, vector<vector<Block> > grid){
  MPI_Request req[NUM_BLOCKS_PER_DIM *  NUM_BLOCKS_PER_DIM];
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
    if (!isValidCluster(i))
      continue;
    //send blocks to the correpsonding processor
    for (int which_row = clusterInfo[i].start_row; which_row <= clusterInfo[i].end_row; which_row++){
      for (int which_col = clusterInfo[i].start_col; which_col <= clusterInfo[i].end_col; which_col++){
        MPI_Isend(&grid[which_row][which_col].particles.front(), grid[which_row][which_col].particles.size(), PARTICLE,
          i, BLOCKS_INITIALIZATION_TAG, MPI_COMM_WORLD, &req[count]);
        count++;
      }
    }
  }
  MPI_Waitall(count, req, MPI_STATUSES_IGNORE);
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
  NUM_BLOCKS_PER_DIM = int(sqrt(ceil(n/64.0)*16)) ;
  BLOCK_SIZE = GRID_SIZE / NUM_BLOCKS_PER_DIM;
  NUM_PARTICLES = n;
  if (BLOCK_SIZE < 0.01){
    NUM_BLOCKS_PER_DIM = max(1,int(GRID_SIZE /  CUT_OFF));
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
  ClusterInfo myCluster = cluster_layout[RANK];
  int real_num_particles;
  MPI_Status status;
  for (int i = 0; i <= myCluster.end_row - myCluster.start_row; i++){
    vector<Block> line;
    for (int j = 0; j <= myCluster.end_col - myCluster.start_col; j++){
      vector<particle_t> buffer(MAX_RECV_BUFFER_SIZE, particle_t());
      MPI_Recv(&buffer.front(), MAX_RECV_BUFFER_SIZE, PARTICLE, source, tag, MPI_COMM_WORLD, &status);
      MPI_Get_count(&status, PARTICLE, &real_num_particles);
      buffer.resize(real_num_particles);
      line.push_back(Block(buffer));
    }
    myBlocks.push_back(line);
  }
}

//request and edge blocks to other processors (currently the one above and below it)
void request_and_feed_edges(int tag){
  ClusterInfo myInfo = cluster_layout[RANK];
  /* ---------------------------------- Send out my edges ----------------------------------*/
  int recepient = -1;
  //send top edge
  if (isValidCluster(RANK) && !isFirstCluster(RANK)){
    recepient = RANK - 1;
    for (int col = 0; col <= myInfo.end_col - myInfo.start_col; col++){
      MPI_Request request;
      MPI_Isend(&myBlocks[0][col].particles.front(), myBlocks[0][col].particles.size(),
                PARTICLE, recepient, tag, MPI_COMM_WORLD, &request);
      if (DEBUG == 3)
        printf("Processor %d: Send out block [%d, %d] to Processor %d, with %d particles \n",
          RANK, myInfo.start_row, col, recepient, myBlocks[0][col].particles.size());
    }
  }
  //send bot edge
  if (isValidCluster(RANK) && !isLastCluster(RANK)){
    recepient = RANK + 1;
    for (int col = 0; col <= myInfo.end_col - myInfo.start_col; col++){
      MPI_Request request;
      MPI_Isend(&myBlocks[myInfo.end_row - myInfo.start_row][col].particles.front(), myBlocks[myInfo.end_row - myInfo.start_row][col].particles.size(),
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
  topEdge.clear();
  botEdge.clear();

  int num_blocks_should_recv = 0;
  if (isValidCluster(RANK) && !isFirstCluster(RANK) && !isLastCluster(RANK)){
    num_blocks_should_recv = 2 * (myInfo.end_col - myInfo.start_col + 1);
  }else if (isFirstCluster(RANK) && !isLastCluster(RANK)){
    num_blocks_should_recv = 1 * (myInfo.end_col - myInfo.start_col + 1);
  }else if (isLastCluster(RANK) && !isFirstCluster(RANK)){
    num_blocks_should_recv = 1 * (myInfo.end_col - myInfo.start_col + 1);
  }

  int topEdge_count = 0;
  int botEdge_count = 0;

  while (topEdge_count + botEdge_count < num_blocks_should_recv){
    vector<particle_t> buffer(MAX_RECV_BUFFER_SIZE, particle_t());
    MPI_Recv(&buffer.front(), MAX_RECV_BUFFER_SIZE, PARTICLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, PARTICLE, &real_num_particles);
    buffer.resize(real_num_particles);

    if (status.MPI_SOURCE == RANK - 1){
        topEdge.push_back(Block(buffer));
        topEdge_count++;
    }else if (status.MPI_SOURCE == RANK + 1){
        botEdge.push_back(Block(buffer));
        botEdge_count++;
    }else{
        printf("WARNING: %d RECEIVE FROM %d during the request_and_feed stage \n", RANK, status.MPI_SOURCE);
    }
  }
}

//send particle to another block on a separate processor
MPI_Request transfer_particle(particle_t& particle, int recepient, int tag){
  MPI_Request request;
  MPI_Isend(&particle, 1, PARTICLE, recepient, tag, MPI_COMM_WORLD, &request);
  return request;
}

//decide memmbership;
/*
1. No membership change.
2. Membership changed to another block within the processor.
3. Membership changed to another block outside the processor
*/
MPI_Request decide_membership(Block& currentBlock, double old_x, double old_y, particle_t& particle){
  ClusterInfo myInfo = cluster_layout[RANK];
  int which_block_x_old = min((int)(old_x / BLOCK_SIZE), NUM_BLOCKS_PER_DIM - 1);
  int which_block_y_old = min((int)(old_y / BLOCK_SIZE), NUM_BLOCKS_PER_DIM - 1);
  int which_block_x = min((int)(particle.x / BLOCK_SIZE), NUM_BLOCKS_PER_DIM - 1);
  int which_block_y = min((int)(particle.y / BLOCK_SIZE), NUM_BLOCKS_PER_DIM - 1);

  if (DEBUG){
    assert((myInfo.start_row <= which_block_y_old) && (myInfo.end_row >= which_block_y_old));
    assert((myInfo.start_col <= which_block_x_old) && (myInfo.end_col >= which_block_x_old));
  }

  MPI_Request request = NULL;
  if (which_block_x_old != which_block_x || which_block_y_old != which_block_y){
    request = transfer_particle(particle, locateRecipient(which_block_x, which_block_y, RANK), TRANSFER_PARTICLE_TAG);
  }else{
    //case 1
    currentBlock.particles.push_back(particle);
  }
  return request;
}


//receive particle from another processor, figure out in which block to put it.
void poll_particles(int tag){
  ClusterInfo myInfo = cluster_layout[RANK];
  MPI_Status status;
  int finished_processes = 0; //indicates how many other processes havve finished sending their particles to this process
  //we wont stop until receive terminate symbols from all  processes (including itself)
  while (finished_processes != NUM_PROC){
    particle_t placeholder;
    //polling any incoming particles
    MPI_Recv(&placeholder, 1, PARTICLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
    //check if it is a terminate symbol //x == y == -1 is our termiante symbols
    if (placeholder.x == -1.0 || placeholder.y == -1.0){
      finished_processes++;
      if (DEBUG == 2)
        printf("Processor: %d: Finish polling from Processor %d \n", status.MPI_SOURCE);
      continue;
    }
    //otherwise, first we need to compute which block it belongs to
    int which_block_x = min((int)(placeholder.x / BLOCK_SIZE), NUM_BLOCKS_PER_DIM - 1);
    int which_block_y = min((int)(placeholder.y / BLOCK_SIZE), NUM_BLOCKS_PER_DIM - 1);
    //assert here
    if (DEBUG){
        assert((myInfo.start_row <= which_block_y) && (myInfo.end_row >= which_block_y));
        assert((myInfo.start_col <= which_block_x) && (myInfo.end_col >= which_block_x));
    }

    //insert
    myBlocks[which_block_y - myInfo.start_row][which_block_x - myInfo.start_col].particles.push_back(placeholder);
  }
  if (DEBUG == 2){
    printf("Processor %d: Finish polling particles from all other processes \n", RANK);
  }
}

//move the particles to appropriate blocks (can be remote) after computing the force
void move_particles(){
  ClusterInfo myInfo = cluster_layout[RANK];
  vector<vector<Block> > oldBlocks = myBlocks;
  vector<MPI_Request> requests;
  //clear everything
  for (int i = 0; i <= myInfo.end_row - myInfo.start_row; i++){
    for (int j = 0; j <= myInfo.end_col - myInfo.start_col; j++){
      myBlocks[i][j].particles.clear();
    }
  }
  //loop through the old block
  for (int i = 0; i <= myInfo.end_row - myInfo.start_row; i++){
    for (int j = 0; j <= myInfo.end_col - myInfo.start_col; j++){
      Block& oldBlock = oldBlocks[i][j];
      for (int k = 0; k < oldBlock.particles.size(); k++){
        double old_x = oldBlock.particles[k].x;
        double old_y = oldBlock.particles[k].y;
        //update position
        move(oldBlock.particles[k]);
        //check if the particle might move to another block
        MPI_Request request = decide_membership(myBlocks[i][j], old_x, old_y, oldBlock.particles[k]);
        if (request != NULL){
          requests.push_back(request);
        }
      }
    }
  }
  //let all the processes including itself knows that it is done
  for (int proc = 0; proc < NUM_PROC; proc++){
    MPI_Request request;
    MPI_Isend(&TERMINATE_SYMBOL, 1, PARTICLE, proc, TRANSFER_PARTICLE_TAG, MPI_COMM_WORLD, &request);
    requests.push_back(request);
  }
  //polling for the particles we need
  poll_particles(TRANSFER_PARTICLE_TAG);

  MPI_Waitall(requests.size(), &requests.front(), MPI_STATUSES_IGNORE);
  if (DEBUG == 2)
    printf("Process %d: Finish moving particles \n", RANK);
}

//Called in compute_force_grid
void compute_force_between_blocks(Block& block_A, Block& block_B, int& navg, double& davg, double& dmin){
  for (vector<particle_t>::iterator it_A = block_A.particles.begin(); it_A != block_A.particles.end(); it_A++){
    for (std::vector<particle_t>::iterator it_B = block_B.particles.begin(); it_B != block_B.particles.end(); it_B++){
      apply_force(*it_A, *it_B, &dmin, &davg, &navg);
    }
  }
}

//Called in compute_force_grid
void compute_force_within_block(Block& block, int& navg, double& davg, double& dmin){
    for (vector<particle_t>::iterator it_1 = block.particles.begin(); it_1 != block.particles.end(); it_1++){
      for (vector<particle_t>::iterator it_2 = block.particles.begin(); it_2 != block.particles.end(); it_2++){
        apply_force(*it_1, *it_2, &dmin, &davg, &navg);
      }
    }
}

//compute force within myblocks
//Called in simulate_particles
void compute_force_grid(int& navg, double& davg, double& dmin){
  ClusterInfo myInfo = cluster_layout[RANK];
  for (int i = 0; i <= myInfo.end_row - myInfo.start_row; i++){
    for (int j = 0; j <= myInfo.end_col - myInfo.start_col; j++){

      int navg_ = 0;
      double davg_ = 0.0;

      //set acceleration to zero
      for (vector<particle_t>::iterator it = myBlocks[i][j].particles.begin(); it != myBlocks[i][j].particles.end(); it++){
        it->ax = it->ay = 0;
      }

      //check top
      if (i != 0){
        compute_force_between_blocks(myBlocks[i][j], myBlocks[i-1][j], navg_, davg_, dmin);
      }else if (i == 0 && myInfo.start_row != 0){
        compute_force_between_blocks(myBlocks[i][j], topEdge[j], navg_, davg_, dmin);
      }
      //check bot
      if (i != myInfo.end_row - myInfo.start_row){
        compute_force_between_blocks(myBlocks[i][j], myBlocks[i+1][j], navg_, davg_, dmin);
      }else if (i == myInfo.end_row - myInfo.start_row && myInfo.end_row != NUM_BLOCKS_PER_DIM - 1){
        compute_force_between_blocks(myBlocks[i][j], botEdge[j], navg_, davg_, dmin);
      }
      //check diagonal left top
      if (j == 0){
      }else if (i == 0 && myInfo.start_row != 0){
        compute_force_between_blocks(myBlocks[i][j], topEdge[j-1], navg_, davg_, dmin);
      }else if (i != 0 && j != 0){
        compute_force_between_blocks(myBlocks[i][j], myBlocks[i-1][j-1], navg_, davg_, dmin);
      }
      //check diagonal left bot
      if (j == 0){
      }else if (i == myInfo.end_row - myInfo.start_row && myInfo.end_row != NUM_BLOCKS_PER_DIM - 1){
        compute_force_between_blocks(myBlocks[i][j], botEdge[j-1], navg_, davg_, dmin);
      }else if (i != myInfo.end_row - myInfo.start_row && j != 0){
        compute_force_between_blocks(myBlocks[i][j], myBlocks[i+1][j-1], navg_, davg_, dmin);
      }
      //check diagonal right top
      if (j == myInfo.end_col - myInfo.start_col){
      }else if (i == 0 && myInfo.start_row != 0){
        compute_force_between_blocks(myBlocks[i][j], topEdge[j+1], navg_, davg_, dmin);
      }else if (i != 0 && j != myInfo.end_col - myInfo.start_col){
        compute_force_between_blocks(myBlocks[i][j], myBlocks[i-1][j+1], navg_, davg_, dmin);
      }
      //check diagonal right bot
      if (j == myInfo.end_col - myInfo.start_col){
      }else if (i == myInfo.end_row - myInfo.start_row && myInfo.end_row != NUM_BLOCKS_PER_DIM - 1){
        compute_force_between_blocks(myBlocks[i][j], botEdge[j+1], navg_, davg_, dmin);
      }else if (i != myInfo.end_row - myInfo.start_row && j != myInfo.end_col - myInfo.start_col){
        compute_force_between_blocks(myBlocks[i][j], myBlocks[i+1][j+1], navg_, davg_, dmin);
      }
      //check left
      if (j != 0){
        compute_force_between_blocks(myBlocks[i][j], myBlocks[i][j-1], navg_, davg_, dmin);
      }
      //check right
      if (j != myInfo.end_col - myInfo.start_col){
        compute_force_between_blocks(myBlocks[i][j], myBlocks[i][j+1], navg_, davg_, dmin);
      }

      //compute within itself
      compute_force_within_block(myBlocks[i][j], navg_, davg_, dmin);
      navg += navg_;
      davg += davg_;
    }
  }
}

/*-----------------Simulate Function---------------------*/
void simulate_particles(char** argv, int argc, particle_t* particles, int n,
    int& navg, double& davg, double& dmin, double& rdavg, double& rdmin,
      int& rnavg, int& nabsavg, double& absavg, double& absmin, FILE* fsave){

  TERMINATE_SYMBOL.x = -1.0;
  TERMINATE_SYMBOL.y = -1.0;

  if(RANK == MASTER){
    master_routine(particles, n);
  }else{
    receive_metaData_from_master(MASTER);
    receive_clusterInfo_from_master(MASTER);
    receive_blocks_from_master(MASTER, BLOCKS_INITIALIZATION_TAG);
  }

  //debug
  // if (DEBUG){
  //   printSummary();
  //   printBlocks();
  // }
  for( int step = 0; step < NSTEPS; step++ ){
    navg = 0;
    dmin = 1.0;
    davg = 0.0;

    if( find_option( argc, argv, "-no" ) == -1 )
      if( fsave && (step%SAVEFREQ) == 0 )
        save( fsave, n, particles );

    //printf("Processor %d: Step %d\n", RANK, step);
    request_and_feed_edges(REQUEST_AND_FEED_EDGES_TAG);
    compute_force_grid(navg, davg, dmin);
    move_particles();

    if( find_option( argc, argv, "-no" ) == -1 ){
      MPI_Reduce(&davg,&rdavg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(&navg,&rnavg,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(&dmin,&rdmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
      if (RANK == MASTER){
        if (rnavg) {
          absavg +=  rdavg/rnavg;
          nabsavg++;
        }
        if (rdmin < absmin) absmin = rdmin;
      }
    }
  }
}

int main( int argc, char **argv ){
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

    FILE *fsave = savename && RANK == 0 ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname && RANK == 0 ? fopen ( sumname, "a" ) : NULL;

    double simulation_time = read_timer();
    simulate_particles(argv, argc, particles, n, navg, davg, dmin, rdavg, rdmin, rnavg, nabsavg, absavg, absmin, fsave);
    simulation_time = read_timer( ) - simulation_time;

    if (RANK == MASTER) {
      printf( "n = %d, simulation time = %g seconds", n, simulation_time);
      if( find_option( argc, argv, "-no" ) == -1 ){
        if (nabsavg){
          absavg /= nabsavg;
        }
        //
        //  -The minimum distance absmin between 2 particles during the run of the simulation
        //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
        //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
        //
        //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
        //
        printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
        if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
        if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
      }
      printf("\n");
      if( fsum)
        fprintf(fsum,"%d %d %g\n",n,NUM_PROC,simulation_time);
    }
    if ( fsum )
        fclose( fsum );
    if( fsave )
        fclose( fsave );

    free( particles );

    MPI_Finalize( );

    return 0;
}
