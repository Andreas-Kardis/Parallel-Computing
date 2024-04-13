#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <unistd.h>
#include <mpi.h>
#include <omp.h>

#define SHIFT_ROW 0
#define SHIFT_COL 1
#define DISP 1
#define NUM_PORTS 3
#define MAX_LOGS 5
#define CURR_MAX_DATA 14
#define CYCLE 1
#define MAX_RUN_COUNT 45
#define NUM_REQ_INFO 4

#define MSG_EXIT 1
#define MSG_FULL_PORTS 2
#define MSG_CN_INFO 3
#define MSG_REQ 4
#define MSG_AVAILABILITIES 5

int server_io(MPI_Comm master_comm, MPI_Comm comm);
int charging_node_io(MPI_Comm master_comm, MPI_Comm comm, int argc, char *argv[]);
struct shared_array_struct;
void init_shared_array(struct shared_array_struct shared_array);
void add(struct shared_array_struct *shared_array, int *data);
int write_to_file(int *CN_log, int **all_neighbours, int **coords, int num_comm, double comm_time, int iteration);
int append_to_file(int node, int **all_neighbours, int *closest_neighbours, int iteration, int **node_availabilities);

struct shared_array_struct {
  int *data[MAX_LOGS];
  int start;
  int end;
  int count;
};

int main(int argc, char *argv[]) {
  int world_rank, size;
  
  // Initialising seperation between server and wireless sensor network
  MPI_Comm new_comm;
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  if (argc == 3) {
    int nrows = atoi(argv[1]);
    int ncols = atoi(argv[2]);
    if (nrows*ncols != (size-1)) {
      if (world_rank == 0) {
        printf("dimensions or size are incorrectly chosen, %d * %d != %d.\n", nrows, ncols, size-1);
      }

      MPI_Finalize();
      return 0;
    }
  }

  MPI_Comm_split(MPI_COMM_WORLD, world_rank == 0, 0, &new_comm);

  if (world_rank == 0) {
    server_io(MPI_COMM_WORLD, new_comm);
  } else {
    charging_node_io(MPI_COMM_WORLD, new_comm, argc, argv);
  }

  MPI_Finalize();
  return 0;
}

int server_io(MPI_Comm master_comm, MPI_Comm comm) {
  int world_size, cart_size, i, j;
  int iteration = 0, receiving_logs = 0;

  MPI_Comm_size(master_comm, &world_size);
  cart_size = world_size - 1;

  int received_logs[cart_size][NUM_REQ_INFO * 4 + CURR_MAX_DATA];
  int **charging_node_neighbours = (int **)malloc(cart_size * sizeof(int *)); // contains neighbours where each row is a subsequent rank starting from 0
  int **coords = (int **)malloc(cart_size * sizeof(int *)); // coords of each charging node

  int available_neighbours[cart_size]; // contains the iteration at which time the neighbour became full
  memset(available_neighbours, 0, sizeof(available_neighbours));

  int full_nodes[cart_size]; // how many times the node has reported full without being reported to file yet
  memset(full_nodes, 0, sizeof(full_nodes));

  int iteration_full[cart_size]; // at which iteration the node reported full
  memset(iteration_full, 0, sizeof(iteration_full));

  int num_comm[cart_size]; // number of times each node has communicated with the server
  for (i = 0; i < cart_size; i++) {
    num_comm[i] = 1;
  }

  int iteration_new[cart_size]; // new iteration if the same node has reported full again without it being written to file yet
  memset(iteration_new, 0, sizeof(iteration_new));

  double communication_time[cart_size];
  memset(communication_time, 0, sizeof(communication_time));
  double start, end;

  MPI_Status CN_info_status[cart_size];

  int tmp_neighbours[6];

  // receive initial data from each node inidicating its coordinates and neighbours
  #pragma omp parallel private(i, j, tmp_neighbours) shared(charging_node_neighbours) num_threads(1) 
  {
    // receive charging node info
    for (i = 0; i < cart_size; i++) {
      MPI_Recv(tmp_neighbours, 6, MPI_INT, i+1, MSG_CN_INFO, master_comm, &CN_info_status[i]);

      charging_node_neighbours[i] = (int *)malloc(4 * sizeof(int));
      for (j = 0; j < 4; j++) {
        charging_node_neighbours[i][j] = tmp_neighbours[j+2];
      }

      coords[i] = (int *)malloc(2 * sizeof(int));
      for (j = 0; j < 2; j++) {
        coords[i][j] = tmp_neighbours[j];
      }
    }
  }

  // while there are still charging nodes to consider
  while(cart_size > 0) {
    iteration++;
    MPI_Request receiving_logs_request[cart_size];
    MPI_Status probe_status[cart_size];

    #pragma omp parallel private(i, j, start, end) num_threads(1)
    {
      // receive logs
      for (i = 0; i < cart_size; i++) {
        MPI_Iprobe(i+1, MSG_FULL_PORTS, master_comm, &receiving_logs, &probe_status[i]); // probe each wireless sensor node
        printf("Probe has reported on cart rank %d the number %d.\n", i, receiving_logs);

        if (receiving_logs) {
          start = MPI_Wtime();
          MPI_Irecv(received_logs[i], NUM_REQ_INFO * 4 + CURR_MAX_DATA, MPI_INT, i+1, MSG_FULL_PORTS, master_comm, &receiving_logs_request[i]);
          MPI_Wait(&receiving_logs_request[i], MPI_STATUS_IGNORE);
          end = MPI_Wtime();

          // recording information about reporting node
          num_comm[i]++;
          communication_time[i] = end - start;

          // update when the node requested for server's help
          if (full_nodes[i] == 0) {
            full_nodes[i]++;
            iteration_full[i] = iteration;
          } else {
            full_nodes[i]++;
            iteration_new[i] = iteration;
          }

          // indicate that this node has reported full to other nodes
          available_neighbours[i] = iteration;

        }
      }      
    }

    // write log to file
    for (i = 0; i < cart_size; i++) {
      if ((full_nodes[i] > 0) && ((iteration_full[i] + 2) == iteration)) { // if there is a node waiting to be written to file and waited two iterations
        // check which neighbours are available if it is closest_neighbours will represent its rank
        int closest_neighbours[4];
        for (j = 0; j < 4; j++) {
          int neighbour = charging_node_neighbours[i][j];
          if (neighbour != -2) {
            if (iteration_full[i] > available_neighbours[neighbour]) {
              closest_neighbours[j] = charging_node_neighbours[i][j];
            } else {
              closest_neighbours[j] = -2;
            }
          } else {
            closest_neighbours[j] = -2;
          }
        }

        MPI_Request availability_req;
        int *node_availabilities = NULL;
        int length;
        
        length = append_to_file(i, charging_node_neighbours, closest_neighbours, iteration-2, &node_availabilities);

        // send node's availability
        #pragma omp parallel num_threads(1)
        {
          MPI_Isend(node_availabilities, length, MPI_INT, i+1, MSG_AVAILABILITIES, master_comm,& availability_req);
        }

        // moving second full report from node i into primary arrays
        if (full_nodes[i] > 1) {
          iteration_full[i] = iteration_new[i];
          iteration_new[i] = 0;
        }

        full_nodes[i]--;
      } else if ((full_nodes[i] > 0) && (iteration_full[i] == iteration)) { // do initial writing into file
        write_to_file(received_logs[i], charging_node_neighbours, coords, num_comm[i], communication_time[i], iteration_full[i]);
      }
    }
    
    sleep(CYCLE);

    // send termination messages.
    #pragma omp parallel private(i) num_threads(1)
    {
      if (iteration == MAX_RUN_COUNT) {
        MPI_Request exit_req[cart_size];

        // send termination messages to each node
        for (i = 0; i < cart_size; i++) {
          MPI_Isend(&iteration, 1, MPI_INT, i+1, MSG_EXIT, master_comm, &exit_req[i]);
        }
        cart_size = 0;
      }
    }
  }

  return 0;
}

int charging_node_io(MPI_Comm master_comm, MPI_Comm comm, int argc, char *argv[]) {
  int ndims = 2, size, rank, cart_rank, nrows, ncols, reorder, ierr;
  int sum_occupied = 0;
  int i, num;
  int coord[ndims];
  int wrap_around[ndims], dims[ndims];
  int nbr_i_lo, nbr_i_hi, nbr_j_lo, nbr_j_hi;
  time_t t;
  struct tm local_time;

  MPI_Comm comm2D;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  // obtaining dimensions of cartesian topology
  if (argc == 3) {
    nrows = atoi(argv[1]);
    ncols = atoi(argv[2]);
    dims[0] = nrows;
    dims[1] = ncols;
  } else {
    nrows = ncols = (int) sqrt(size);
    dims[0] = dims[1] = 0;
    MPI_Dims_create(size, ndims, dims);
  }

  // create cartesian topology
  wrap_around[0] = wrap_around[1] = 0;
  reorder = 1;
  ierr = 0;
  ierr = MPI_Cart_create(comm, ndims, dims, wrap_around, reorder, &comm2D);
  // check if mpi_cart_create was done correctly
  if (ierr != 0) printf("error[%d] creating cart\n", ierr);

  // retreiving cart coords and cart rank
  MPI_Cart_coords(comm2D, rank, ndims, coord);
  MPI_Cart_rank(comm2D, coord, &cart_rank);

  // get neighbour's coords and rank
  MPI_Cart_shift(comm2D, SHIFT_ROW, DISP, &nbr_i_lo, &nbr_i_hi);
  MPI_Cart_shift(comm2D, SHIFT_COL, DISP, &nbr_j_lo, &nbr_j_hi);
  int neighbours[4] = {nbr_i_lo, nbr_i_hi, nbr_j_lo, nbr_j_hi};

  int neighbour_info[6] = {coord[0], coord[1], nbr_i_lo, nbr_i_hi, nbr_j_lo, nbr_j_hi};
  MPI_Send(neighbour_info, 6, MPI_INT, 0, MSG_CN_INFO, master_comm);

  // begin construction of the shared array 
  struct shared_array_struct shared_array;
  //init_shared_array(shared_array);
  shared_array.start = 0;
  shared_array.end = -1;
  shared_array.count = 0;

  for (int i = 0; i < MAX_LOGS; i++) {
    shared_array.data[i] = (int *)malloc(CURR_MAX_DATA * sizeof(int));
  }

  MPI_Datatype shared_array_type;
  MPI_Datatype type[4] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  int blocklen[4] = {MAX_LOGS, 1, 1, 1};
  MPI_Aint displacement[4];

  // obtain the address of each element in the struct
  MPI_Get_address(shared_array.data, &displacement[0]);
  MPI_Get_address(&shared_array.start, &displacement[1]);
  MPI_Get_address(&shared_array.end, &displacement[2]);
  MPI_Get_address(&shared_array.count, &displacement[3]);

  // make the addresses relative to displacement[0]
  displacement[1] -= displacement[0];
  displacement[2] -= displacement[0];
  displacement[3] -= displacement[0];
  displacement[0] = 0;

  // create mpi struct
  MPI_Type_create_struct(4, blocklen, displacement, type, &shared_array_type);
  MPI_Type_commit(&shared_array_type);

  unsigned int seed = time(NULL) * cart_rank;

  int num_neighbours = 0;
  int real_neighbours[4];
  for (int i = 0; i < 4; i++) {
    if (neighbours[i] != -2) {
      real_neighbours[num_neighbours] = neighbours[i];
      num_neighbours++;
    }
  }

  // request and status statements to be used in varies send and recv functions
  MPI_Request recv_log_request[num_neighbours*2];
  MPI_Status recv_log_status[num_neighbours*2];
  MPI_Request send_log_request[num_neighbours];
  MPI_Request send_log[num_neighbours];
  MPI_Status send_log_status[num_neighbours*2];

  MPI_Request send_server_request;

  MPI_Request availability_req;
  int availabilities[16] = {-1};

  // intiialising variables
  int req_cart_rank = -1, receive_completed = 0, exit_message = 0, num_full_neighbours = 0, iteration = 0;
  int closest_neighbour;
  int closest_neighbour_event[MAX_RUN_COUNT] = {0};
  int requested_logs[num_neighbours][NUM_REQ_INFO];
  bool requesting_neighbour_log = false, exit = false;

  // while the sever has not sent a termination message
  while (!exit) {
    iteration++;
    // simulate ports for this node
    #pragma omp parallel for private(i) schedule(static) reduction(+ : sum_occupied) num_threads(NUM_PORTS)
    for (i=0; i < NUM_PORTS; i++) {
      num = rand_r(&seed) % 2;
      sum_occupied += num; // num of 1, means the port is occupied
    }

    t = time(NULL);
    local_time = *localtime(&t);

    MPI_Barrier(comm2D);

    // add to a log
    int log[CURR_MAX_DATA] = {NUM_PORTS-sum_occupied, cart_rank, coord[0], coord[1], local_time.tm_year+1900, local_time.tm_mon+1, local_time.tm_mday, local_time.tm_hour, local_time.tm_min, local_time.tm_sec, -1, -1, -1, -1};

    // if ports are at least almost full, prompt for neighbour data.
    if (sum_occupied >= NUM_PORTS - 1) {
      requesting_neighbour_log = true;
      // prompt for neighbour data
      for (int i = 0; i < num_neighbours; i++) {
        //printf("%d requesting %d's info\n", cart_rank, real_neighbours[i]);
        MPI_Isend(&cart_rank, 1, MPI_INT, real_neighbours[i], MSG_REQ, comm2D, &send_log_request[i]); // send cart_rank to each neighbour
      }
    }

    MPI_Barrier(comm2D);

    // test to see if a receive call is produced from each neighbour
    int neighbour_data[NUM_REQ_INFO] = {log[0], log[1], log[2], log[3]};
    for (int i = 0; i < num_neighbours; i++) {

      MPI_Iprobe(real_neighbours[i], MSG_REQ, comm2D, &receive_completed, &send_log_status[i]);

      if (receive_completed) {
        MPI_Irecv(&req_cart_rank, 1, MPI_INT, real_neighbours[i], MSG_REQ, comm2D, &recv_log_request[i]); // receive its cart rank
        MPI_Isend(&neighbour_data, NUM_REQ_INFO, MPI_INT, real_neighbours[i], 0, comm2D, &send_log[i]); // send this node's logs
        //printf("%d receiving data from %d and sending to %d. AND REQ_CART_RANK IS %d\n", cart_rank, real_neighbours[i], real_neighbours[i], req_cart_rank);
      }
    }

    MPI_Barrier(comm2D);

    // receive log from neighbour with vacant ports
    if (requesting_neighbour_log) {
      for (int i = 0; i < num_neighbours; i++) {
        MPI_Recv(&requested_logs[i], NUM_REQ_INFO, MPI_INT, real_neighbours[i], 0, comm2D, &recv_log_status[i]); // receive neighbour log
        printf("%d received requested log info: rank %d, coord (%d, %d), availability %d\n", cart_rank, requested_logs[i][1], requested_logs[i][2], requested_logs[i][3], requested_logs[i][0]);
      }
    }

    printf("\n%d's number of available ports are: %d.\n", cart_rank, log[0]);

    MPI_Barrier(comm2D);

    // check if neighbours have vacant ports
    if (requesting_neighbour_log) {
      for (int i = 0; i < num_neighbours; i++) {
        //printf("%d's neighbour %d has %d vacant ports.\n", cart_rank, real_neighbours[i], requested_logs[i][0]);
        if (requested_logs[i][0] <= 1) { // there is not enough ports in this neighbour
          num_full_neighbours++;
        } else {
          closest_neighbour = i;
          break;
        }
      }

      // neighbours don't have vacant ports
      if (num_neighbours == num_full_neighbours) {
        printf("\n%d's neighbours have no vacancies\n\n", cart_rank);

        // send logs to server to request other charging nodes
        // first make logs into 1D list
        int flattened_server_logs[num_neighbours * NUM_REQ_INFO + CURR_MAX_DATA];
        for (int i = 0; i < CURR_MAX_DATA; i++) {
          flattened_server_logs[i] = log[i];
        }

        for (int i = 0; i < num_neighbours; i++) {
          for (int j = 0; j < NUM_REQ_INFO; j++) {
            flattened_server_logs[i * NUM_REQ_INFO + j + CURR_MAX_DATA] = requested_logs[i][j];
          }
        }

        MPI_Isend(flattened_server_logs, num_neighbours * NUM_REQ_INFO + CURR_MAX_DATA, MPI_INT, 0, MSG_FULL_PORTS, master_comm, &send_server_request);

        MPI_Irecv(availabilities, 16, MPI_INT, 0, MSG_AVAILABILITIES, master_comm, &availability_req);
      } else {
        // add closest neighbour data indicating which port to go to
        for (int i = 0; i < 4; i++) {
          log[10+i] = requested_logs[closest_neighbour][i];
        }
        
        closest_neighbour_event[iteration]++;
        printf("%d's closest neighbour has characteristics %d ports available, rank: %d, coord: (%d, %d)\n", cart_rank, log[10], log[11], log[12], log[13]);
      }
    }
    
    MPI_Barrier(comm2D);

    // adding the whole log into the shared_array
    add(&shared_array, log);

    // reset requested_logs
    for (int row = 0; row < num_neighbours; row++) {
      for (int col = 0; col < CURR_MAX_DATA; col++) {
        requested_logs[row][col] = -1;
      }
    }
    
    // reset variables
    num_full_neighbours = 0;
    receive_completed = 0;
    requesting_neighbour_log = false;
    sum_occupied = 0;

    MPI_Barrier(comm2D);

    // wait to re-simulate node
    sleep(CYCLE);

    if (cart_rank == 0){
      printf("\n########################################\n\n");
    }

    // check if there is a termination message from server
    MPI_Iprobe(0, MSG_EXIT, master_comm, &exit_message, MPI_STATUS_IGNORE);

    if(exit_message) {
      exit = true;
      MPI_Status status;
      int last_message;
      MPI_Recv(&last_message, 1, MPI_INT, 0, MSG_EXIT, master_comm, &status);
      
      char file_name[32];
      sprintf(file_name, "node_%d_event1", cart_rank);
      FILE *file = fopen(file_name, "w");
      for (int i = 0; i < MAX_RUN_COUNT; i++) {
        fprintf(file, "%d, ", closest_neighbour_event[i]);
      }
      fclose(file);
    }

    MPI_Barrier(comm2D);
  }

  return 0;
}

/*
* adds data to the shared array
*/
void add(struct shared_array_struct *shared_array, int *data) {
  // if array is not full
  if (shared_array->count < MAX_LOGS) {
    shared_array->end = (shared_array->end + 1) % MAX_LOGS;
    for (int i = 0; i < CURR_MAX_DATA; i++) {
      shared_array->data[shared_array->end][i] = data[i];
    }

    shared_array->count++;
  } else { // array is full
    shared_array->start = (shared_array->start + 1) % MAX_LOGS;
    shared_array->end = (shared_array->end + 1) % MAX_LOGS;   

    for (int i = 0; i < CURR_MAX_DATA; i++) {
      shared_array->data[shared_array->end][i] = data[i];
    }
  }
}

/*
* does the initial writing to the file, hereby creates the file.
*/
int write_to_file(int *CN_log, int **all_neighbours, int **coords, int num_comm, double comm_time, int iteration) {
  char file_name[32];
  int *neighbours = all_neighbours[CN_log[1]];
  time_t t;
  struct tm local_time;

  int num_neighbours = 0;
  for (int i = 0; i < 4; i++) {
    if (neighbours[i] != -2) {
      num_neighbours++;
    }
  }

  t = time(NULL);
  local_time = *localtime(&t);
  
  sprintf(file_name, "KPM_iter_%d_node_%d", iteration, CN_log[1]);
  FILE *file = fopen(file_name, "w");

  fprintf(file, "Iteration: %d\n", iteration);
  fprintf(file, "Communication Time (seconds): %.8f\n", comm_time);
  fprintf(file, "Time charging station became full: %d-%d-%d %d:%d:%d\n", CN_log[4], CN_log[5], CN_log[6], CN_log[7], CN_log[8], CN_log[9]);
  fprintf(file, "Reported time:                     %d-%d-%d %d:%d:%d\n", local_time.tm_year+1900, local_time.tm_mon+1, local_time.tm_mday, local_time.tm_hour, local_time.tm_min, local_time.tm_sec);
  fprintf(file, "Total messages sent by the reporting node to the server: %d\n", num_comm);
  fprintf(file, "Number of adjacent node(s): %d\n\n", num_neighbours);
  fprintf(file, "Reporting node      Coord       Port Value      Available Ports\n");
  fprintf(file, "%d                  (%d, %d)       %d               %d\n\n", CN_log[1], CN_log[2], CN_log[3], NUM_PORTS, CN_log[0]);
  fprintf(file, "Adjacent node(s)    Coord       Port Value      Available Ports\n");

  int real_neighbours[num_neighbours];
  for (int i = 0; i < num_neighbours; i++) {
    int idx = CURR_MAX_DATA + NUM_REQ_INFO * i;
    real_neighbours[i] = CN_log[idx+1];

    fprintf(file, "%d                  (%d, %d)       %d               %d\n", CN_log[idx+1], CN_log[idx+2], CN_log[idx+3], NUM_PORTS, CN_log[idx]);
  }
  fprintf(file, "\n");

  fprintf(file, "Nearby node(s)      Coord\n");

  int *considered_nodes = (int *)malloc(num_neighbours*4 * sizeof(int));
  considered_nodes[0] = -2;
  considered_nodes[1] = CN_log[1];
  int idx = 2;
  for (int i = 0; i < num_neighbours; i++) {
    int *node_neighbour = all_neighbours[real_neighbours[i]];
    for (int j = 0; j < 4; j++) {
      int node = node_neighbour[j];
      bool found = false;

      for (int k = 0; k < idx; k++) {
        if (node == considered_nodes[k]) {
          found = true;
        }
      }

      if (!found) {
        considered_nodes[idx] = node;
        idx++;
        fprintf(file, "%d                  (%d, %d)\n", node, coords[node][0], coords[node][1]);
      }
    }
  }

  fprintf(file, "\n");

  fclose(file);

  return 0;
}


/*
* Appends to the already created file.
*/
int append_to_file(int node, int **all_neighbours, int *closest_neighbours, int iteration, int **node_availabilities) {
  char file_name[32];
  int *available_neighbours = (int *)malloc(16 * sizeof(int)); // saves the available neighbours identified
  available_neighbours[0] = -1;
  int AN_idx = 0;

  int *considered_nodes = (int *)malloc(16 * sizeof(int)); // saves considered nodes
  considered_nodes[0] = -2;
  considered_nodes[1] = node;
  int CN_idx = 2;

  // identifies all available nodes by checking if they are already considered and a valid neighbour
  for (int i = 0; i < 4; i++) {
    int neighbour = closest_neighbours[i];

    if (neighbour != -2) {
      for (int j = 0; j < 4; j++) {
        int distant_neighbour = all_neighbours[neighbour][j];
        bool found = false;

        for (int k = 0; k < CN_idx; k++) {
          if (distant_neighbour == considered_nodes[k]) {
            found = true;
          }
        }

        if (!found) {
          available_neighbours[AN_idx] = distant_neighbour;
          AN_idx++;

          considered_nodes[CN_idx] = distant_neighbour;
          CN_idx++;
        }
      }
    }
  }

  sprintf(file_name, "KPM_iter_%d_node_%d", iteration, node);
  FILE *file = fopen(file_name, "a");

  if (available_neighbours[0] != -1) {
    fprintf(file, "Available charging station(s) nearby (no report received in last 2 iterations): %d", available_neighbours[0]);

    *node_availabilities = (int *)malloc(AN_idx * sizeof(int));
    (*node_availabilities)[0] = available_neighbours[0];
    for (int i = 1; i < AN_idx; i++) {
      fprintf(file, ", %d", available_neighbours[i]);
      (*node_availabilities)[i] = available_neighbours[i];
    }
  } else {
    fprintf(file, "Available charging station(s) nearby (no report received in last 2 iterations): NULL");
  }

  fprintf(file, "\n");

  fclose(file);
  free(available_neighbours);
  return AN_idx;
}
