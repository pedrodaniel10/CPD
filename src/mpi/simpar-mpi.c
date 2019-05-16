#include <mpi.h>
#include "array_list.h"
#include "structures.h"

#define NUM_PARTICLES_BUFFER 10000

#define CONVERT_TO_LOCAL(id, p, n, x) (x - BLOCK_LOW(id, p, n) + 1)

enum { TAG_INIT_PARTICLES };

array_list_t* particles;
cell_t** cells;

int myRank, number_processors;
int size_processor_grid[2] = {0, 0};
int my_coordinates[2];
node_t* adjacent_processes[8] = {0};

int size_local_cell_matrix[2];
MPI_Comm cart_comm;

void create_cartesian_communicator(int grid_size) {
	// Calculate best grid size for the number of processes
	int periodic[2] = {1, 1};

	if (number_processors <= 3) {
		// Generate the best fit
		MPI_Dims_create(number_processors, 2, size_processor_grid);
	} else {
		size_processor_grid[1] = sqrt(number_processors);
		if (size_processor_grid[1] >= grid_size) {
			size_processor_grid[0] = size_processor_grid[1] = grid_size;
		} else {
			size_processor_grid[0] = number_processors / size_processor_grid[1];
		}
	}

	MPI_Cart_create(MPI_COMM_WORLD, 2, size_processor_grid, periodic, 1, &cart_comm);

	if (size_processor_grid[0] * size_processor_grid[1] <= myRank) {
		printf("[%d] EXITING...\n", myRank);
		fflush(stdout);
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		exit(0);
	}

	MPI_Comm_rank(cart_comm, &myRank);
	MPI_Cart_coords(cart_comm, myRank, 2, my_coordinates);
	printf("[%d] My rank is %d...\n", myRank, myRank);

	// Calcualte rank adjacent processes
	int adjacent_processes_rank[8];
	int counter = 0;
	for (int i = -1; i <= 1; i++) {
		for (int j = -1; j <= 1; j++) {
			if (i == 0 && j == 0) {
				continue;
			}
			int coords[2];
			coords[0] = my_coordinates[0] + i;
			coords[1] = my_coordinates[1] + j;
			MPI_Cart_rank(cart_comm, coords, &adjacent_processes_rank[counter]);
			counter++;
		}
	}

	// Set strcutre for adjacent cells
	for (int i = 0; i < 8; i++) {
		if (adjacent_processes[i] != NULL) {
			continue;
		}
		adjacent_processes[i] = (node_t*)calloc(1, sizeof(node_t));
		adjacent_processes[i]->rank = adjacent_processes_rank[i];
		for (int j = i + 1; j < 8; j++) {
			if (adjacent_processes_rank[j] == adjacent_processes_rank[i]) {
				adjacent_processes[j] = adjacent_processes[i];
			}
		}
	}
}

void init_cells(int grid_size) {
	size_local_cell_matrix[0] = BLOCK_SIZE(my_coordinates[0], size_processor_grid[0], grid_size) + 2;
	size_local_cell_matrix[1] = BLOCK_SIZE(my_coordinates[1], size_processor_grid[1], grid_size) + 2;

	cells = (cell_t**)malloc(sizeof(cell_t*) * size_local_cell_matrix[0]);
	cell_t* cells_chunk = (cell_t*)calloc(size_local_cell_matrix[0] * size_local_cell_matrix[1], sizeof(cell_t));

	for (int i = 0; i < size_local_cell_matrix[0]; i++) {
		cells[i] = &cells_chunk[i * size_local_cell_matrix[1]];
	}
}

void init_particles(long seed, long grid_size, long long number_particles) {
	long long i;
	int number_processors_grid = size_processor_grid[0] * size_processor_grid[1];
	int* counters = (int*)calloc(number_processors_grid, sizeof(int));
	particle_t* buffer_space =
	    (particle_t*)malloc(sizeof(particle_t) * NUM_PARTICLES_BUFFER * (number_processors_grid));
	particle_t** buffers = (particle_t**)malloc(sizeof(particle_t*) * (number_processors_grid));

	for (int i = 1; i < number_processors_grid; i++) {
		buffers[i - 1] = &buffer_space[NUM_PARTICLES_BUFFER * (i - 1)];
	}

	srandom(seed);

	for (i = 0; i < number_particles; i++) {
		particle_t particle;

		particle.index = i;
		particle.position.x = RND0_1;
		particle.position.y = RND0_1;
		particle.velocity.x = RND0_1 / grid_size / 10.0;
		particle.velocity.y = RND0_1 / grid_size / 10.0;

		particle.mass = RND0_1 * grid_size / (G * 1e6 * number_particles);

		int cell_coordinate_x = particle.position.x * grid_size;
		int cell_coordinate_y = particle.position.y * grid_size;
		int coords_proc_grid[2];
		coords_proc_grid[0] = BLOCK_OWNER(cell_coordinate_x, size_processor_grid[0], grid_size);
		coords_proc_grid[1] = BLOCK_OWNER(cell_coordinate_y, size_processor_grid[1], grid_size);

		int proc_id_to_send;
		MPI_Cart_rank(cart_comm, coords_proc_grid, &proc_id_to_send);

		if (proc_id_to_send == 0) {
			append(particles, particle);
			continue;
		}

		// Add to process' buffer
		buffers[proc_id_to_send - 1][counters[proc_id_to_send - 1]] = particle;
		counters[proc_id_to_send - 1]++;

		// Buffer size reach => send
		if (counters[proc_id_to_send - 1] == NUM_PARTICLES_BUFFER) {
			MPI_Send(buffers[proc_id_to_send - 1], SIZEOF_PARTICLE(NUM_PARTICLES_BUFFER), MPI_DOUBLE, proc_id_to_send,
			         TAG_INIT_PARTICLES, cart_comm);
			counters[proc_id_to_send - 1] = 0;
		}
	}

	// Send eveything that is not been sent yet.
	for (int i = 1; i < number_processors_grid; i++) {
		if (counters[i - 1] != 0) {
			MPI_Send(buffers[i - 1], SIZEOF_PARTICLE(counters[i - 1]), MPI_DOUBLE, i, TAG_INIT_PARTICLES, cart_comm);
		} else {
			double endFlag[1] = {-1};
			MPI_Send(endFlag, 1, MPI_DOUBLE, i, TAG_INIT_PARTICLES, cart_comm);
		}
	}

	free(counters);
	free(buffers[0]);
	free(buffers);
}

void receiveParticles() {
	int number_elements_received;
	MPI_Status status;
	while (1) {
		MPI_Probe(0, TAG_INIT_PARTICLES, cart_comm, &status);
		MPI_Get_count(&status, MPI_DOUBLE, &number_elements_received);

		if (number_elements_received == 1) {
			break;
		}

		MPI_Recv((void*)allocate_array(particles, GET_NUMBER_PARTICLE(number_elements_received)), NUM_PARTICLES_BUFFER,
		         MPI_DOUBLE, 0, TAG_INIT_PARTICLES, cart_comm, &status);

		if (GET_NUMBER_PARTICLE(number_elements_received) != NUM_PARTICLES_BUFFER) {
			break;
		}
	}
}

void calculate_centers_of_mass(int grid_size) {
	for (int i = 0; i < particles->length; i++) {
		particle_t* particle = list_get(particles, i);
		int global_cell_index_x = particle->position.x * grid_size;
		int global_cell_index_y = particle->position.y * grid_size;
		int local_cell_index_x =
		    CONVERT_TO_LOCAL(my_coordinates[0], size_processor_grid[0], grid_size, global_cell_index_x);
		int local_cell_index_y =
		    CONVERT_TO_LOCAL(my_coordinates[1], size_processor_grid[1], grid_size, global_cell_index_y);

		cell_t* cell = &cells[local_cell_index_x][local_cell_index_y];

		cell->mass_sum += particle->mass;

		cell->center_of_mass.x += particle->mass * particle->position.x;
		cell->center_of_mass.y += particle->mass * particle->position.y;
	}

	for (int i = 1; i <= size_local_cell_matrix[0] - 2; i++) {
		for (int j = 1; j <= size_local_cell_matrix[1] - 2; j++) {
			cell_t* cell = &cells[i][j];
			if (cell->mass_sum != 0) {
				cell->center_of_mass.x /= cell->mass_sum;
				cell->center_of_mass.y /= cell->mass_sum;
			}
		}
	}
}

void send_recv_centers_of_mass() {
	int cells_size_x = size_local_cell_matrix[0] - 2;
	int cells_size_y = size_local_cell_matrix[0] - 2;
	// Prepare send buffers
	for (int i = 0; i < 8; i++) {
		if (adjacent_processes[i]->cells_buffer_send == NULL) {
			adjacent_processes[i]->cells_buffer_send =
			    (cell_t*)calloc(cells_size_x * 2 + cells_size_y * 2 + 4, sizeof(cell_t));
		}
		if (adjacent_processes[i]->cells_buffer_recv == NULL) {
			adjacent_processes[i]->cells_buffer_recv =
			    (cell_t*)malloc(sizeof(cell_t) * (cells_size_x * 2 + cells_size_y * 2 + 4));
		}
	}

	for (int i = 7; i >= 0; i--) {
		switch (i) {
			case DIAGONAL_UP_LEFT:
				adjacent_processes[i]->cells_buffer_send[adjacent_processes[i]->length_send_buffer] = cells[1][1];
				adjacent_processes[i]->length_send_buffer++;
				break;
			case DIAGONAL_UP_RIGHT:
				adjacent_processes[i]->cells_buffer_send[adjacent_processes[i]->length_send_buffer] =
				    cells[1][cells_size_y];
				adjacent_processes[i]->length_send_buffer++;
				break;
			case DIAGONAL_DOWN_LEFT:
				adjacent_processes[i]->cells_buffer_send[adjacent_processes[i]->length_send_buffer] =
				    cells[cells_size_x][1];
				adjacent_processes[i]->length_send_buffer++;
				break;
			case DIAGONAL_DOWN_RIGHT:
				adjacent_processes[i]->cells_buffer_send[adjacent_processes[i]->length_send_buffer] =
				    cells[cells_size_x][cells_size_y];
				adjacent_processes[i]->length_send_buffer++;
				break;
			case LEFT:
				for (int i = 1; i <= cells_size_x; i++) {
					adjacent_processes[i]->cells_buffer_send[adjacent_processes[i]->length_send_buffer] = cells[i][1];
					adjacent_processes[i]->length_send_buffer++;
				}
				break;
			case RIGHT:
				for (int i = 1; i <= cells_size_x; i++) {
					adjacent_processes[i]->cells_buffer_send[adjacent_processes[i]->length_send_buffer] =
					    cells[i][cells_size_y];
					adjacent_processes[i]->length_send_buffer++;
				}
				break;
			case UP:
				for (int i = 1; i <= cells_size_y; i++) {
					adjacent_processes[i]->cells_buffer_send[adjacent_processes[i]->length_send_buffer] = cells[1][i];
					adjacent_processes[i]->length_send_buffer++;
				}
				break;
			case DOWN:
				for (int i = 1; i <= cells_size_y; i++) {
					adjacent_processes[i]->cells_buffer_send[adjacent_processes[i]->length_send_buffer] =
					    cells[cells_size_x][i];
					adjacent_processes[i]->length_send_buffer++;
				}
				break;

			default:
				printf("[%d] Default case in send send_recv_centers_of_mass\n", myRank);
				fflush(stdout);
				break;
		}
	}
}

int main(int argc, char* argv[]) {
	if (argc != 5) {
		printf("Expected 5 arguments, but %d were given\n", argc);
		exit(1);
	}
	int grid_size = atoi(argv[2]);
	int number_particles = atoi(argv[3]);
	int n_time_steps = atoi(argv[4]);

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &number_processors);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	particles = create_array_list(2048);

	create_cartesian_communicator(grid_size);
	init_cells(grid_size);

	if (myRank == 0) {
		init_particles(atoi(argv[1]), grid_size, number_particles);
	} else {
		receiveParticles();
	}

	for (int n = 0; n < n_time_steps; n++) {
		calculate_centers_of_mass(grid_size);
		send_recv_centers_of_mass();
		// calculate_new_iteration(particles, cells, grid_size, number_particles);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();
}