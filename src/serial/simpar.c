#include "../include/common.h"

cell_t** init_cells(int grid_size) {
    // Allocate Cells
    cell_t **cells_matrix = (cell_t **) malloc(sizeof(cell_t *) * grid_size);
    cell_t *cells_chunk = (cell_t *) calloc(grid_size * grid_size, sizeof(cell_t));

    for (int i = 0; i < grid_size; i++) {
        cells_matrix[i] = cells_chunk + (i * grid_size);
    }

    // Allocate Adjacent Cells Memory
    cell_t **adjacent_cells_chunk = (cell_t **) malloc(sizeof(cell_t*) * grid_size * grid_size * ADJACENT_CELLS_NUMBER);
    cell_t ***adjacent_cells_rows = (cell_t ***) malloc(sizeof(cell_t **) * grid_size * grid_size);
    adjacent_cells = (cell_t ****) malloc(sizeof(cell_t ***) * grid_size);

    // Initialize Adjacent Cells Memory
    for(int i = 0; i < grid_size; i++) {
        adjacent_cells[i] = &adjacent_cells_rows[i * grid_size];

        for(int j = 0; j < grid_size; j++) {
            adjacent_cells[i][j] = &adjacent_cells_chunk[(i * grid_size + j) * ADJACENT_CELLS_NUMBER];
            int index_adjacent_cells = 0;
            cell_t** adjacent_cell = adjacent_cells[i][j];

            for (int x = -1; x <= 1; x++) {
                for(int y = -1; y <= 1; y++) {
                    adjacent_cell[index_adjacent_cells] = &cells_matrix[(i + x + grid_size) % grid_size][(j + y + grid_size) % grid_size];
                    index_adjacent_cells++;
                }
            }
        }        
    }
    return cells_matrix;
}

void calculate_centers_of_mass(particle_t *particles, cell_t **cells, int grid_size, int number_particles) {
    for (int i = 0; i < number_particles; i++) {
        particle_t *particle = &particles[i];
        cell_t *cell = &cells[particle->cell.x][particle->cell.y];

        cell->mass_sum += particle->mass;

        cell->center_of_mass.x += particle->mass * particle->position.x;
        cell->center_of_mass.y += particle->mass * particle->position.y;
    }
    
    for (int i = 0; i < grid_size; i++) {
        for (int j = 0; j < grid_size; j++) {
            cell_t *cell = &cells[i][j];
            if (cell->mass_sum != 0) {
                cell->center_of_mass.x /= cell->mass_sum;
                cell->center_of_mass.y /= cell->mass_sum;
            }
        }
    }
}

void calculate_new_iteration(particle_t *particles, cell_t **cells, int grid_size, int number_particles) {
    for (int i = 0; i < number_particles; i++) {
        particle_t *particle = &particles[i];
        coordinate_t force = {0}, acceleration = {0};
        cell_t** particle_adjacent_cell = adjacent_cells[particle->cell.x][particle->cell.y];
        
        // Calculate force
        for (int i = 0; i < ADJACENT_CELLS_NUMBER; i++) {
            cell_t* adjacent_cell = particle_adjacent_cell[i];
            coordinate_t force_a_b;  
            
            force_a_b.x = adjacent_cell->center_of_mass.x - particle->position.x;
            force_a_b.y = adjacent_cell->center_of_mass.y - particle->position.y;
            double distance_squared = force_a_b.x * force_a_b.x + force_a_b.y * force_a_b.y;

            if (distance_squared < EPSLON * EPSLON) {
                continue;
            }

            double scalar_force = G * particle->mass * adjacent_cell->mass_sum / (distance_squared * sqrt(distance_squared));

            force.x += force_a_b.x * scalar_force;
            force.y += force_a_b.y * scalar_force;
        }

        // Calculate acceleration
        acceleration.x = force.x / particle->mass;
        acceleration.y = force.y / particle->mass;

        // Calculate new position
        particle->position.x += particle->velocity.x + acceleration.x * 0.5 + 1;
        particle->position.y += particle->velocity.y + acceleration.y * 0.5 + 1;

        particle->position.x -= (int) particle->position.x;
        particle->position.y -= (int) particle->position.y;

        // Calculate new velocity
        particle->velocity.x += acceleration.x;
        particle->velocity.y += acceleration.y;
         
        // Calculate new cell position
        particle->cell.x = particle->position.x * grid_size;
        particle->cell.y = particle->position.y * grid_size;
    }
}

coordinate_t calculate_overall_center_of_mass(particle_t* particles, int number_particles) {
    coordinate_t center_of_mass = {0};
    double total_mass = 0;
    for (int i = 0; i < number_particles; i++) {
        particle_t *particle = &particles[i];
        
        total_mass += particle->mass; 
        center_of_mass.x += particle->mass * particle->position.x;
        center_of_mass.y += particle->mass * particle->position.y;
    }
    center_of_mass.x /= total_mass;
    center_of_mass.y /= total_mass;

    return center_of_mass;
}

int main(int argc, const char** argv) {
    if (argc != 5) {
        printf("Expected 5 arguments, but %d were given\n", argc);
        exit(1);
    }

    int grid_size= atoi(argv[2]);
    int number_particles = atoi(argv[3]);
    int n_time_steps = atoi(argv[4]);

    // Initialize
    particle_t *particles = init_particles(atoi(argv[1]), grid_size, number_particles);
    cell_t **cells = init_cells(grid_size);

    for (int n = 0; n < n_time_steps; n++) {
        calculate_centers_of_mass(particles, cells, grid_size, number_particles);
        calculate_new_iteration(particles, cells, grid_size, number_particles);

        memset(cells[0], 0, sizeof(cell_t) * grid_size * grid_size);
    }
    
    coordinate_t center_of_mass = calculate_overall_center_of_mass(particles, number_particles);
    printf("%.2f %.2f\n", particles->position.x, particles->position.y);
    printf("%.2f %.2f\n", center_of_mass.x, center_of_mass.y);

    // Free resources
    free(particles);
    free(cells[0]);
    free(cells);
    free(adjacent_cells[0][0]);
    free(adjacent_cells[0]);
    free(adjacent_cells);
    return 0;
}