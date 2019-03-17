#include "../include/common.h"
#include <omp.h>

void init_cells(long grid_size, cell_t** cells) {
    #pragma omp parallel for
    for(int i = 0; i < grid_size; i++) {
        for(int j = 0; j < grid_size; j++) {
            int index_adjacent_cells = 0;
            for (int x = -1; x <= 1; x++) {
                for(int y = -1; y <= 1; y++) {
                    cells[i][j].adjacent_cells[index_adjacent_cells].x = (i + (x + grid_size)) % grid_size;
                    cells[i][j].adjacent_cells[index_adjacent_cells].y = (j + (y + grid_size)) % grid_size;
                    index_adjacent_cells++;
                }
            }
        }        
    }
}

void calculate_centers_of_mass(particle_t *particles, cell_t **cells, int grid_size, int number_particles) {
    for (int i = 0; i < number_particles; i++) {
        particle_t *particle = &particles[i];
        cell_t *cell = &cells[particle->cell.x][particle->cell.y];

        cell->mass_sum += particle->mass;

        cell->center_of_mass.x += particle->mass * particle->position.x;
        cell->center_of_mass.y += particle->mass * particle->position.y;
    }
    
    #pragma omp parallel for if(grid_size >= 10)
    for (int i = 0; i < grid_size; i++) {
        for (int j = 0; j < grid_size; j++) {
            cell_t *cell = &cells[i][j];
            if (cell->mass_sum != 0) {
                cell->center_of_mass.x = cell->center_of_mass.x / cell->mass_sum;
                cell->center_of_mass.y = cell->center_of_mass.y / cell->mass_sum;
            }
        }
    }
}

void calculate_new_iteration(particle_t *particles, cell_t **cells, int grid_size, int number_particles) {
    #pragma omp parallel for
    for (int i = 0; i < number_particles; i++) {
        particle_t *particle = &particles[i];
        coordinate_t force, acceleration, velocity;
        cell_t cell_particle = cells[particle->cell.x][particle->cell.y];
        force.x = force.y = 0;
        acceleration.x = acceleration.y = 0;
        velocity.x = velocity.y = 0;

        // Calculate force
        for (int i = 0; i < ADJACENT_CELLS_NUMBER; i++) {
            coordinate_cell_t adjacent_cell = cell_particle.adjacent_cells[i];
            cell_t cell = cells[adjacent_cell.x][adjacent_cell.y];

            coordinate_t force_a_b;
            
            force_a_b.x = cell.center_of_mass.x - particle->position.x;
            force_a_b.y = cell.center_of_mass.y - particle->position.y;
            double distance = pow(pow(force_a_b.x, 2) + pow(force_a_b.y, 2), 1/2);

            if (distance < EPSLON) {
                double distance_cubed = pow(distance, 3);
                force_a_b.x *= G * particle->mass * cell.mass_sum / distance_cubed;
                force_a_b.y *= G * particle->mass * cell.mass_sum / distance_cubed;

                force.x += force_a_b.x;
                force.y += force_a_b.y;
            }  
        }

        // Calculate acceleration
        acceleration.x = force.x / particle->mass;
        acceleration.y = force.y / particle->mass;

        // Calculate new velocity
        particle->velocity.x += acceleration.x;
        particle->velocity.y += acceleration.y;
         
        // Calculate new position
        particle->position.x += particle->velocity.x + acceleration.x * 1/2;
        particle->position.y += particle->velocity.y + acceleration.y * 1/2;

        if (particle->position.x >= 1){
            particle->position.x -= 1;
        } else if (particle->position.x < 0) {
            particle->position.x += 1;
        }
        if (particle->position.y >= 1){
            particle->position.y -= 1;
        } else if (particle->position.y < 0) {
            particle->position.y += 1;
        }

        // Calculate new cell position
        particle->cell.x = particle->position.x * grid_size;
        particle->cell.y = particle->position.y * grid_size;
    }
}

coordinate_t calculate_overall_center_of_mass(particle_t* particles, int number_particles) {
    double center_of_mass_x, center_of_mass_y;
    double total_mass = center_of_mass_x = center_of_mass_y = 0;

    #pragma omp parallel for reduction(+: total_mass, center_of_mass_x, center_of_mass_y)
    for (int i = 0; i < number_particles; i++) {
        particle_t *particle = &particles[i];
        
        total_mass += particle->mass; 
        center_of_mass_x += particle->mass * particle->position.x;
        center_of_mass_y += particle->mass * particle->position.y;
    }
    center_of_mass_x /= total_mass;
    center_of_mass_y /= total_mass;

    return (coordinate_t) {center_of_mass_x, center_of_mass_y};
}

int main(int argc, const char** argv) {
    if (argc != 5) {
        printf("Expected 5 arguments, but %d were given\n", argc);
        exit(1);
    }

    int grid_size= atoi(argv[2]);
    int number_particles = atoi(argv[3]);
    int n_time_steps = atoi(argv[4]);

    // Allocate resources
    particle_t *particles = (particle_t *) malloc (sizeof(particle_t) * number_particles);
    cell_t **cells = (cell_t **) malloc(sizeof(cell_t *) * grid_size);
    
    for (int i = 0; i < grid_size; i++) {
        cells[i] = (cell_t *) calloc(sizeof (cell_t), grid_size);
    }
    
    init_particles(atoi(argv[1]), grid_size, number_particles, particles);
    init_cells(grid_size, cells);

    for (int n = 0; n < n_time_steps; n++) {
        calculate_centers_of_mass(particles, cells, grid_size, number_particles);
        calculate_new_iteration(particles, cells, grid_size, number_particles);

        for (int i = 0; i < grid_size; i++) {
            for(int j = 0; j < grid_size; j++) {
                cells[i][j].center_of_mass = (const coordinate_t) {0};
                cells[i][j].mass_sum = 0;
            }            
        }
    }
    
    coordinate_t center_of_mass = calculate_overall_center_of_mass(particles, number_particles);
    printf("%.2f %.2f\n", particles->position.x, particles->position.y);
    printf("%.2f %.2f\n", center_of_mass.x, center_of_mass.y);

    // Free resources
    free(particles);
    for(int i = 0; i < grid_size; i++) {
        free(cells[i]);
    }
    free(cells);

    return 0;
}