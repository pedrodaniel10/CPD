#include "../include/common.h"

#define DEBUG 0

void compute_centers_of_mass(particle_t *particles, cell_t **cells, int n_part, int grid_size) {
    int index_x;
    int index_y;
    double grid_unit = 1 / (double)grid_size;
    
    for (int i = 0; i < n_part; i++) {
        index_x = particles[i].position.x / grid_unit;
        index_y = particles[i].position.y / grid_unit;
        
        #if DEBUG
        printf("Particle %d\n", i);
        printf("Mass: %f\t", particles[i].mass);
        printf("Coordinates: (%f, %f)\t", particles[i].position.x, particles[i].position.y);
        printf("Index: (%d, %d)\n\n", index_x, index_y);
        #endif
        
        cells[index_x][index_y].mass_sum += particles[i].mass;

        cells[index_x][index_y].center_of_mass.x += particles[i].mass * particles[i].position.x;
        cells[index_x][index_y].center_of_mass.y += particles[i].mass * particles[i].position.y;
    }
    
    for (int i = 0; i < grid_size; i++) {
        for (int j = 0; j < grid_size; j++) {
            if (cells[i][j].mass_sum != 0) {
                cells[i][j].center_of_mass.x = cells[i][j].center_of_mass.x / cells[i][j].mass_sum;
                cells[i][j].center_of_mass.y = cells[i][j].center_of_mass.y / cells[i][j].mass_sum;
            }
        }
    }
}

int main(int argc, const char** argv) {
    if (argc != 5) {
        printf("Expected 5 arguments, but %d were given\n", argc);
        exit(1);
    }

    int grid_size= atoi(argv[2]);
    int n_part = atoi(argv[3]);
    //int n_time_steps = atoi(argv[4]);
    particle_t *particles = (particle_t *) malloc (sizeof(particle_t) * n_part);
    cell_t **cells = (cell_t **) malloc(sizeof(cell_t *) * grid_size);
    
    for (int i = 0; i < grid_size; i++) {
        cells[i] = (cell_t *) calloc(sizeof (cell_t), grid_size);
    }
    
    init_particles(atoi(argv[1]), grid_size, n_part, particles);
    compute_centers_of_mass(particles, cells, n_part, grid_size);

    #if DEBUG
    printf("----------------------------------------------------\n\n");
    for (int i = 0; i < grid_size; i++) {
        for (int j = 0; j < grid_size; j++) {
            printf("Cell [%d,%d]\n", i, j);
            printf("Mass sum: %f\t", cells[i][j].mass_sum);
            printf("Center of mass: (%f, %f)\n\n", cells[i][j].center_of_mass.x, cells[i][j].center_of_mass.y);
        }
    }
    #endif

    return 0;
}