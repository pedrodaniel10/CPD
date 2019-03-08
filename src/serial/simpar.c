#include "../include/common.h"

void compute_centers_of_mass(particle_t *particles, cell_t **cells, int n_part, int grid_size) {
    int index_x;
    int index_y;
    double grid_size_d = grid_size;
    double grid_unit = 1/grid_size_d;
    
    for (int i = 0; i < n_part; i++) {
        index_x = particles[i].position.x/grid_unit;
        index_y = particles[i].position.y/grid_unit;
        
        /*
        printf("Particle %d\n", i);
        printf("Mass: %f\t", particles[i].mass);
        printf("Coordinates: (%f, %f)\t", particles[i].position.x, particles[i].position.y);
        printf("Index: (%d, %d)\n\n", index_x, index_y);
        */
        
        cells[index_x][index_y].mass_sum += particles[i].mass;

        cells[index_x][index_y].center_of_mass.x += particles[i].mass*particles[i].position.x;
        cells[index_x][index_y].center_of_mass.y += particles[i].mass*particles[i].position.y;
    }
    
    for (int i = 0; i < grid_size; i++) {
        for (int j = 0; j < grid_size; j++) {
            cells[i][j].center_of_mass.x = cells[i][j].center_of_mass.x/cells[i][j].mass_sum;
            cells[i][j].center_of_mass.y = cells[i][j].center_of_mass.y/cells[i][j].mass_sum;
        }
    }
}

int main(int argc, const char** argv) {
    int seed;
    int grid_size;
    int n_part;
    //int n_time_steps;
    particle_t *particles;
    cell_t **cells;
    

    seed = atoi(argv[1]);
    grid_size = atoi(argv[2]);
    n_part = atoi(argv[3]);
    //n_time_steps = atoi(argv[4]);

    particles = (particle_t *) malloc (sizeof(particle_t) * n_part);
    
    cells = (cell_t **) malloc(sizeof(cell_t *) * grid_size);
    for (int i = 0; i < grid_size; i++) {
        cells[i] = (cell_t *) malloc(sizeof (cell_t) * grid_size);
    }
    for (int i = 0; i < grid_size; i++) {
        for (int j = 0; j < grid_size; j++) {
            cells[i][j].mass_sum = 0;
            cells[i][j].center_of_mass.x = 0;
            cells[i][j].center_of_mass.y = 0;
        }
    }
    
    init_particles(seed, grid_size, n_part, particles);
    compute_centers_of_mass(particles, cells, n_part, grid_size);

    /*
    printf("----------------------------------------------------\n\n");
    for (int i = 0; i < grid_size; i++) {
        for (int j = 0; j < grid_size; j++) {
            printf("Cell [%d,%d]\n", i, j);
            printf("Mass sum: %f\t", cells[i][j].mass_sum);
            printf("Center of mass: (%f, %f)\n\n", cells[i][j].center_of_mass.x, cells[i][j].center_of_mass.y);
        }
    }
    */
    return 0;
}