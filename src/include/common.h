#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RND0_1 ((double) random() / ((long long)1<<31))
#define G 6.67408e-11
#define EPSLON 0.0005
#define ADJACENT_CELLS_NUMBER 9

typedef struct {
    double x;
    double y;
} coordinate_t;

typedef struct {
    int x;
    int y;
} coordinate_cell_t;

typedef struct {
    coordinate_cell_t cell;
    coordinate_t position;
    coordinate_t velocity;
    double mass;
} particle_t;

typedef struct {
    double mass_sum;
    coordinate_t center_of_mass;
} cell_t;

static coordinate_cell_t*** adjacent_cells;

particle_t * init_particles(long seed, long grid_size, long long number_particles) {
    long long i;
    particle_t *particles = (particle_t *) malloc (sizeof(particle_t) * number_particles);

    srandom(seed);

    for(i = 0; i < number_particles; i++) {
        particle_t *particle = &particles[i];
        
        particle->position.x = RND0_1;
        particle->position.y = RND0_1;
        particle->velocity.x = RND0_1 / grid_size / 10.0;
        particle->velocity.y = RND0_1 / grid_size / 10.0;

        particle->mass = RND0_1 * grid_size / (G * 1e6 * number_particles);

        particle->cell.x = particle->position.x * grid_size;
        particle->cell.y = particle->position.y * grid_size;
    }
    return particles;
}
#endif