#ifndef PARTICLES_H
#define PARTICLES_H
#include "common.h"

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

static cell_t ****adjacent_cells;

particle_t *init_particles(long seed, long grid_size, long long number_particles) {
	long long i;
	particle_t *particles = (particle_t *)malloc(sizeof(particle_t) * number_particles);

	srandom(seed);

	for (i = 0; i < number_particles; i++) {
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