#ifndef STRUCTURES_H
#define STRUCTURES_H

#include "../include/common.h"

#define PARTICLE_SIZE 6
#define SIZEOF_PARTICLE(n) (n * PARTICLE_SIZE)
#define GET_NUMBER_PARTICLE(n) (n / PARTICLE_SIZE)

// Where n is the number of elements, p the number of processes and id is the rank of the process.
#define BLOCK_LOW(id, p, n) ((id) * (n) / (p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW((id) + 1, p, n) - 1)
#define BLOCK_SIZE(id, p, n) (BLOCK_HIGH(id, p, n) - BLOCK_LOW(id, p, n) + 1)
#define BLOCK_OWNER(index, p, n) (((p) * ((index) + 1) - 1) / (n))

typedef struct {
	double x;
	double y;
} coordinate_t;

typedef struct {
	int x;
	int y;
} coordinate_cell_t;

typedef struct {
	double index;
	coordinate_t position;
	coordinate_t velocity;
	double mass;
} particle_t;

typedef struct {
	double mass_sum;
	coordinate_t center_of_mass;
} cell_t;

#endif