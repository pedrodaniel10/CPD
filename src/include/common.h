#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>
#include <stdlib.h>

#define RND0_1 ((double) rand() / ((long long)1<<31))
#define G 6.67408e-11
#define EPSLON 0.01

typedef struct {
    double x;
    double y;
} coordinate_t;

typedef struct {
    coordinate_t position;
    coordinate_t velocity;
    double mass;
} particle_t;

typedef struct {
    double mass_sum;
    coordinate_t center_of_mass;
} cell_t;

void init_particles(long seed, long ncside, long long n_part, particle_t *par) {
    long long i;

    srand(seed);

    for(i = 0; i < n_part; i++) {
        par[i].position.x = RND0_1;
        par[i].position.y = RND0_1;
        par[i].velocity.x = RND0_1 / ncside / 10.0;
        par[i].velocity.y = RND0_1 / ncside / 10.0;

        par[i].mass = RND0_1 * ncside / (G * 1e6 * n_part);
    }
}

#endif