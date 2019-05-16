#ifndef ARRAY_LIST_H
#define ARRAY_LIST_H

#include "structures.h"

enum { DIAGONAL_UP_LEFT, UP, DIAGONAL_UP_RIGHT, LEFT, RIGHT, DIAGONAL_DOWN_LEFT, DOWN, DIAGONAL_DOWN_RIGHT };

typedef struct {
	int length;
	int size;
	particle_t* particles;
} array_list_t;

typedef struct {
	int rank;
	array_list_t* particles_buffer_send;
	array_list_t* particles_buffer_recv;
	int length_send_buffer;
	cell_t* cells_buffer_send;
	cell_t* cells_buffer_recv;
	int sent;
	int received;
	int index;
} node_t;

array_list_t* create_array_list(int initial_size) {
	if (initial_size <= 0) {
		return NULL;
	}

	array_list_t* list = (array_list_t*)malloc(sizeof(array_list_t));
	list->length = 0;
	list->size = initial_size;
	list->particles = (particle_t*)malloc(sizeof(particle_t) * initial_size);

	return list;
}

void append(array_list_t* list, particle_t particle) {
	if (list->size == list->length) {
		list->size *= 2;
		list->particles = (particle_t*)realloc(list->particles, sizeof(particle_t) * list->size);
	}
	list->particles[list->length] = particle;
	list->length++;
}

particle_t* allocate_array(array_list_t* list, int size) {
	while (list->size - list->length < size) {
		list->size *= 2;
		list->particles = (particle_t*)realloc(list->particles, sizeof(particle_t) * list->size);
	}
	particle_t* returnPointer = &list->particles[list->length];
	list->length += size;
	return returnPointer;
}

particle_t* list_get(array_list_t* list, int index) {
	if (index >= list->length) {
		return NULL;
	}
	return &list->particles[index];
}
#endif