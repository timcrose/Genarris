#ifndef READ_INPUT_H
#define READ_INPUT_H

#include <stdio.h>
#include <stdlib.h>


typedef struct
{ 
	char *atoms;
	float *X;
	float *Y;
	float *Z;
	int num_of_atoms;
}molecule;

void read_control(int* num_structures, int* Z, float* Zp_max,
	float* volume_mean, float* volume_std, float *sr, long* max_attempts);

void read_geometry(molecule* mol);

void print_input_geometry(molecule* mol);

void print_molecule(molecule *mol);

void print_input_settings(int* num_structures, int* Z, float* Zp_max, 
	float* volume_mean, float* volume_std, float *sr, long *max_attempts);

#endif 
