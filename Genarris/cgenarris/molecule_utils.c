#include <stdio.h>
#include <math.h>
#include "read_input.h"
#include "algebra.h"

extern float TOL;
//function to recenter the molecule to origin
//TODO : compute and align the molecule with principal axis
void recenter_molecule(molecule* mol)
{
	int N = (*mol).num_of_atoms;
	float xcom = 0, ycom = 0 , zcom = 0;
	for(int i = 0; i < N; i++)
	{
		xcom += (*mol).X[i];
		ycom += (*mol).Y[i];
		zcom += (*mol).Z[i];
	}
	
	xcom /= N;
	ycom /= N;
	zcom /= N;
	
	//printf("xcom = %f, ycom = %f, zcom =%f \n", xcom, ycom, zcom);
	
	for(int i = 0; i < N; i++)
	{
		(*mol).X[i] -= xcom;
		(*mol).Y[i] -= ycom;
		(*mol).Z[i] -= zcom;
	}
	
}


//print a molecule information contained in structure
void print_molecule(molecule *mol)
{
	int N = (*mol).num_of_atoms;
	printf("#total atoms in the molecule = %d \n\n", N);
	
	for(int i =0; i < N; i++)
	{
		printf("atom %12f %12f %12f %4c%c \n",
			(*mol).X[i], (*mol).Y[i], (*mol).Z[i], (*mol).atoms[2*i],(*mol).atoms[2*i+1]);
	}
}

//to rotate a molecule 
void molecule_rotate(molecule* mol, float rotation_matrix[3][3])
{
	int N = (*mol).num_of_atoms;
	for(int i = 0; i < N; i++)
	{
		float temp[3] = {(*mol).X[i], (*mol).Y[i], (*mol).Z[i]};
		vector3_mat3b3_multiply(rotation_matrix, temp, temp);
		(*mol).X[i] = temp[0];
		(*mol).Y[i] = temp[1];
		(*mol).Z[i] = temp[2];
	}
}

void molecule_translate(molecule* mol, float translation[3])
{
	int N = (*mol).num_of_atoms;
	for(int i = 0; i < N; i++)
	{
		(*mol).X[i] += translation[0];
		(*mol).Y[i] += translation[1];
		(*mol).Z[i] += translation[2];
	}
}

void rotation_matrix_from_vectors(float rotmat[3][3],
								  float a[3],
								  float b[3])
{
	//find axis and angle
	normalise_vector3(a);
	normalise_vector3(b);
	float axis[3];
	cross_vector3_vector3(axis,a,b);
	float angle;
	angle = acos(dot_vector3_vector3(a,b));
	
	//if and b are the same, return identity
	if ( check_vec3_isNull(axis, TOL) )
	{
		rotmat[0][0] = 1;
		rotmat[0][1] = 0;
		rotmat[0][2] = 0;
	
		rotmat[1][0] = 0;
		rotmat[1][1] = 1;
		rotmat[1][2] = 0;
		
		rotmat[2][0] = 0;
		rotmat[2][1] = 0;
		rotmat[2][2] = 1;
		
		return;
	}
	//find rotation matrix from axis and angle
	normalise_vector3(axis);
	rotation_mat_around_axis(rotmat, axis, angle);
	
	//print_mat3b3(rotmat);
}

void copy_positions_to_array(molecule *mol,
							 float *X,
							 float *Y,
							 float *Z       )
{
	int N = mol->num_of_atoms;
	for(int i = 0; i < N; i++)
	{
		X[i] = mol->X[i];
		Y[i] = mol->Y[i];
		Z[i] = mol->Z[i];
	}
}

void copy_positions_to_mol(molecule *mol,
						   float *X,
						   float *Y,
						   float *Z		 )
{
	int N = mol->num_of_atoms;
	for(int i = 0; i < N; i++)
	{
		mol->X[i] = X[i];
		mol->Y[i] = Y[i];
		mol->Z[i] = Z[i];
	}
}
