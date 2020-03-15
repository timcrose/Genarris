#ifndef MOLECULE_UTILS_H
#define MOLECULE_UTILS_H

void recenter_molecule(molecule* mol);
void print_molecule(molecule *mol);
void molecule_rotate(molecule* mol, float rotation_matrix[3][3]);
void molecule_translate(molecule* mol, float translation[3]);
void rotation_matrix_from_vectors(float rotmat[3][3],
								  float a[3],
								  float b[3]);

void copy_positions_to_array(molecule *mol,
							 float *X,
							 float *Y,
							 float *Z       );
							 
void copy_positions_to_mol(molecule *mol,
						   float *X,
						   float *Y,
						   float *Z		 );
							 
  
#endif

