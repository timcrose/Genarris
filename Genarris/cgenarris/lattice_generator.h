#ifndef LATTICE_GENERATOR_H
#define LATTICE_GENERATOR_H

void gen_triclinic_lattice(float lattice_vector[3][3], 
	float target_volume, float max_angle, float min_angle);
void gen_monoclinic_lattice(float lattice_vector[3][3], 
	float target_volume, float max_angle, float min_angle);
void gen_orthorhombic_lattice(float lattice_vector[3][3], 
	float target_volume);
void gen_tetragonal_lattice(float lattice_vector[3][3],
	float target_volume);
void gen_hexagonal_lattice(float lattice_vector[3][3], 
	float target_volume);
void gen_cubic_lattice(float lattice_vector[3][3],
	float target_volume);
void generate_lattice(float lattice_vector[3][3], int spg,
 float max_angle, float min_angle, float target_volume);
void generate_fake_lattice(float lattice_vector[3][3], int spg);

void standardise_lattice(float lattice[3][3], int spg);

#endif

