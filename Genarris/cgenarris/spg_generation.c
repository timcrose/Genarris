#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include "read_input.h"
#include "spg_generation.h"
#include "combinatorics.h"
#include "lattice_generator.h"
#include "algebra.h"
#include "crystal_utils.h"
#include "molecule_utils.h"
#include "molecule_placement.h"
#include "spglib.h"

#define PI 3.141592653

extern unsigned int *seed2;
#pragma omp threadprivate(seed2)


int generate_crystal(crystal* random_crystal, molecule* mol,float volume,
	float Z, float Zp_max, int spg, COMPATIBLE_SPG compatible_spg[],
	int len_compatible_spg, int compatible_spg_index)
{
	
	//crystal random_crystal;
	float max_angle = 30 * PI/180;
	float min_angle = 150 * PI/180;

	random_crystal->Z = Z;
	int N = mol->num_of_atoms;

	//copy molecules to an array to save it. molecule might deform
	//upon many rotations.	
	float Xm[N]; //molecule X coordinate
	float Ym[N];
	float Zm[N];
	copy_positions_to_array(mol, Xm, Ym, Zm);

	int hall_number;
	hall_number = hall_number_from_spg(spg);
	
	generate_lattice(random_crystal->lattice_vectors, spg, max_angle, min_angle, volume);
	
	//find a random pos
	int pos_index = rand_r(seed2) % compatible_spg[compatible_spg_index].num_allowed_pos;
	int pos = compatible_spg[compatible_spg_index].allowed_pos[pos_index];
	random_crystal->wyckoff_position = pos;
	
	//place, align and attempt to generate crystal at position pos
	int result = auto_align_and_generate_at_position(random_crystal,
							mol,
							hall_number,
							spg, 
							pos_index,
							compatible_spg[compatible_spg_index]);
	//copy back to mol
	copy_positions_to_mol(mol, Xm, Ym, Zm);
	
	if(!result)
	{
		
		return 0;
	}
	else
		return 1;
}

