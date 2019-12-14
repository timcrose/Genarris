#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "read_input.h"
#include "spg_generation.h"
#include "algebra.h"
#include "spglib.h"
#include "combinatorics.h"
#include "randomgen.h"
#include "molecule_placement.h"
#include "crystal_utils.h"
#include "check_structure.h"
#include "molecule_utils.h"
#include "point_group.h"
#include "lattice_generator.h"

extern float TOL;
#define ZMAX 192
#define CONST_TOL 0.001

float viewing_directions[16][3] = {
			{1, 0, 0}, 
			{0, 1, 0},
			{0, 0, 1},
			{1, -1, 0},
			{1, 1, 0},
			{-1, -1, 0},
			{1, 2, 0},
			{-2, -1, 0},
			{1, 1, 1},
			{0, 1, -1},
			{-1, 0, 1},
			{1, -1, -1},
			{-1, 1, -1},
			{-1, -1, 1},
			{0, 1, 1},
			{1, 0, 1}
								};
#pragma omp threadprivate(viewing_directions)

const int num_viewing_direction = 16;



/*
 * TODO: finding compatible spg and positions:
 * 1)choose an appropriate lattice consistent with the space group
 * 2)do the second degree alignment. find all possible rotations
 * about the primary axis to see if special position is satisfied. 
 * 3)find overlapping molecules without actually placing the molecule
 * 
 * improve speed: replace pdist with normal dist after moving molecules 
 * to origin

*/


/*
 * Find compatible general positions. doesnot care about symmetry
 */

void find_allowed_spg(int allowed_spg[230], int* num_allowed_spg, int Z)
{
	//just for calling spglib; not used
	double translations[192][3];
	int rotations[192][3][3];
	*num_allowed_spg = 0;
	int counter = 0;
	//loop over space groups
	for(int i = 0; i < 230; i++)
	{
		int	hall_number; 
		hall_number = hall_number_from_spg(i+1);
		int num_of_operations = spg_get_symmetry_from_database(rotations,
			translations, hall_number);
		if (num_of_operations == Z)
		{
			allowed_spg[counter] = i+1;
			counter++;
		}
	}
	
	*num_allowed_spg = counter;
}

/* find compatible positions once the symmetry is given
 */

void find_allowed_positions_using_molecular_symmetry(char mol_sym[6],
									int Z,
									int Zpp)
{
	
	int index = 0;
	int spg_counter = 0;
	int position_counter = 0;
	printf("molecular symmetry = %s\n", mol_sym);
	while (index < 31)
	{
		if (strcmp(mol_sym, point_group_dbs[index].point_group) == 0)
		{
			break;
		}
		index++;
	}
	// index has the location of point group
	
	for (int spg = 0; spg < 230; spg++ )
	{
		int max = spg_positions[spg].num_of_positions;
		int temp = position_counter;
		for(int pos = 0; pos < max; pos++ )
		{
			if (spg_positions[spg].multiplicity[pos] == Z)
			{
				for(int i = 0; i < point_group_dbs[index].len_sub_groups; i++)
				{
					if (strcmp(spg_positions[spg].site_symmetry[pos], point_group_dbs[index].sub_groups[i]) == 0)
						{	printf("spg:%d wyckoff position:%d%c site symmetry:%s allowed \n", spg+1,
							spg_positions[spg].multiplicity[pos],
							spg_positions[spg].wyckoff_letter[pos],
							spg_positions[spg].site_symmetry[pos]);
							position_counter++;  }
				}
			}
		}
		
		if(position_counter > temp)
			{spg_counter++;}	
	}
	
	printf("Total allowed spacegroup types : %d \nTotal allowed positions: %d\n", spg_counter, position_counter);
	
}

/*find all compatible space groups and positions without any knowledge 
 * of molecular symmetry
 */
void find_compatible_spg_positions(molecule *mol,
									int Z,
									COMPATIBLE_SPG compatible_spg[],
									int *num_compatible_spg,
									int thread_num)
{
	crystal xtal;
	int N = mol->num_of_atoms;
	xtal.Xcord = (float *)malloc(ZMAX*N*sizeof(float));
	xtal.Ycord = (float *)malloc(ZMAX*N*sizeof(float));
	xtal.Zcord = (float *)malloc(ZMAX*N*sizeof(float));
	xtal.atoms = (char *)malloc(ZMAX*N*2*sizeof(char));
	xtal.num_atoms_in_molecule = mol->num_of_atoms;
	
	//for storing information of allowed position and spg
	unsigned int pos_list[47];
	unsigned int len_pos_list = 0;
	
	//possibly equivalent atoms
	int eq_atoms[N];
	int len_eq_atoms = 0 ;
	find_farthest_equivalent_atoms(mol, eq_atoms, &len_eq_atoms);
	

	//allocate memory for mol_axes
	float *mol_axes = (float *)malloc(3*(len_eq_atoms*len_eq_atoms+3)*sizeof(float));
	int num_axes = 0;	
	find_possible_mol_axis(mol, mol_axes, &num_axes, eq_atoms, len_eq_atoms);
	//for (int i = 0; i < num_axes; i++)
	//	printf("%f %f %f \n", *((*mol_axes)+3*i+0), *((*mol_axes)+3*i+1), *((*mol_axes)+3*i+2) );
	
	for (int spg = 0; spg < 230; spg++ )
	{
		//make a fake lattice with large volume
		//generate_lattice(xtal.lattice_vectors, spg+1, 60, 120, 1000);
		generate_fake_lattice(xtal.lattice_vectors, spg+1);
		int max = spg_positions[spg].num_of_positions; 
		int hall_number = hall_number_from_spg(spg+1);
		for(int pos = 0; pos < max; pos++ )
		{		
			//check multiplicity
			if (spg_positions[spg].multiplicity[pos] != Z)
				{continue;}
			
			//general position
			if (pos == 0)
			{	
				if(thread_num == 1)
					printf("spg %3d Wyckoff position %d%c with site symmetry %3s is compatible.\n",
						spg+1,
						spg_positions[spg].multiplicity[pos],
						spg_positions[spg].wyckoff_letter[pos],
						spg_positions[spg].site_symmetry[pos]);
						
				pos_list[len_pos_list] = pos;
				compatible_spg[*num_compatible_spg].pos_overlap_list[len_pos_list] = (int *)\
				malloc(sizeof(int));
		
				*(compatible_spg[*num_compatible_spg].pos_overlap_list[len_pos_list])\
					= 0;
				
				
				len_pos_list++;
				continue;
			}
			
			
			//get rot and trans from database
			int rot[3][3];
			float trans[3];
			for(int i = 0; i < 3; i++)
			{
				trans[i] = spg_positions[spg].first_position_trans[pos][i];
				for(int j = 0; j < 3; j++)
				{
					rot[i][j] = spg_positions[spg].first_position_rot[pos][i][j];
				}
			}

			float rand_frac_array[3]={0.19316, 0.74578, 0.26245};
			vector3_intmat3b3_multiply(rot,rand_frac_array,
									rand_frac_array);
			vector3_add(trans,
						rand_frac_array,
						rand_frac_array);
			
			//for converting to fractional coordinates and
			// place the first molecule
			float inverse_lattice_vectors[3][3];
			inverse_mat3b3(inverse_lattice_vectors,
						xtal.lattice_vectors);
			mat3b3_transpose(inverse_lattice_vectors,
							inverse_lattice_vectors);
			float mol_Xfrac[N]; //stores fractional ...
			float mol_Yfrac[N]; //coordinates of first mol
			float mol_Zfrac[N];
			
			for(int i = 0; i < N; i++)
			{
				float atom_i[3] = {(*mol).X[i],
								(*mol).Y[i],
								(*mol).Z[i]};
				vector3_mat3b3_multiply(inverse_lattice_vectors,
										atom_i,
										atom_i);
				vector3_add(atom_i, rand_frac_array, atom_i); 
				mol_Xfrac[i] = atom_i[0];
				mol_Yfrac[i] = atom_i[1];
				mol_Zfrac[i] = atom_i[2];
			} 
			
			//apply symm operations
			apply_all_symmetry_ops(&xtal,
								mol,
								mol_Xfrac,
								mol_Yfrac,
								mol_Zfrac,
								N,
								hall_number);
			xtal.spg = spg+1;
			
			/*find the molecules which overlap with the first 
			* molecule com which is rand_frac_array.
			* Check ifits length compatible with the multiplicity of 
			* the special position.
			* find the atoms farthest from com
			* average over pairs of atoms of different molecule
			* test if molecule overlaps again after symm operations
			*/
			float lattice_vectors_trans[3][3];
			mat3b3_transpose(lattice_vectors_trans,
							xtal.lattice_vectors);
			
			//multiplicity of general position 
			int Z_gen = spg_positions[spg].multiplicity[0];
			int order = Z_gen / Z;  
			int overlap_list[order];
			int len_overlap_list = 0;
			//com of first molecule in cartesian
			float first_com[3];
			vector3_mat3b3_multiply(lattice_vectors_trans,
									rand_frac_array,
									first_com);
			//transpose again to get back
			mat3b3_transpose(inverse_lattice_vectors,\
							inverse_lattice_vectors);
							
			//find overlap list. The overlap list is not specfic to the 
			//molecule. Hence, techically it needs to be calculated only 
			//once and can be saved. However, it is calculated everytime
			//cgenarris executes on the fly.
			find_overlap_list (xtal,
								first_com,
								inverse_lattice_vectors,
								overlap_list,
								&len_overlap_list,
								Z_gen, 
								N);

			// if the overlap isn't the multiplicity, its a problem
			if(len_overlap_list != order)
			{	
				if(thread_num == 1)
				printf("***ERROR: spg = %d, pos = %d overlap error."
					"Detected %d order, order should be %d\n",
					spg+1, pos, len_overlap_list, order);
				//print_crystal(&xtal);	
				continue;
			}
			
			//bring all molecules to the origin to check overlap of molecules
			//using the list of overlap molecules
			bring_molecules_to_origin(&xtal);

			//debug 
			compatible_spg[*num_compatible_spg].compatible_axes[pos].num_combinations = 0;
			//check compatiblility of a wyckoff position using standard 
			//viewing directions
			int result = check_pos_compatibility_using_std_orientations(&xtal,
									&compatible_spg[*num_compatible_spg].compatible_axes[pos],
									mol,
									hall_number,
									mol_axes,
									num_axes,
									rand_frac_array,
									overlap_list,
									len_overlap_list);
					
			//not compatible? then skip
			if( !result)
				{continue;}
			
			//if compatible
			if(thread_num == 1)
				printf("spg %3d Wyckoff position %d%c with site symmetry %3s is compatible.\n",
					spg+1,
					spg_positions[spg].multiplicity[pos],
					spg_positions[spg].wyckoff_letter[pos],
					spg_positions[spg].site_symmetry[pos]);
				
			pos_list[len_pos_list] = pos;
			compatible_spg[*num_compatible_spg].pos_overlap_list[len_pos_list] = (int *)\
				malloc(len_overlap_list*sizeof(int) );
				
			for(int i =0; i < len_overlap_list; i++)
			{
				*(compatible_spg[*num_compatible_spg].pos_overlap_list[len_pos_list]+i)\
					= overlap_list[i];
			}
				
			len_pos_list++;
		}//end of pos loop
	
		// pos_list contains information about all compatible wyckoff 
		//positions. len_pos_list is its length. Now allocate memory to 
		// allowed pos inside the compatible spg structure and copy 
		//the pos_list for this space group. 
		//The number of allowed pos is also stored in the structure.
		if(len_pos_list != 0)
		{
			compatible_spg[*num_compatible_spg].allowed_pos =\
			(unsigned int *)malloc(len_pos_list*(sizeof(unsigned int)));
			compatible_spg[*num_compatible_spg].spg = spg+1;
			compatible_spg[*num_compatible_spg].num_allowed_pos = len_pos_list;
			for (int i = 0; i < len_pos_list; i++)
			{
				compatible_spg[*num_compatible_spg].allowed_pos[i] = pos_list[i];
				pos_list[i] = 0;
			}  
			(*num_compatible_spg)++;
		}
		
		len_pos_list = 0;
	}//end of spg loop
	
	
}



int check_pos_compatibility_using_std_orientations(crystal* xtal_1,	
							  COMPATIBLE_AXES *comp_axes,
							  molecule* mol,		
							  int hall_number,		
							  float *mol_axes,	
							  int num_axes,	
							  float first_com[3],
							  int overlap_list[],
							  int len_overlap_list)
{
	int N = mol->num_of_atoms;
	float rotation_matrix[3][3];
	float mol_Xfrac[N]; //stores fractional ...
	float mol_Yfrac[N]; //coordinates of first mol
	float mol_Zfrac[N];
	float lattice_vectors[3][3];
	float inv_lattice_vectors[3][3];
	copy_mat3b3_mat3b3(lattice_vectors, xtal_1->lattice_vectors);
	inverse_mat3b3(inv_lattice_vectors, lattice_vectors);
	mat3b3_transpose(inv_lattice_vectors, inv_lattice_vectors);
	
	//create a temp xtal
	crystal *xtal = (crystal*)malloc(sizeof(crystal));
	copy_mat3b3_mat3b3(xtal->lattice_vectors, xtal_1->lattice_vectors);
	allocate_xtal(xtal, ZMAX, N);
	xtal->num_atoms_in_molecule = mol->num_of_atoms;
	
	//debug
	int spg = xtal_1->spg;

	(*comp_axes).usable_mol_axes = NULL;
	(*comp_axes).usable_view_dir = NULL;
	
	
	//pick a mol axis
	for (int i = 0; i < num_axes; i++)
	{
		float mol_axis[3] = {*(mol_axes+3*i + 0),
							 *(mol_axes+3*i + 1),
							 *(mol_axes+3*i + 2) };
		//pick a viewing direction 
		for (int j = 0; j < num_viewing_direction; j++)
		{
			float view_dir[3] = {viewing_directions[j][0],
								 viewing_directions[j][1],
								 viewing_directions[j][2] };
			normalise_vector3(view_dir);
						
			//create rotation matrix
			rotation_matrix_from_vectors(rotation_matrix,mol_axis, view_dir);
			//rotate and store first molecule in fractional coord 
			//in mol_frac
			for(int z = 0; z < N; z++)
			{
				float temp[3] = {xtal_1->Xcord[z] ,
								 xtal_1->Ycord[z] ,
								 xtal_1->Zcord[z] };
				vector3_mat3b3_multiply(rotation_matrix, temp, temp);
				//vector3_add(temp,com,temp);
				//convert to frac
				vector3_mat3b3_multiply(inv_lattice_vectors,
										temp,
										temp );
				vector3_add(temp, first_com, temp);
				mol_Xfrac[z] = temp[0];
				mol_Yfrac[z] = temp[1];
				mol_Zfrac[z] = temp[2];
			}
			
			apply_all_symmetry_ops(xtal,
						mol,
						mol_Xfrac,
						mol_Yfrac,
						mol_Zfrac,
						N,
						hall_number);
			
			//bring_molecules_to_origin(xtal);

			int result = check_overlap_xtal(xtal,
										overlap_list,
										len_overlap_list,
										N);
					
			if (result)
			{
				int Z_gen = xtal->Z;
				combine_close_molecules(xtal);
				int Z_return = xtal->Z;
				//len_overlap_list == order, Z/len = multiplicity
				if(Z_return != Z_gen/len_overlap_list)
					{continue;}

				int comb = (*comp_axes).num_combinations;
				
				(*comp_axes).usable_mol_axes = (float *) 
						realloc((*comp_axes).usable_mol_axes, 3*(comb+1)*sizeof(float));
				(*comp_axes).usable_view_dir = (float *) 
						realloc((*comp_axes).usable_view_dir, 3*(comb+1)*sizeof(float));

				*((*comp_axes).usable_mol_axes + 3*comb + 0) = mol_axis[0];
				*((*comp_axes).usable_mol_axes + 3*comb + 1) = mol_axis[1];
				*((*comp_axes).usable_mol_axes + 3*comb + 2) = mol_axis[2];
				
				*((*comp_axes).usable_view_dir + 3*comb + 0) = view_dir[0];
				*((*comp_axes).usable_view_dir + 3*comb + 1) = view_dir[1];
				*((*comp_axes).usable_view_dir + 3*comb + 2) = view_dir[2];
				
				(*comp_axes).num_combinations++;
			}
		}
	}
	
	free_xtal(xtal);
	if((*comp_axes).num_combinations)
		return 1;
	else
		return 0;
}


//find all possible axes of the single molecule. return an array after 
//removing duplicates
void find_possible_mol_axis(molecule *mol, float *mol_axes, int* num_axes,\
							int eq_atoms[], int len_eq_atoms)
{
	//add principle axes
	*num_axes = 0;
	/*
	float p1[3] = {1, 0, 0};
	float p2[3] = {0, 1, 0};
	float p3[3] = {0, 0, 1};
	
	add_axis_to_mol_axes(mol_axes, num_axes, p1);
	add_axis_to_mol_axes(mol_axes, num_axes, p2);
	add_axis_to_mol_axes(mol_axes, num_axes, p3);
	*/
	int N = mol->num_of_atoms;
	for(int i = 0; i < len_eq_atoms; i++)
	for(int j = i; j < len_eq_atoms; j++)
	{
		int atomi_index = eq_atoms[i];
		int atomj_index = eq_atoms[j];
		
		float axis_i[3]= {mol->X[atomi_index], mol->Y[atomi_index], mol->Z[atomi_index]};
		float axis_j[3]= {mol->X[atomj_index], mol->Y[atomj_index], mol->Z[atomj_index]};
		float avg[3];
		float cross[3];
				
		if( !compute_average_axis(axis_i, axis_j, avg) )
			{avg[0] = 1; avg[1] = 0; avg[2] = 0;}
		if( !compute_cross_axis(axis_i, axis_j, cross) )
			{cross[0] = 1; cross[1] = 0; cross[2] = 0;} 
		
		//assume max memory of N^2 for axes
		add_axis_to_mol_axes(mol_axes, num_axes, avg);
		//do the same for cross
		add_axis_to_mol_axes(mol_axes, num_axes, cross);
				
	}
}


void add_axis_to_mol_axes(float *mol_axes, int *num_axes, float axis[3])
{
	int check = 1;
	for(int z = 0; z < *num_axes; z++ )
	{
		float temp_axis[3] = {*(mol_axes+3*z + 0) ,
							  *(mol_axes+3*z + 1) ,
							  *(mol_axes+3*z + 2) };
		float is_parallel[3];
		cross_vector3_vector3(is_parallel, temp_axis, axis);
		if (check_vec3_isNull(is_parallel, TOL) )
		{
			check = 0;
			break;
		}
	}
		
	//if axis is not in list then add
	if (check)
	{
		*(mol_axes + *num_axes*3 + 0) = axis[0];
		*(mol_axes + *num_axes*3 + 1) = axis[1];
		*(mol_axes + *num_axes*3 + 2) = axis[2];
		(*num_axes)++;
	}
}






void find_farthest_equivalent_atoms(molecule *mol, 	 	\
									int atom_index[],	\
									int* len_atom_index	)
{
	
	//find farthest atom
	int N = mol->num_of_atoms;
	float max = 0;
	int max_atom_index = 0;
	
	recenter_molecule(mol);
	
	for (int i = 0; i < N; i++)
	{
		float atom_dist[3] = {mol->X[i], mol->Y[i], mol->Z[i]};
		if (vector3_norm(atom_dist) > max)
		{
			max = vector3_norm(atom_dist);
			max_atom_index = i ;
		}
	}
	
	int count = 0; 
	//find all equivalent atoms
	for (int i = 0; i < N; i++)
	{
		float atom_dist[3] = {mol->X[i], mol->Y[i], mol->Z[i]};
		if ( fabs (vector3_norm(atom_dist) - max ) < TOL && \
			 mol->atoms[2*i] == mol->atoms[2*max_atom_index] )
			{
				atom_index[count] = i;
				count++;
			}
	}
	*len_atom_index = count;
}




int check_overlap_xtal(crystal* xtal,
					   int overlap_list[],
					   int len_overlap_list,
					   int N)
{
	float lattice_vectors[3][3];
	float inv_lattice_vectors[3][3];
	copy_mat3b3_mat3b3(lattice_vectors, xtal->lattice_vectors);
	inverse_mat3b3(inv_lattice_vectors, lattice_vectors);
	
	for (int i = 1; i < len_overlap_list; i++)
	{
		int mol_index = overlap_list[i];
		for(int j = 0; j < N; j++)
		{
			float atom1[3] = {xtal->Xcord[mol_index+j], \
							  xtal->Ycord[mol_index+j], \
							  xtal->Zcord[mol_index+j]  };
							  
			float min_dist = 1000*TOL ;
			
			for(int k = 0; k < N; k++)
			{
				float atom2[3] = { xtal->Xcord[k], \
								   xtal->Ycord[k], \
								   xtal->Zcord[k]	};
				float dist = pdist_appx(lattice_vectors,		\
								   inv_lattice_vectors,	\
								   atom1[0],
								   atom1[1],
								   atom1[2],
								   atom2[0],
								   atom2[1],
								   atom2[2]				);
				
				//float dist = cart_dist(atom1, atom2);
				if (dist < min_dist)
					{min_dist = dist ;}
			}
			
			if (min_dist > 1*TOL)
				{return 0;}
		}
	}
	
	return 1;
}

/*
*function to check overlaps in a crystal using simple cartesian distance
* uses the overlap list to know location of symmetrically equivalent 
* molecules
*/
int check_overlap_xtal_cartesian(crystal* xtal,
								int overlap_list[],
								int len_overlap_list,
								int N)
{
	float lattice_vectors[3][3];
	float inv_lattice_vectors[3][3];
	copy_mat3b3_mat3b3(lattice_vectors, xtal->lattice_vectors);
	inverse_mat3b3(inv_lattice_vectors, lattice_vectors);
	
	float first_com[3], second_com[3];
	compute_molecule_COM(*xtal, first_com, 0);
	
	for (int i = 1; i < len_overlap_list; i++)
	{
		int mol_index = overlap_list[i];
		compute_molecule_COM(*xtal, second_com, mol_index);
		
		if ( (first_com[0] - second_com[0])*(first_com[0] - second_com[0]) +
			 (first_com[1] - second_com[1])*(first_com[1] - second_com[1]) +
			 (first_com[2] - second_com[2])*(first_com[2] - second_com[2])		
								> TOL 										)
		{
			float displace[3];
			vector3_subtract(second_com, first_com, displace);
			
			for(int l = mol_index; l < mol_index + N; l++)
			{
				xtal->Xcord[l] -= displace[0];
				xtal->Ycord[l] -= displace[1];
				xtal->Zcord[l] -= displace[2];  
			}
			
		}
		
		for(int j = 0; j < N; j++)
		{
			float atom1[3] = {xtal->Xcord[mol_index+j], 
							  xtal->Ycord[mol_index+j], 
							  xtal->Zcord[mol_index+j]  };
							  
			float min_dist = 1000*TOL ;
			
			for(int k = 0; k < N; k++)
			{
				float atom2[3] = { xtal->Xcord[k], 
								   xtal->Ycord[k], 
								   xtal->Zcord[k]	};
				float dist = (atom1[0] - atom2[0]) * (atom1[0] - atom2[0])+
							 (atom1[1] - atom2[1]) * (atom1[1] - atom2[1])+
							 (atom1[2] - atom2[2]) * (atom1[2] - atom2[2]);

				//float dist = cart_dist(atom1, atom2);
				if (dist < min_dist)
					{min_dist = dist ;}
			}
			
			if ( sqrt(min_dist) > 1*TOL)
				{return 0;}
		}
	}
	
	return 1;
}



/*compute avg of two atoms in a molecule. return 0 if the resultant is 
 * too small
 */
int compute_atom_average_xtal(crystal* xtal, int index1, int index2, float avg[3])
{
	//safety
	int total = xtal->num_atoms_in_molecule * xtal->Z;
	if (index1 > total || index2 > total)
		{printf("compute_atom_avg: out of bounds %d, %d \n", index1, index2), print_crystal(xtal); exit(0);}
	
	float atom1[3] = {xtal->Xcord[index1], \
					  xtal->Ycord[index1], \
					  xtal->Zcord[index1]  };
	float atom2[3] = {xtal->Xcord[index2], \
					  xtal->Ycord[index2], \
					  xtal->Zcord[index2]  };
	
	normalise_vector3(atom1);
	normalise_vector3(atom2);

	avg[0] = (atom1[0]+atom2[0])/2 ;
	avg[1] = (atom1[1]+atom2[1])/2 ;
	avg[2] = (atom1[2]+atom2[2])/2 ;
	
	if( check_vec3_isNull(avg, TOL) )
		{return 0;}
	
	normalise_vector3(avg);
	return 1;
	
}

int compute_atom_cross_xtal(crystal* xtal, int index1, int index2, float avg[3])
{
	//safety
	int total = xtal->num_atoms_in_molecule * xtal->Z;
	if (index1 > total || index2 > total)
		{printf("compute_atom_cross: out of bounds, %d, %d \n", index1, index2); exit(0);}
	
	float atom1[3] = {xtal->Xcord[index1], \
					  xtal->Ycord[index1], \
					  xtal->Zcord[index1]  };
	float atom2[3] = {xtal->Xcord[index2], \
					  xtal->Ycord[index2], \
					  xtal->Zcord[index2]  };
					  
	normalise_vector3(atom1);
	normalise_vector3(atom2);
	cross_vector3_vector3(avg, atom1, atom2);
	
	if( check_vec3_isNull(avg, TOL) )
		{return 0;}
	
	normalise_vector3(avg);
	return 1;
	
} 



int compute_average_axis(float axis1[3],float axis2[3], float avg[3])
{
	normalise_vector3(axis1);
	normalise_vector3(axis2);

	avg[0] = (axis1[0] + axis2[0]) / 2;
	avg[1] = (axis1[1] + axis2[1]) / 2;
	avg[2] = (axis1[2] + axis2[2]) / 2;
	
	if( check_vec3_isNull(avg, TOL) )
		{return 0;}
	
	normalise_vector3(avg);
	return 1;
	
}

void bring_molecules_to_origin(crystal* xtal)
{
	int N = xtal->num_atoms_in_molecule;
	int Z = xtal->Z;
	float com[3];
	
	for(int i = 0; i < N*Z; i = i + N)
	{
		compute_molecule_COM(*xtal, com, i);
		for(int j = 0; j < N; j++)
		{
			xtal->Xcord[i+j] -= com[0];
			xtal->Ycord[i+j] -= com[1];
			xtal->Zcord[i+j] -= com[2];
		}
		
	}
}

int compute_cross_axis(float axis1[3],float axis2[3], float avg[3])
{
	cross_vector3_vector3(avg , axis1, axis2);
	
	if( check_vec3_isNull(avg, TOL) )
		{return 0;}
	
	normalise_vector3(avg);
	return 1;
	
}



void find_overlap_list (crystal xtal,
						float first_com[3],
						float inverse_lattice_vectors[3][3],
						int overlap_list[],
						int *len_overlap_list,
						int Z_gen, 
						int N)
{
	//find molecules which overlap with the first molecule
	//The index of first atoms of those molecules are stored
	//in overlap list.
	
	int count = 0;
	
	for(int i = 0; i < Z_gen * N; i = i + N)
	{
		float com[3];
		compute_molecule_COM(xtal, com, i);

		if ( pdist_appx(xtal.lattice_vectors,
				inverse_lattice_vectors,
				com[0],
				com[1],
				com[2],
				first_com[0],
				first_com[1],
				first_com[2]) < CONST_TOL)
				
		{	
			overlap_list[count] = i;
			count++;
			//not tested; check if + or - !
			if (cart_dist(com, first_com) > CONST_TOL) 
			{
				float displace[3];
				vector3_subtract(com, first_com, displace);
				for(int k = i; k < i + N; k++)
				{
					xtal.Xcord[k] -= displace[0];
					xtal.Ycord[k] -= displace[1];
					xtal.Zcord[k] -= displace[2];  
				}
			}
		}
	}
	*len_overlap_list = count;
}



