#ifndef MOLECULE_PLACEMENT_H
#define MOLECULE_PLACEMENT_H



void apply_all_symmetry_ops(crystal *xtal,
							molecule *mol,
							float* mol_Xfrac,
							float *mol_Yfrac,
							float *mol_Zfrac,
							int N,
							int hall_number);

int auto_align_and_generate_at_position(crystal *Xtal,
								molecule *mol,
								int hall_number,
								int spg,
								int wyckoff_pos,
								COMPATIBLE_SPG compatible_spg);
								
int align_using_std_orientations(crystal* xtal_1,
								molecule* mol,
								int hall_number,
								float first_com[3],
								int overlap_list[],
								int len_overlap_list,
								int dof,
								COMPATIBLE_AXES compatible_axes);
							  
int get_degrees_of_freedom(int spg, int pos);

#endif

