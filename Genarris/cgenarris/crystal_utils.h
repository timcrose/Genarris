#ifndef CRYSTAL_UTILS_H
#define CRYSTAL_UTILS_H

void print_crystal(crystal* xtal);
void print_crystal2file(crystal* xtal, FILE* out_file);
void compute_molecule_COM(crystal xtal, float com[3], int i);
void remove_close_molecules(crystal* xtal);
void remove_close_molecules(crystal* xtal);
void remove_close_molecules_mirror_plane(crystal *xtal );
void convert_xtal_to_fractional(crystal *xtal);
void bring_all_molecules_to_first_cell(crystal *xtal);
void convert_xtal_to_cartesian(crystal *xtal);
void combine_close_molecules(crystal* xtal);
void average_positions_using_list(float 	*tempx,
								  float 	*tempy,
								  float 	*tempz,
								  char	 	*tempatom,
								  int   	 molecule_counter,
								  int		 N,
								  int	 	 Z,
								  crystal 	*xtal,
								  int 		*averaging_list);
int detect_spg_using_spglib(crystal* xtal);
void allocate_xtal(crystal* xtal, int Z, int N);
void free_xtal(crystal* xtal);
void copy_xtal(crystal* xtal1, crystal* xtal2);


#endif
