#ifndef PYGENARRIS_H
#define PYGENARRIS_H



void generate_molecular_crystals(char *filename,int num_structures, int Z,
	double volume_mean1, double volume_std1, double sr1, double tol1, 
	int max_attempts);

int c_check_structure(crystal xtal, double sr);

void create_crystal_from_array(crystal *xtal, double lattice_vector[3][3], double *Xc,int total_atoms1,
		double *Yc,int total_atoms2, double *Zc, int total_atoms3,char* atoms, int total_atoms, int Z, int spg);

//Generator with van der waal cutoff distance matrix
void generate_molecular_crystals_with_vdw_cutoff_matrix(char *filename,
	int seedstate,
	float *vdw_matrix,
	int dim1,
	int dim2,
	int num_structures,
	int Z,
	double volume_mean1,
	double volume_std1,
	double tol1, 
	int max_attempts);

int num_compatible_spacegroups(int Z, double tolerance);

#endif
