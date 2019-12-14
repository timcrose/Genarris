
%module pygenarris
%{
#include "read_input.h"
#include "spg_generation.h"
#include "pygenarris.h"
#include "combinatorics.h"
#include "check_structure.h"
#include "crystal_utils.h"
#include "randomgen.h"
#include "spglib.h"
#include <omp.h>

%}


%{
#define SWIG_FILE_WITH_INIT
%}
%include "numpy.i"
%init %{
import_array();
%}

void generate_molecular_crystals(char *filename,int num_structures, int Z,
	double volume_mean1, double volume_std1, double sr1, double tol1, 
	int max_attempts);
	
void find_allowed_positions_using_molecular_symmetry(char mol_sym[6],
	int Z, int Zpp);
	
void allocate_xtal(crystal* xtal, int Z, int N);

%include "spg_generation.h"


%apply (double INPLACE_ARRAY2[ANY][ANY]) {(double lattice_vector[3][3])};
%apply (double* IN_ARRAY1, int DIM1) {(double *Xc, int total_atoms1)};
%apply (double* IN_ARRAY1, int DIM1) {(double *Yc, int total_atoms2)};
%apply (double* IN_ARRAY1, int DIM1) {(double *Zc, int total_atoms3)};
void create_crystal_from_array(crystal *xtal, double lattice_vector[3][3], double *Xc,int total_atoms1,
		double *Yc,int total_atoms2, double *Zc, int total_atoms3,char *atoms,  int total_atoms, int Z, int spg);


%apply (double INPLACE_ARRAY1[ANY]) {(double a[3])};

void print_crystal(crystal* xtal);
void free_xtal(crystal* xtal);
int c_check_structure(crystal xtal, double sr);


%apply ( float* IN_ARRAY2, int DIM1, int DIM2) {(float *vdw_matrix, int dim1, int dim2)};
int check_structure_with_vdw_matrix(crystal random_crystal,
	float *vdw_matrix,
	int dim1,
	int dim2);

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

