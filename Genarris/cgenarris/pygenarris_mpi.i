%module pygenarris_mpi
%{
#include "read_input.h"
#include "spg_generation.h"
#include "pygenarris.h"
#include "pygenarris_mpi.h"
#include "combinatorics.h"
#include "check_structure.h"
#include "crystal_utils.h"
#include "randomgen.h"
#include "spglib.h"
#include "mpi.h"
#include "cgenarris_mpi.h"

%}

%{
#define SWIG_FILE_WITH_INIT
%}
%include "numpy.i"
%init %{
import_array();
%}

%include "mpi4py/mpi4py.i"
%mpi4py_typemap(Comm, MPI_Comm);

%apply ( float* IN_ARRAY2, int DIM1, int DIM2) {(float *vdw_matrix, int dim1, int dim2)};

void mpi_generate_molecular_crystals_with_vdw_cutoff_matrix(
	float *vdw_matrix,
	int dim1,
	int dim2,
	int num_structures,
	int Z,
	double volume_mean1,
	double volume_std1,
	double tol1, 
	long max_attempts, 
	MPI_Comm world_comm);

void send_xtal(MPI_Comm comm, int destination, crystal* xtal, int total_atoms);
void receive_xtal(MPI_Comm comm, int source, crystal* xtal, int total_atoms);
