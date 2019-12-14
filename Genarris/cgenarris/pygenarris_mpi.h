#ifndef PYGENARRIS_MPI_H
#define PYGENARRIS_MPI_H

#include "mpi.h"

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


#endif
