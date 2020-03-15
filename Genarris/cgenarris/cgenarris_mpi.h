
#ifndef CGENARRIS_MPI_H
#define CGENARRIS_MPI_H


void send_xtal(MPI_Comm comm, int destination, crystal* xtal, int total_atoms);
void receive_xtal(MPI_Comm comm, int source, crystal* xtal, int total_atoms);

#endif
