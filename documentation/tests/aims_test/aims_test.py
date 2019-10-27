from mpi4py import MPI
import numpy
import aims_w

def main():
    comm = MPI.COMM_WORLD
    commf = comm.py2f()
    aims_w.aims_w(commf)

main()
