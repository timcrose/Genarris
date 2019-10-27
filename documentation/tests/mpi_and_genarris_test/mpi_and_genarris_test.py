from mpi4py import MPI
import Genarris

comm = MPI.COMM_WORLD
print("comm.rank", comm.rank)
