#!/bin/bash
#
#SBATCH -J job # Job name
#SBATCH -n 112 # Number of total cores

#SBATCH -N 2 # Number of nodes

#SBATCH --mem-per-cpu=2142
#SBATCH -o j_%j.out # File to which STDOUT will be written %j is the job #
#SBATCH -p cpu

mpirun -np 112 ./aims.160328_3_intel_mkl.scalapack.mpi.x > aims.out
