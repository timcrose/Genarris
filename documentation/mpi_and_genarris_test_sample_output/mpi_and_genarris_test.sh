#!/bin/bash
#SBATCH -J genarris_test # Job name
#SBATCH -n 2  # Number of total cores
#SBATCH -N 1 # Number of nodes                                          
#SBATCH --mem-per-cpu=1100
#SBATCH -o j_%j.out # File to which STDOUT will be written %j is the job #                                       
#SBATCH -p debug 

mpirun -n 2 python mpi_and_genarris_test.py
