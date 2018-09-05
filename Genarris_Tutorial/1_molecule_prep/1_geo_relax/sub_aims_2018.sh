#!/bin/bash       
#SBATCH -J DUX # Job name
#SBATCH -n 56 # Number of total cores
#SBATCH -N 1 # Number of nodes                                          
#SBATCH --mem=0 # Memory pool for all cores in MB (see also --mem-per-cpu)                        
#SBATCH -o j_%j.out # File to which STDOUT will be written %j is the job                                       
#SBATCH -p cpu


module load genarris
mpirun -n 56 ${aims_bin} > aims.out
