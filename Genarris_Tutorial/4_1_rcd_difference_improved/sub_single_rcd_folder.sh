#!/bin/bash       
#SBATCH -J test_job # Job name
#SBATCH -n 1 # Number of total cores
#SBATCH -N 1 # Number of nodes                                          
#SBATCH --mem=2142 # Memory pool for all cores in MB (see also --mem-per-cpu)                        
#SBATCH -o j_%j.out # File to which STDOUT will be written %j is the job                                       
#SBATCH -p cpu

module load genarris
python ${genarris_master} rcd_folder_inner.conf
