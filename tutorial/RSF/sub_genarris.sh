#!/bin/bash                                                                                      
#                                                                                                
#SBATCH -J genarris_test # Job name
#SBATCH -n 24  # Number of total cores
#SBATCH -N 1 # Number of nodes                                          
#SBATCH --mem=0
#SBATCH -o j_%j.out # File to which STDOUT will be written %j is the job #                                       
#SBATCH -p Manny_is_supercool

#source deactivate
ulimit -s unlimited
ulimit -v unlimited
export PYTHONUNBUFFERED=TRUE
echo "PATH"
echo $PATH
echo " "
echo "PYTHONPATH"
echo $PYTHONPATH
echo " "
echo "LD_LIBRARY_PATH"
echo $LD_LIBRARY_PATH
echo `which python`
echo `which mpirun`
export OMP_NUM_THREADS=1
module list

mpirun -n 24 python -u ../../Genarris/genarris_master.py ui.conf
echo " "
echo "Job Ended at `date`"
