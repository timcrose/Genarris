#!/bin/bash                                                                                      
#                                                                                                
#SBATCH -J genarris_test # Job name
#SBATCH -n 56  # Number of total cores
#SBATCH -N 1 # Number of nodes                                          
#SBATCH --mem-per-cpu=1000
#SBATCH -o j_%j.out # File to which STDOUT will be written %j is the job #                                       
#SBATCH -p idle

#source deactivate
module load impi/2018_Update_3 intel/18.0.3.222

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

mpirun -n 56 python -u ../../Genarris/genarris_master.py ui.conf
echo " "
echo "Job Ended at `date`"
