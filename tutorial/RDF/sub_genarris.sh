#!/bin/bash                                                                                      
#                                                                                                
#SBATCH -J genarris_test # Job name
#SBATCH -n 56  # Number of total cores

#SBATCH -N 1 # Number of nodes                                          

# SBATCH --mem=0
# SBATCH --mem-per-cpu=2285
#SBATCH --mem-per-cpu=1100 # Memory pool for all cores in MB (see also --mem-per-cpu)#2285
#SBATCH -o j_%j.out # File to which STDOUT will be written %j is the job #                                       
# SBATCH -p cpu
#SBATCH -p idle
# SBATCH -p debug 

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
export LD_PRELOAD=${MKLROOT}/lib/intel64_lin/libmkl_intel_lp64.so:${MKLROOT}/lib/intel64_lin/libmkl_sequential.so:${MKLROOT}/lib/intel64_lin/libmkl_core.so:${MKLROOT}/lib/intel64_lin/libmkl_blacs_intelmpi_lp64.so:${MKLROOT}/lib/intel64_lin/libmkl_scalapack_lp64.so:${I_MPI_ROOT}/intel64/lib/libmpi.so.12

mpirun -n 5 python -u ${HOME}/genarris_mpi4py/merged_4_17_2019/Genarris/Genarris/genarris_master.py ui.conf
echo " "
echo "Job Ended at `date`"
