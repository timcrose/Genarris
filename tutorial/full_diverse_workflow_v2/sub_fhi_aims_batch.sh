#!/bin/bash                                                                                      
#                                                                                                
#SBATCH -J genarris_test # Job name
#SBATCH -n 56  # Number of total cores

#SBATCH -N 1 # Number of nodes                                          

#SBATCH --mem=0
# SBATCH --mem-per-cpu=2285
# SBATCH --mem-per-cpu=1070 # Memory pool for all cores in MB (see also --mem-per-cpu)#2285
#SBATCH -o j_%j.out # File to which STDOUT will be written %j is the job #                                       
#SBATCH -p cpu
# SBATCH -p debug 

#source deactivate
ulimit -s unlimited
ulimit -v unlimited
#ulimit -l unlimited
export PYTHONUNBUFFERED=TRUE
#./clean.sh
#rm -r SPE_structure
#./clean_aimspy.sh
#mpirun -n 14 python aimspy.py aimsconf.conf
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
#mpirun -n 14 python parallelaims_standalone.py > aims.out
#/opt/ohpc/pub/intel/intel16u4/compilers_and_libraries_2016.4.258/linux/mpi/intel64/bin/mpirun -np 32 python parallelaims.py > aims.out 
#mpirun -n 56 python -u  ${HOME}/test/mpi4py_test/test_comm.py
#mpirun -n 56 python -u /home/trose/genarris20_mpi4py/genarris_mpi4py/src/evaluation/run_aims.py fhi_aims_batch_run.conf
module list
mpirun -n 57 ./preload_scripts.sh
echo " "
echo "Job Ended at `date`"
