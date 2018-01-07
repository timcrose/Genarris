#!/bin/bash
#SBATCH --qos=normal
#SBATCH --job-name=noname
#SBATCH --time=0-06:00:00
#SBATCH -o aims.log
#SBATCH -e aims.err
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=20

#module unload tulane/intel-psxe/2015
#module load intel-psxe/2015-update1

module load anaconda
module load gcc/4.8.2

############ THE JOB ITSELF #############################
echo `date`
echo Start Job
echo $SLURM_NODELIST
export OMP_NUM_THREADS=1
export LD_PRELOAD=/share/apps/gcc/4.7.4/lib64/libgfortran.so.3

python ../../src/genarris_master.py ./ui.conf 

echo End Job
echo `date`
