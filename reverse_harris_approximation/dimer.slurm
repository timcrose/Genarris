#!/bin/bash
#SBATCH --qos=normal
#SBATCH --job-name=noname
#SBATCH --time=0-06:00:00
#SBATCH -o aims.log
#SBATCH -e aims.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20

module unload tulane/intel-psxe/2015
module load intel-psxe/2015-update1
module load anaconda
module load gcc/4.8.2

############ THE JOB ITSELF #############################
echo `date`
echo Start Job
echo $SLURM_NODELIST
export OMP_NUM_THREADS=1

python /home/xli20/isgep_dev/isgep/src/isgep_master.py /home/xli20/isgep_dev/isgep/examples/reverse_harris_approximation/dimer.conf 

echo End Job
echo `date`
