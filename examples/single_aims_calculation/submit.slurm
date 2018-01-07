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

############ THE JOB ITSELF #############################
echo `date`
echo Start Job
echo $SLURM_NODELIST
export OMP_NUM_THREADS=1

python /home/xli20/isgep_dev/isgep/src/isgep_master.py /home/xli20/isgep_dev/isgep/examples/single_aims_calculation/ui.conf 

echo End Job
echo `date`
