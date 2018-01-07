#!/bin/bash
#SBATCH --qos=normal
#SBATCH --job-name=noname
#SBATCH --time=0-06:00:00
#SBATCH -o run.log
#SBATCH -e run.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20


############ THE JOB ITSELF #############################
echo `date`
echo Start Job
echo $SLURM_NODELIST
export OMP_NUM_THREADS=1

python /path/to/genarris/src/genarris_master.py ui.conf 

echo End Job
echo `date`
