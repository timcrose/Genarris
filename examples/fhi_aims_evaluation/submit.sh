#!/bin/bash
#COBALT -n 10
#COBALT -t 60
#COBALT -q default
#COBALT -A HybridPV_tesp

#proc=64
#bin=/home/fcurtis/fhi-aims.160328_3/bin-fhik/aims.160328_3.scalapack.mpi.x

echo "Running Cobalt Job $COBALT_JOBID."

python ../../src/genarris_master.py ./theta.conf 
#aprun -n $(($COBALT_PARTSIZE*$proc)) -N $proc -e OMP_NUM_THREADS=1 $bin > aims.out

#EXIT_STATUS=$?
echo "Job $COBALT_JOBID completed."
exit $EXIT_STATUS
