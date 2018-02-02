#!/bin/bash
#COBALT -t 05:30:00
#COBALT -n 8
#COBALT -A HybridPV_tesp
#COBALT -q cache-quad

proc=64
thr=1
bin=/home/trose/bin-o3intel-theta/aims.170418.theta
module add atp
export  OMP_NUM_THREADS=$thr
ntot=$(($COBALT_PARTSIZE*$proc))
aprun  -n $ntot -N $proc -cc depth -d $thr    $bin  >& aims.out

