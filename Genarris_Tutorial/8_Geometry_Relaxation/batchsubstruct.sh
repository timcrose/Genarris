#!/bin/bash

#
# - This script will sub a aims job for each of the geometry files found in structdir
# - Makes sure to change the dirname name slicing to fit your current naming scheme
# - structdir =  the directory where all the geometry files you want to run are. Usually, 
#                 this directory will contain a copy of all the geomtry files you want to run
# - controlpath = path to the control file you want to use
# - slurmpath = path to the slurm submission script you want to use
# - slurmname = name of the slurm submission script 

structdir="/home/ibier/genarris-runs/FUQJIK/8mpc/larger_sr_generation/8_Geometry_Relaxation/batch_eval/"
controlpath="/home/ibier/genarris-runs/FUQJIK/8mpc/larger_sr_generation/8_Geometry_Relaxation/control.in.FULL.k333"
slurmpath="/home/ibier/genarris-runs/FUQJIK/8mpc/larger_sr_generation/8_Geometry_Relaxation/sub_aims_2018.sh"
slurmname="sub_aims_2018.sh"

cd "$structdir"

for struct in *
do

dirname=`echo init_$struct | sed 's/FUQJIK_8mpc_//' | sed 's/_geometry//'`
name=`echo $dirname | sed 's/init_//'`
mkdir "$dirname"

mv "$struct" "$dirname"/geometry.in

cd "$dirname"
cp "$controlpath" ./control.in
cp "$slurmpath" .
sed -i "2i\#SBATCH -J aims$name # Job name" "$slurmname";
sbatch "$slurmname"

cd "$structdir"
echo "$dirname"

done
