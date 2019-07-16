

import json
from ase.calculators.aims import Aims
from ibslib.calculators import Slurm,AimsBatchCalc, \
                               tin_arguments,tier1_SPE_settings

# Directory with structure files to calculate
struct_dir = "struct_dir"
# Directory where calculations will take place
calc_dir = "batch_calc"
# Initialize settings for Slurm and Aims
aims = Aims(**tier1_SPE_settings)
# Setup command for calculations
#tin_arguments["pre-command"] = "module load intel"
tin_arguments["command"] = "mpirun -np 12 /opt/marom/bin/aims.171221_1.scalapack.mpi.x > aims.out"
slurm = Slurm(tin_arguments)

calc = AimsBatchCalc(struct_dir, Aims=aims, Slurm=slurm, calc_dir=calc_dir)
calc.calc(overwrite=True)
