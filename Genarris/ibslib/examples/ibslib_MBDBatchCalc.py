

import json
from ase.calculators.aims import Aims
from ibslib.calculators import Slurm,MBDBatchCalc, \
                               tin_arguments,mbd_settings

# Directory with structure files to calculate
struct_dir = "struct_dir"
# Directory where calculations will take place
calc_dir = "batch_calc"
# Setup command for calculations
tin_arguments["command"] = "mpirun -np 12 /home/ibier/Software/MBD/DFT_MBD_AT_rsSCS.x geometry.xyz setting.in > mbd.out" 
slurm = Slurm(tin_arguments)

calc = MBDBatchCalc(struct_dir, settings=mbd_settings, Slurm=slurm, calc_dir=calc_dir)
calc.calc(overwrite=True)
