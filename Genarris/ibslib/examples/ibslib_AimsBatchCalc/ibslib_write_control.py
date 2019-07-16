

import json

from ibslib.io import read,write
from ibslib.calculators import tier1_SPE_settings,tier1_relaxed_settings

from ase.calculators.aims import Aims

struct = read("C60_Eval.json")
atoms = struct.get_ase_atoms()

calc = Aims(**tier1_SPE_settings)
calc.write_input(atoms)