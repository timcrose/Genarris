

import os
import copy

# Base is to use the Aims class for calculations
from ase.calculators.aims import Aims

from ibslib.io import read,write,check_dir


tier1_SPE_settings = \
    {
        # Settings for ASE 
        "restart": None, 
        "ignore_bad_restart_file": False,
        "atoms": None, 
        "cubes": None, 
        "radmul": None,
        # aims_command not specified because it will be handled by Slurm
        "aims_command": "None",

        # Folder settings 
        "directory": "./",
        "label": "", 
        "outfilename": "aims.out", 
        "species_dir": "/home/ibier/aims/fhi-aims.171221_1/species_defaults/light/",

        # Calculation settings
        "tier": None, 
        "xc": "pbe",
        "spin": "none",
        "relativistic": "atomic_zora scalar",
        "charge": 0,

        "occupation_type": "gaussian 0.01",
        "mixer": "pulay",
        "n_max_pulay": 8,
        "charge_mix_param": 0.02,
        "sc_accuracy_rho": 1e-4,
        "sc_accuracy_eev": 1e-2,
        "sc_accuracy_etot": 1e-6,
        "sc_iter_limit": 10000,

        "KS_method": "parallel",
        "empty_states": 6,
        "basis_threshold": 1e-5,
        "k_grid": [3,3,3],
        "vdw_correction_hirshfeld": True,
    }

tier1_relaxed_settings = copy.deepcopy(tier1_SPE_settings)
tier1_relaxed_settings["sc_accuracy_forces"] = 1e-4
tier1_relaxed_settings["relax_geometry"] = "trm 1e-2"
tier1_relaxed_settings["relax_unit_cell"] = "full"
tier1_relaxed_settings["hessian_to_restart_geometry"] = False
tier1_relaxed_settings["harmonic_length_scale"] = 0.01
tier1_relaxed_settings["energy_tolerance"] = 5e-4





class AimsBatchCalc():
    """
    Class for starting calculations of an entire directory using Slurm.
    """
    def __init__(self, struct_dir, Aims=None, Slurm=None, calc_dir=""):
        """
        Arguments
        ---------
        struct_dir: str 
            Path to directory of Structures to be calculated
        Aims: ase.calculators.aims.Aims
            Aims class initialized with all the parameters desired for the 
            calculation. 
        Slurm: ibslib.calculators.slurm.Slurm
            Slurm class initialized with all parameters
        calc_dir: str
            Path to directory where calculations should take place. 
            Default is to calculate in batch_calc folder. calc_dir
            is not created until BatchCalc.calc or BatchCalc.calc_struct
            is called. 

        """
        if Aims == None:
            raise Exception("Must supply initialized ase.calculators.aims.Aims"+ 
                        " class to ibslib.calculators.aims.AimsBatchCalc.")
        if Slurm == None:
            raise Exception("Must supply initialized ibslib.calculators.slurm.Slurm"+ 
                        " class to ibslib.calculators.aims.AimsBatchCalc.")
        self.Aims = Aims
        self.Slurm = Slurm
        self.struct_dir = struct_dir

        if len(calc_dir) == 0:
            self.calc_dir = "batch_calc"
        else:
            self.calc_dir = calc_dir

        self.struct_dict = read(struct_dir)


    def calc(self, calc_dir="",overwrite=False):
        """
        Arguments
        ---------
        calc_dir: str
            Optional to augment the calc dir in calc method. 
        overwrite: bool
            Only for overwrite control of the Slurm submission script. 
            Overwriting control.in or geometry.in is handled by ASE.
        """
        if len(calc_dir) != 0:
            self.calc_dir = calc_dir

        check_dir(self.calc_dir)
        
        for struct_id,struct in self.struct_dict.items():
            self.calc_struct(struct_id, overwrite=overwrite)
        
    
    def calc_struct(self, struct, calc_dir="",overwrite=False):
        """
        Calculate a single structure in the struct_dict

        Arguments
        ---------
        struct: Structure
            Structure to calculate
        """
        if len(calc_dir) != 0:
            self.calc_dir = calc_dir
        
        struct_id = struct.struct_id

        # Make calculation directory
        check_dir(self.calc_dir)
        # Make structure directory
        struct_path = os.path.join(self.calc_dir, struct_id)
        check_dir(struct_path)

        atoms = struct.get_ase_atoms()

        # Move to calculation directory
        cwd_temp = os.path.abspath(os.getcwd())
        os.chdir(struct_path)

        self.Aims.directory = os.curdir
        self.Aims.write_input(atoms)
        self.Slurm.calc("Submit.sh", overwrite=overwrite)

        os.chdir(cwd_temp)



    
        

        
