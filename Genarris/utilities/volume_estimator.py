from copy import deepcopy
import os
import numpy as np
from ase.data import vdw_radii,atomic_numbers
from ibslib.io import read
from ibslib.molecules.volume_estimator import MoleculeVolumeEstimator
from Genarris.core.instruct import get_molecule_path

def estimate_unit_cell_volume(inst, comm):
    """
    Using the provided inst, if section pygenarris_structure_generation exists
    then create or overwrite option 'volume_mean' and 'volume_std'. If section
    pygenarris_structure_generation exists, then create or overwrite option
    'ucv_target', 'ucv_std', and 'ucv_ratio_range'.
    
    
    
    """
    print("---------------------------------------------------------------------------------------------------")
    print("Begin Unit Cell Volume Estimation",flush=True)
    print("---------------------------------------------------------------------------------------------------")
    sname = 'estimate_unit_cell_volume'
    verbose = inst.get_boolean(sname, 'verbose')
    Z = float(inst.get_eval(sname, 'Z'))
    if comm.rank == 0:
        mve = MoleculeVolumeEstimator()
        sname = 'estimate_unit_cell_volume'
        molecule_path = get_molecule_path(inst, sname)
        if verbose:
            print('molecule path :', molecule_path, flush=True)

        molecule = read(molecule_path)

        vol_estimate = mve.calc(molecule)
    else:
        vol_estimate = None
    vol_estimate = comm.bcast(vol_estimate, root=0)
    struct_vol_estimate = Z * vol_estimate
    # Make an option but default is from underlying distribution

    std_of_predicted_errors = inst.get_with_default(sname, 
                                'std_of_predicted_errors', 0.025, eval=True)
    std_to_use = 3.0 * std_of_predicted_errors * struct_vol_estimate
    vol_lower_bound = struct_vol_estimate - std_to_use
    vol_upper_bound = struct_vol_estimate + std_to_use
    ucv_ratio_range = [vol_lower_bound / struct_vol_estimate, vol_upper_bound / struct_vol_estimate]
    
    if inst.has_section('pygenarris_structure_generation'):
        if not inst.has_option('pygenarris_structure_generation', 'volume_mean'):
            inst.set('pygenarris_structure_generation', 'volume_mean', str(struct_vol_estimate))
        if not inst.has_option('pygenarris_structure_generation', 'volume_std'):
            inst.set('pygenarris_structure_generation', 'volume_std', str(std_to_use))
        if comm.rank == 0:
            print('volume mean:', inst.get('pygenarris_structure_generation', 'volume_mean'), flush=True)
            print('volume standard deviation:', inst.get('pygenarris_structure_generation', 'volume_std'), flush=True)
    if inst.has_section('structure_generation_batch'):
        if not inst.has_option('structure_generation_batch', 'ucv_target'):
            inst.set('structure_generation_batch', 'ucv_target', str(struct_vol_estimate))
        if not inst.has_option('structure_generation_batch', 'ucv_std'):
            inst.set('structure_generation_batch', 'ucv_std', str(std_to_use))
        if not inst.has_option('structure_generation_batch', 'ucv_ratio_range'):
            inst.set('structure_generation_batch', 'ucv_ratio_range', str(ucv_ratio_range))
        if comm.rank == 0:
            print('ucv_target', inst.get('structure_generation_batch', 'ucv_target'), flush=True)
            print('ucv_std', inst.get('structure_generation_batch', 'ucv_std'), flush=True)
            print('ucv_ratio_range', inst.get('structure_generation_batch', 'ucv_ratio_range'), flush=True)
    
    return inst

    
if __name__ == "__main__":
    from ibslib.io import read,write
    import matplotlib.pyplot as plt
    struct_dict = read("Single_Molecules/json")
    
    x = []
    y = []
    mve = MoleculeVolumeEstimator()
    for struct_id,struct in struct_dict.items():
        exp_volume = struct.get_property("molecule_volume")
        pred_volume = mve.calc(struct)
        print("{} {} {}".format(exp_volume, pred_volume, exp_volume-pred_volume))
        x.append(exp_volume)
        y.append(pred_volume)
    
    plt.scatter(x,y)
