from copy import deepcopy
import os
import numpy as np
from ase.data import vdw_radii,atomic_numbers
from ibslib.io import read
from Genarris.core.instruct import get_molecule_path

class MoleculeVolumeEstimator():
    """
    Class for performing molecular volume estimation using vdW radii
    
    Arguments
    ---------
    tol: float
        Tolerance to converge the estimation.
    iterations: int
        Maximum number of iterations to converge.
    batch: int
        Number of samples to process in MC method at a single time. 
    vdW: numpy.array
        Array of vdW volumes index by the atomic number.
        
    """
    def __init__(self, tol=1e-2, iterations=1e7, batch=1000, vdW=vdw_radii):
        # Change vdW radii
        self.iterations = int(iterations)
        self.batch = int(batch)
        self.vdW = vdW
        self.tol = tol
        
    
    def calc(self, molecule_struct):
        """
        Calculates the predicted of the molecular volume for a Structure object
        of a molecule. 
        
        """
        volume = self.MC(molecule_struct)
        # Adjust volume based on CSD analysis
        volume = self.linear_correction(volume)
        return volume
    
    
    def linear_correction(self, volume):
        """
        Applies linear correction to the input volume from CSD analysis.
        """
        slope = 1.3357455
        intercept = -0.00037431
        return volume*slope + intercept
        
    
    def MC(self, molecule_struct):
        """
        Monte Carlo method for volume estimation using vdW radii
        
        """
        self.struct = molecule_struct
        self.geo = molecule_struct.get_geo_array()
        
        # Get radii array for each element in geometry
        self.radii = np.array([self.vdW[atomic_numbers[ele]] for ele in 
                      molecule_struct.geometry["element"]])[:,None]
        
        # Get simulation sample region
        self._set_region()
        
        # Building iteration counter
        n_batchs = int(self.iterations / self.batch) 
        last_batch = self.iterations - n_batchs*self.batch
        batchs = [self.batch for x in range(n_batchs)]
        # Add a last batch if necessary to reach exactly self.iterations
        if last_batch != 0:
            batchs.append(last_batch)
        
        # Keep track of in an out
        self.vdW_in = 0
        self.vdW_out = 0
        self.volume_tracking = []
        for batch_size in batchs:
            points = self._sample_points(batch_size)
            # Pairwise distance between points and all atoms in geometry
            distances = points - self.geo[:,None]
            distances = np.linalg.norm(distances, axis=-1)
            
            # Find which samples were within a vdW radii from any atom
            vdW_dist = distances - self.radii
            vdW_in = vdW_dist <= 0
            vdW_in = np.sum(vdW_in, axis=0)
            vdW_in = np.where(vdW_in >= 1)[0]
            vdW_in = vdW_in.shape[0]
            self.vdW_in += vdW_in
            self.vdW_out += batch_size - vdW_in
            self.volume_tracking.append(
                    (self.vdW_in / self.vdW_out) * self.region_volume)
            
            # Check for convergence
            if len(self.volume_tracking) > 2:
                converged = self._check_convergence()
                if converged == False:
                    pass
                else:
                    break
        
        self.molecule_volume = (self.vdW_in / self.vdW_out) * self.region_volume
        return self.molecule_volume

    
    def _set_region(self):
        """ 
        Sets the sample region for the simulation
        
        """
        # Minimum and Max in each direction bias of 5
        min_region = np.min(self.geo,axis=0) - 5
        max_region = np.max(self.geo,axis=0) + 5
        self.region = np.array(list(zip(min_region,max_region)))
        
        self.region_volume = (self.region[0][1] - self.region[0][0])* \
                             (self.region[1][1] - self.region[1][0])* \
                             (self.region[2][1] - self.region[2][0])
        
    
    def _sample_points(self, batch_size):
        """
        Returns a sampling of the MC region for the given batch_size
        """
        x = np.random.uniform(self.region[0][0], self.region[0][1],
                              size=(batch_size,1))
        y = np.random.uniform(self.region[1][0], self.region[1][1],
                              size=(batch_size,1))
        z = np.random.uniform(self.region[2][0], self.region[2][1],
                              size=(batch_size,1))
        return np.concatenate((x,y,z),axis=1)
    
    
    def _check_convergence(self):
        """
        Checks is calculation has converged the volume of the system
        """
        # Difference of last two volume predictions
        err = abs(self.volume_tracking[-1] - self.volume_tracking[-2])
        if err > self.tol:
            return False
        else:
            return True


def estimate_unit_cell_volume(inst, comm):
    """
    Using the provided inst, if section pygenarris_structure_generation exists
    then create or overwrite option 'volume_mean' and 'volume_std'. If section
    pygenarris_structure_generation exists, then create or overwrite option
    'ucv_target', 'ucv_std', and 'ucv_ratio_range'.
    """

    sname = 'estimate_unit_cell_volume'
    verbose = inst.get_boolean(sname, 'verbose')
    Z = float(inst.get_eval(sname, 'Z'))
    if comm.rank == 0:
        mve = MoleculeVolumeEstimator()
        sname = 'estimate_unit_cell_volume'
        molecule_path = get_molecule_path(inst, sname)
        if verbose:
            print('molecule_path', molecule_path, flush=True)

        molecule = read(molecule_path)

        vol_estimate = mve.calc(molecule)
    else:
        vol_estimate = None
    vol_estimate = comm.bcast(vol_estimate, root=0)
    struct_vol_estimate = Z * vol_estimate
    std_of_predicted_errors = 0.062 * struct_vol_estimate
    std_to_use = 3.0 * std_of_predicted_errors
    vol_lower_bound = struct_vol_estimate - std_to_use
    vol_upper_bound = struct_vol_estimate + std_to_use
    ucv_ratio_range = [vol_lower_bound / struct_vol_estimate, vol_upper_bound / struct_vol_estimate]
    
    if inst.has_section('pygenarris_structure_generation'):
        if not inst.has_option('pygenarris_structure_generation', 'volume_mean'):
            inst.set('pygenarris_structure_generation', 'volume_mean', str(struct_vol_estimate))
        if not inst.has_option('pygenarris_structure_generation', 'volume_std'):
            inst.set('pygenarris_structure_generation', 'volume_std', str(std_to_use))
        if comm.rank == 0:
            print('volume_mean', inst.get('pygenarris_structure_generation', 'volume_mean'), flush=True)
            print('volume_std', inst.get('pygenarris_structure_generation', 'volume_std'), flush=True)
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
