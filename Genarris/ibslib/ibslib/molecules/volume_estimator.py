

import numpy as np
from ase.data import vdw_radii,atomic_numbers
from ase.data.colors import jmol_colors

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from ibslib.descriptor import BondNeighborhood


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
        This is also the number for which convergence is tested. 
    vdW: numpy.array
        Array of vdW volumes index by the atomic number.
        
    """
    def __init__(self, tol=1e-2, iterations=1e8, batch=100000, vdW=vdw_radii):
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
    
    
    def linear_correction(self,volume):
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
        self.total_samples = 0
        self.volume_tracking = []
        self.var_tracking = []
        self.points_inside = []
        self.c = []
        for batch_size in batchs:
            points = self._sample_points(batch_size)
            # Pairwise distance between points and all atoms in geometry
            distances = points - self.geo[:,None]
            distances = np.linalg.norm(distances, axis=-1)
            
            # Find which samples were within a vdW radii from any atom
            vdW_dist = distances - self.radii
            
            self.vdW_dist = vdW_dist
            r_idx,in_idx = np.where(vdW_dist <= 0)
            self.points_inside.append(points[in_idx])
            self.c.append(self._define_colors(self.struct.geometry["element"][r_idx]))
            
            vdW_in = vdW_dist <= 0
            vdW_in = np.sum(vdW_in, axis=0)
            vdW_in = np.where(vdW_in >= 1)[0]
            vdW_in = vdW_in.shape[0]
            self.vdW_in += vdW_in
            self.total_samples += batch_size
            self.volume_tracking.append(
                    (self.vdW_in / self.total_samples) * self.region_volume)
            self.var_tracking.append(self.volume_tracking[-1] * \
                            (self.region_volume - self.volume_tracking[-1]) \
                            / self.total_samples)
            
            # Check for convergence
            if len(self.volume_tracking) > 2:
                converged = self._check_convergence()
                if converged == False:
                    pass
                else:
                    break
        
        self.molecule_volume = (self.vdW_in / self.total_samples) * self.region_volume
        
        # Saving points and colors for plotting 
        self.points_inside = np.vstack(self.points_inside).reshape(-1,3)
        self.c = np.vstack(self.c)
        
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
        

    def plot_MC(self, elev=0, azim=0, figname=""):
        """
        Plot the points which were found to be inside the molecule. Useful
        for visualizing the volume which is found by MC method.
        
        """
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111, projection='3d')
        ax.w_xaxis.set_pane_color((0., 0., 0., 1.0))
        ax.w_yaxis.set_pane_color((0., 0., 0., 1.0))
        ax.w_zaxis.set_pane_color((0., 0., 0., 1.0))
        ax.scatter(self.points_inside[:,0],
                   self.points_inside[:,1],
                   self.points_inside[:,2], alpha=0.1,
                   c = self.c)
        ax.grid(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
#        for ii in xrange(0,360,10):
        ax.view_init(elev=elev, azim=azim) 
        plt.tight_layout()
        if len(figname) > 0:
            fig.savefig(figname, bbox_inches='tight',
                        transparent=True, pad_inches=0)
        
    
    def _define_colors(self,elements):
        idx = [atomic_numbers[x] for x in elements]
        return jmol_colors[idx]
    
    
    
if __name__ == "__main__":
    pass
#    from ibslib.io import read,write
#    import matplotlib.pyplot as plt
#    struct_dict = read("Single_Molecules/json")
#    
#    x = []
#    y = []
#    mve = MoleculeVolumeEstimator()
#    for struct_id,struct in struct_dict.items():
#        exp_volume = struct.get_property("molecule_volume")
#        pred_volume = mve.calc(struct)
#        print("{} {} {}".format(exp_volume, pred_volume, exp_volume-pred_volume))
#        x.append(exp_volume)
#        y.append(pred_volume)
#    
#    plt.scatter(x,y)