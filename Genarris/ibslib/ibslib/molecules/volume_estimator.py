

import numpy as np
from ase.data import vdw_radii,atomic_numbers


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
    
    
    def linear_correction(volume):
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


"""
From: https://blog.levilentz.com/integrated-charge-density-from-cube-file/
"""

class Cube_Class():
    """
    Class for running molecule volume estimation using a CUBE file
    """
    def __init__(self, cube_path, charge_tol=1e-3):
        """
        charge_tol: float
            Value which determines separation between molecule and vacuum
        """
        self.charge_tol = charge_tol

        ######  These are all set in _read   ######
        # N = number of incremenets in each direction
        self.N = np.zeros(3)
        # N_matrix = 3x3 of increment vectors stored in rows
        self.N_matrix = np.zeros((3,3))
        # Atom coordinates from cube file
        self.atom_coords = np.array([])
        # Atom elements as atomic number
        self.atom_ele = []
        # Atom charges
        self.atom_charge = []
        # 3D grid will be constructed using these values in _read
        self.charge_grid = np.array([])
        # Total volume of cube file
        self.total_volume = 1

        # Read in cube file info
        self._read(cube_path)

    
    def _read(self, cube_path):
        """
        From: https://blog.levilentz.com/integrated-charge-density-from-cube-file/
        """
        with open(cube_path,'r') as f:
            # Skip first two lines
            next(f)
            next(f)

            # NAtoms, X-Origin, Y-Origin, Z-Origin
            line = next(f)
            split_line = line.split()
            self.n_atoms = int(split_line[0])

            # Next three lines have format: 
            # N1, X1, Y1, Z1   increments in first direction
            for i in range(0,3):
                line = next(f).split()
                self.N[i] = int(line[0])
                for j in range(1,4):
                    self.N_matrix[i][j-1] = float(line[j])
            
            # Initialize atom_coords and read in positions
            self.atom_coords = np.zeros((self.n_atoms,3))
            for i in range(0,self.n_atoms):
                line = next(f).split()
                self.atom_ele.append(int(line[0]))
                self.atom_charge.append(float(line[1]))
                for j in range(0,3):
                    self.atom_coords[i][j] = float(line[j+2])
            
            # Reading grid in
            lines = f.readlines()
            lines = [x.split() for x in lines]
            cube_grid = []
            for line in lines:
                for value in line:
                    cube_grid.append(float(value))
            cube_grid = np.array(cube_grid)
            cube_grid = cube_grid.ravel()
            # Reshaping based on N 
            columns = np.arange(0,self.N[2])
            rows = np.arange(0,self.N[1])*self.N[2]
            depth = np.arange(0,self.N[0])*(self.N[1]*self.N[2])
            broadcast = depth[:,None,None] + rows[:,None] + columns
            self.cube_grid = np.take(cube_grid, broadcast.astype(int))
        
        for i in range(0,3):
            self.total_volume *= self.N[i]*np.linalg.norm(self.N_matrix[i,:])
    
    
    def run_MC(self, n_iter=1e7, n_batch=100000):
        """
        MC iterations can be performed in batches to speed up calculation. 

        Arguments
        ---------
        n_inter: int
            Number of iterations to perform
        n_batch: int
            Number of batches to perform at a single time. Limiting factor 
            will be the memory of the machine being used. Allocated matrix
            of samples has to fit in memory. 

        """
        # Iteration counter
        batch_iter = int(n_iter / n_batch) + 1
        self.n_in = 0
        self.n_out = 0
        for i in range(batch_iter):
            samples = self._sample_point(n_batch)
            results = self.cube_grid[samples[:,0],
                                     samples[:,1],
                                     samples[:,2]]

            # updating counts based on results
            n_in = np.sum(results > self.charge_tol)
            n_out = n_batch - n_in
            self.n_in += n_in
            self.n_out += n_out
        
        self.molecule_volume = (self.n_in / (self.n_in + self.n_out)) * self.total_volume
        return self.molecule_volume

    
    def _sample_point(self, n_batch=1000):
        """ Perform batch sampling
        """
        x = np.random.random_integers(0,int(self.N[0])-1, size=(n_batch,1))
        y = np.random.random_integers(0,int(self.N[1])-1, size=(n_batch,1))
        z = np.random.random_integers(0,int(self.N[2])-1, size=(n_batch,1))
        return np.concatenate((x,y,z),axis=1)
    
    
if __name__ == "__main__":
   