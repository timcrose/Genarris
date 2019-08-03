
import os
import torch
import numpy as np

from ibslib import Structure
from ibslib.io import read,write
from ibslib.descriptor.R_descriptor import calc_R,ele_R,struct_R
from ibslib.motif.utils import get_molecules,reconstruct_with_whole_molecules
from ibslib.molecules import LabelAtoms


class UniqueMolecules():
    """
    Class for identifying unique molecules in a structure
    """
    def __init__(self, struct, mult=1):
        """
        For now, struct is argument for class and will automatically
        identify molecules and save as parameter of class.
        
        Arguments
        ---------
        mult: float
            Multiplicative factor to multiply ase.neighborlist.natural_cutoff.
        """
        self.struct = struct
        # Reconstruct with whole molecules because it will make all further 
        # calculations faster
        self.molecule_struct_list = get_molecules(self.struct, mult=mult)
        self.unique_molecule_list = find_unique_molecules(self.molecule_struct_list)
        self.rstruct = self._reconstruct_structure()
        
    
    def write_unique(self,folder,file_format="json",overwrite=False):
        """
        Writes unique single molecules found for the structure to a single 
        folder.
        """
        struct_id = self.struct.struct_id
        volume = self.struct.get_unit_cell_volume() 
        molecule_volume = volume / len(self.molecule_struct_list)
        
        for i,molecule_struct in enumerate(self.unique_molecule_list):
            if len(self.unique_molecule_list) == 1:
                # Add molecular volume if there's only one type of molecule
                molecule_struct.set_property("molecule_volume", 
                                             molecule_volume)
            molecule_name = struct_id + "_molecule_{}".format(i)
            path = os.path.join(folder,molecule_name)
            write(path, molecule_struct, file_format=file_format,
                  overwrite=overwrite)
           
    
    def write_rstruct(self,folder,file_format="json",overwrite=False):
        """
        Writes the reconstructed geometry from the smallest molecule
        representation.
        """
        struct_id = self.struct.struct_id
        self.rstruct.set_property("num_molecules", 
                                  len(self.molecule_struct_list))
        self.rstruct.set_property("num_unique_molecules", 
                                  len(self.unique_molecule_list))
        path = os.path.join(folder,struct_id)
        write(path,self.rstruct,file_format=file_format,
              overwrite=overwrite)
        
    
    def _reconstruct_structure(self):
        """
        Reconstructs the system using the smallest molecule representation.
        If there's only one unique molecule, then reorders the molecules to 
        all have the same atom indexing.

        """
        rstruct = Structure()
        rstruct.set_lattice_vectors(self.struct.get_lattice_vectors())
        if len(self.unique_molecule_list) == 1:
            self._reorder_atoms()
        for molecule_struct in self.molecule_struct_list:
            geo_array = molecule_struct.get_geo_array()
            ele = molecule_struct.geometry['element']
            for i,coord in enumerate(geo_array):
                rstruct.append(coord[0],coord[1],coord[2],ele[i])
        return rstruct
    
    
    def _reorder_atoms(self):
        """
        Reorders atoms in each molecule based off the ordering of the unique 
        molecule using LabelAtoms from ibslib using Pymatgen functions.
        """
        unique_molecule = self.unique_molecule_list[0]
        la  = LabelAtoms()
        for molecule in self.molecule_struct_list:
            la.uniform_labels(unique_molecule, molecule, 
                              hold_first_order=True)
        


def find_unique_molecules(molecule_struct_list, tol=1):
    """
    Bad wrapper for now
    """
    molecule_groups = find_unique_groups(molecule_struct_list,tol=1)
    unique_molecule_idx = [x[0] for x in molecule_groups]
    return [molecule_struct_list[x] for x in unique_molecule_idx]


def reconstruct_from_unique_molecules(molecule_struct_list, tol=1):
    """
    Reconstruct the geometry from unique molecules identified in the structure
    """
    molecule_groups = find_unique_groups(molecule_struct_list,tol=1)
    unique_molecule_idx = [x[0] for x in molecule_groups]
    print("Not finished yet")


def find_unique_groups(molecule_struct_list, tol=1):
    """
    Returns a list of indices for unique molecule groups where every molecule 
    in each group are indentical within the tolerance.
    
    Arguments
    ---------
    struct: Structure
    tol: float
        Tolerance for identifying individual molecules
        
    """   
    difference_matrix = p_molecule_distance(molecule_struct_list)
    molecule_groups = unique_groups(difference_matrix, tol=tol)
    return molecule_groups
    

def unique_groups(difference_matrix, tol=1):
    """ 
    Breaks difference matrix into groups of similar molecules.
    Returns list molecules which make up unique groups for indexing into the 
    original molecule_struct_list.
    """
    # List of molecules which need to be sorted
    unsorted = np.array([x for x in range(difference_matrix.shape[0])])
    # List of groups of the same molecules
    molecule_groups = []
    
    while len(unsorted) != 0:
        # Pick first molecule which hasn't been sorted
        current_idx = unsorted[0]
        # Take its row of the difference matrix
        row = difference_matrix[current_idx,:]
        # Take positions of other molecules yet to be sorted
        row = row[unsorted]
        # Find those same by tolerance value
        same_tol = row < tol
        # Must be greater than 0 so value cannot be -1
        same_nonzero = row >= 0
        # Combine
        same = np.logical_and(same_tol, same_nonzero)
        # Gather results
        same_idx = np.where(same == True)[0]
        # Reference original molecule index which is obtained from the 
        # unsorted list for the final unique groups
        molecule_groups.append(unsorted[same_idx])
        # Delete indexs in unsorted which have now been sorted
        unsorted = np.delete(unsorted, same_idx)              
    
    return molecule_groups


def p_molecule_distance(molecule_list):
    """ 
    Finds minimum pairwise distance between all molecules in molecule list.
    Pairwise distance computed by taking the difference between the internal 
    coordinates of each molecule in the molecule list.
    
    Returns -1 if the elements or number of elements are not the same.
    Otherwise, returns positive value of the difference between internal 
    coordinates of the two different molecules. 
    
    """
    num_mol = len(molecule_list)
    
    # Initialize to zeros
    difference_matrix = np.zeros((num_mol,num_mol))
    
    # Loop over upper diagonal of pairwise matrix
    for i in range(num_mol):
        for j in range(i+1,num_mol):
            molecule_0 = molecule_list[i]
            molecule_1 = molecule_list[j]
            
            # Compute internal coordinates of 1
            R_0 = calc_R(molecule_0)
            R_0 = ele_R(R_0)
            R_0 = struct_R(R_0)
            
            # Compute internal coordinates of 2
            R_1 = calc_R(molecule_1)
            R_1 = ele_R(R_1)
            R_1 = struct_R(R_1)
            
            # Check if elements are the same
            if [key for key in R_0.keys()] != [key for key in R_1.keys()]:
                 # If elements are not all the same then the structures
                 # are certainly different so return -1
                 difference_matrix[i,j] = -1
                 continue
            
            # Computer difference w.r.t. different interactions in system
            difference = 0
            diff = 0
            for element,inter_dict in R_0.items():
                for inter,value_tensor in inter_dict.items():
                    R_1_value_tensor = R_1[element][inter]
                    
                    # Check if number of interactions are the same
                    if len(value_tensor) != len(R_1_value_tensor):
                         # If iteractions are not the same length then they
                         # are certainly different so return -1
                         diff = -1
                         continue
                    
                    temp = torch.sum(value_tensor - R_1_value_tensor)
                    difference += torch.abs(temp) / len(value_tensor.view(-1))
            
            # Put value in difference matrix symmetric across diagonal
            if diff < 0:
                difference_matrix[i,j] = -1
                difference_matrix[j,i] = -1
            else:
                difference_matrix[i,j] = difference
                difference_matrix[j,i] = difference
    
    return difference_matrix


if __name__ == "__main__":
    pass