# -*- coding: utf-8 -*-

import numpy as np
from scipy import sparse

from ase.neighborlist import NeighborList,natural_cutoffs,get_distance_indices

from .inter_dist_values import inter_dist, \
                               construct_pair_keys


class MoleculeBonding():
    """
    Identifies all bonded atoms in a molecule or structure using ASE NeighborList
    
    """
    def __init__(self, struct):
        self.struct = struct
        self.ele = struct.geometry["element"]
        self.bonding_list = self._get_bonding(struct)
        
        # Donor and acceptor elements for identifying hydrogen bonds
        self.donor_elements = ["N", "O"]
        self.acceptor_elements = ["N", "O"]
        
    
    def get_cutoff_matrix(self,vdw_mult=0.85):
        """
        Returns NxN matrix expressing a cutoff distance between intermolecular
        contacts in a crystal system of the given molecule.
        
        Arguments
        ---------
        vdw_mult: float
            Multiplicative factor for the cutoff distance for vdw type contacts. 
            A cutoff value of 0.85 is well supported by statistical analysis 
            for these types of contacts.
        """
        # Get pair_key_matrix and construct initial vdw cutoff
        pair_key_matrix = construct_pair_keys(self.struct.geometry["element"])
        cutoff_matrix = np.zeros(pair_key_matrix.shape)
        for index,value in np.ndenumerate(pair_key_matrix):
            cutoff_matrix[index] = inter_dist["vdw"][value]
        cutoff_matrix = cutoff_matrix*vdw_mult
        
        # Add hydrogen bond cutoff distances
        donor_idx,acceptor_idx = self._get_hydrogen_bond_idx()
        if len(donor_idx) > 0:
            acceptor_ele = self.ele[acceptor_idx]
            donor_ele = self.ele[np.concatenate(self.bonding_list[donor_idx])]
            hbond_key = np.char.add(donor_ele,"H-")
            hbond_key = np.char.add(hbond_key[:,None],acceptor_ele)
            
            hbond_values = np.zeros(hbond_key.shape)
            for index,value in np.ndenumerate(hbond_key):
                hbond_values[index] = inter_dist["h_bond"][value]
                
            # Create pairwise index grid for indexing into cutoff_matrix
            ixgrid1 = np.ix_(donor_idx,acceptor_idx)
            ixgrid2 = np.ix_(acceptor_idx,donor_idx)
            cutoff_matrix[ixgrid1] = hbond_values
            cutoff_matrix[ixgrid2] = hbond_values.T
            
        return cutoff_matrix
    
    
    def get_cutoff_matrix_vdw(self,vdw_mult=0.85):
        """
        Returns NxN matrix expressing a cutoff distance between intermolecular
        contacts in a crystal system of the given molecule. Only uses vdw 
        contact distances.
        
        Arguments
        ---------
        vdw_mult: float
            Multiplicative factor for the cutoff distance for vdw type contacts. 
            A cutoff value of 0.85 is well supported by statistical analysis 
            for these types of contacts.
        """
        # Get pair_key_matrix and construct initial vdw cutoff
        pair_key_matrix = construct_pair_keys(self.struct.geometry["element"])
        cutoff_matrix = np.zeros(pair_key_matrix.shape)
        for index,value in np.ndenumerate(pair_key_matrix):
            cutoff_matrix[index] = inter_dist["vdw"][value]
        cutoff_matrix = cutoff_matrix*0.85
        
        return cutoff_matrix
    
    
    def get_crystal_cutoff_matrix(self,nmpc,vdw_mult=0.85):
        """
        Copies the intermolecular distance matrix from 
        MoleculeBonding.get_cutoff_matrix into the correct size for a specific
        number of molecules in the unit cell.
        """
        cutoff_matrix = self.get_cutoff_matrix(vdw_mult=vdw_mult)
        return np.tile(cutoff_matrix,(nmpc,nmpc))
        
    
    def _get_bonding(self,struct):
        # Use ASE neighborlist class to identify bonding of atoms
        atoms = struct.get_ase_atoms()
        cutOff = natural_cutoffs(atoms)
        neighborList = NeighborList(cutOff, self_interaction=False,
                                            bothways=True)
        neighborList.update(atoms)
        
        # Construct bonding list indexed by atom in struct
        bonding_list = [[] for x in range(struct.geometry.shape[0])]
        for i in range(struct.geometry.shape[0]):
            bonding_list[i] = neighborList.get_neighbors(i)[0]
        
        return np.array(bonding_list)
    
    
    def _get_hydrogen_bond_idx(self):
        """
        For all hydrogen atoms in the system, identify the index of any which 
        can participate in hydrogen. This is defined as hydrogens which are 
        bonded to polar elements: nitrogen, oxygen. Returns a list of 
        indices for donors and acceptors of hydrogen bonds
        """
        self.donor_idx = []
        self.acceptor_idx = []
        
        for i,ele in enumerate(self.struct.geometry["element"]):
            # Donors definition
            if ele == "H":
                bond_idx = self.bonding_list[i]
                elements = self.struct.geometry["element"][bond_idx]
                for h_ele in self.acceptor_elements:
                    if h_ele in elements:
                        self.donor_idx.append(i)
                        break
            # Acceptor definition
            elif ele in self.acceptor_elements:
                bonding = self.bonding_list[i]
                # Check for terminal oxygen or bridging oxygen
                if ele == "O":
                    if len(bonding) == 1:
                        self.acceptor_idx.append(i)
                    elif len(bonding) == 2:
                        bond_ele = self.ele[bonding]
                        unique_ele = np.unique(bond_ele)
                        if len(unique_ele) == 1 and unique_ele[0] == "C":
                            self.acceptor_idx.append(i)
                        
                # Check for terminal nitrogen
                if ele == "N":
                    if len(bonding) <= 2:
                        self.acceptor_idx.append(i)
        
        return self.donor_idx,self.acceptor_idx



if __name__ == "__main__":
    pass
        
        
