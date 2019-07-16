# -*- coding: utf-8 -*-



from scipy import sparse
from ase.neighborlist import NeighborList,natural_cutoffs,get_distance_indices


class MoleculeBonding():
    """
    Identifies all bonded atoms in a molecule or structure using ASE NeighborList
    
    """
    def __init__(self, struct):
        self.struct = struct
        self.bonding_list = self._get_bonding(struct)
        
    
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
        
        return bonding_list
    
    
    def get_hydrogen_bonded_index(self):
        """
        For all hydrogen atoms in the system, identify the index of any which 
        can participate in hydrogen. This is defined as hydrogens which are 
        bonded to polar elements: nitrogen, oxygen, or sulfur. 
        """
        self.h_bonded_idx = []
        h_bonded_elements = ["N", "O"]
        for i,ele in enumerate(self.struct.geometry["element"]):
            if ele != "H":
                continue
            bond_idx = self.bonding_list[i]
            elements = self.struct.geometry["element"][bond_idx]
            for h_ele in h_bonded_elements:
                if h_ele in elements:
                    self.h_bonded_idx.append(i)
                    break
        
        return self.h_bonded_idx



if __name__ == "__main__":
    pass
    