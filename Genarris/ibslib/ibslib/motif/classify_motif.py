
from ibslib import Structure,StructDict
from ibslib.motif.utils import motif_definitions, \
                               construct_orientation_supercell, \
                               compute_orientation_difference 


class MotifClassifier():
    def __init__(self, molecule_finder, supercell=3, include_negative=True, 
                 num_mol=12):
        """

        Arguments
        ---------
        molecule_finder: callable function
            Function which will perform the task of identifying the molecules
            in each structure. Can be general algorithm or can use the 
            index of each atom if the molecules are indexed systematically.
            For example, use ibslib.motif.get_molecules for general case or
            ibslib.motif.molecule_by_index for indexed case. 
        supercell: int
            Value of the supercell dimension. Example 3x3x3
        include_negative: bool
            False: Only supercells in the positive direction will be constructed
            True: Supercells in the positive and negative direction will be 
                    constructed. This will double the number constructed. 
         num_mol: int >= 4
            Number of nearest neighbor molecules to be used for motif 
            identification. Should be at least four.

        """
        self.molecule_finder = molecule_finder
        self.supercell = supercell
        self.include_negative = include_negative
        self.num_mol = num_mol

    
    def calc(self, struct_obj):
        """
        Calculate motifs for either a Structure or StructDict
        """
        obj_type = type(struct_obj)
        if obj_type == dict or obj_type == StructDict:
            return self._calc_dict(struct_obj)
        elif obj_type == Structure:
            return self._calc_struct(struct_obj)
    

    def _calc_struct(self, struct):
        """
        Calculates motif for a single molecule
        """
        # Get molecules
        molecule_struct_list = self.molecule_finder(struct)
        # Construct orientations and COM positions in supercell
        orientation_tensor,COM_array = construct_orientation_supercell(struct, 
                                        self.supercell,self.include_negative,
                                        molecule_struct_list)
        # Compute orientation difference from a central molecule
        deg_array,plane_deg_min = compute_orientation_difference(orientation_tensor,
                                        COM_array,molecule_struct_list,
                                        num_mol=self.num_mol)
        # Use these differences to identify the motif
        return motif_definitions(deg_array,plane_deg_min)
        

    def _calc_dict(self, struct_dict):
        motif_list = []
        for struct_id,struct in struct_dict.items():
            motif_list.append(self._calc_struct(struct))
        return motif_list
    

    
