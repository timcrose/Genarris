# -*- coding: utf-8 -*-

__author__ = 'Manny Bier'

from ibslib.molecules.find_bonding import MoleculeBonding
from ibslib.molecules.label_molecule_atoms import LabelAtoms
from ibslib.molecules.identify_single_molecule import UniqueMolecules
from ibslib.molecules.volume_estimator import MoleculeVolumeEstimator
__all__ = ["MoleculeBonding","UniqueMolecules","MoleculeVolumeEstimator",
           "LabelAtoms"]

def check_molecule(struct):
    # Check for valid molecule_struct
    if len(struct.get_lattice_vectors_better()) > 0:
        raise Exception("Structure with lattice vectors {} was passed "
                .format(struct.get_lattice_vectors_better())+
                "into MoleculeBonding class. Molecule structure "+
                "without lattice vectors should be passed to "+
                "MoleculeBonding.")
    else:
        return