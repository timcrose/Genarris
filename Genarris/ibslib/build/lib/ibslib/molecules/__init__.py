# -*- coding: utf-8 -*-

__author__ = 'Manny Bier'

from .volume_estimator import MoleculeVolumeEstimator

__all__ = ["MoleculeVolumeEstimator"]

def check_molecule(struct):
    # Check for valid molecule_struct
    if len(struct.get_lattice_vectors_better()) > 0:
        raise Exception("Structure with lattice vectors {} was passed "
                .format(struct.get_lattice_vectors_better())+
                "into MoleculeBonding class. Molecule structure "+
                "without lattice vectors should be passed to "+
                "MoleculeBonding.")
    else:
        return True
