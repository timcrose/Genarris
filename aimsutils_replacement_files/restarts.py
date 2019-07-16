# This file is part of aimsutils.
# (C) 2015 Christoph Schober
import os
from aimsutils.parser import parse_aimsout
import restartutils
import numpy as np
import scipy.linalg


def read_restart(outfile):
    """
    Read a FHIaims restart file in binary format.

    Parameters
    ----------
    outfile : str
        FHIaims output file.

    Returns
    -------
    KS_eigenvector :
        The KS eigenvector
    KS_eigenvalue :
        The KS eigenvalues
    occupations :
        The occupations.
    """
    base, aimsout = os.path.split(outfile)
    meta = parse_aimsout(outfile)
    #print('meta',meta)
    if base != "":
        meta["restartfile"] = os.path.join(base, meta["restartfile"])

    if meta["periodic"]:
        ks_ev, ks_e, occ = \
            restartutils.read_restart_periodic(meta["n_states"],
                                               meta["n_basis"],
                                               meta["n_spin"],
                                               meta["restartfile"])

        KS_eigenvector = ks_ev.reshape(meta["n_basis"],
                                       meta["n_states"],
                                       meta["n_spin"])

    else:
        ks_ev, ks_e, occ = \
            restartutils.read_restart_cluster(meta["n_states"],
                                              meta["n_basis"],
                                              meta["n_spin"],
                                              meta["restartfile"])

        KS_eigenvector = ks_ev.reshape(meta["n_basis"],
                                       meta["n_states"],
                                       meta["n_spin"])

    KS_eigenvalue = list()
    for item in ks_e:
        KS_eigenvalue.append([float(x) for x in item])
    KS_eigenvalue = np.array(KS_eigenvalue)

    occupations = list()
    for item in occ:
        occupations.append([float(x) for x in item])
    occupations = np.array(occupations)

    return KS_eigenvector, KS_eigenvalue, occupations


def write_restart(filename, periodic,
                  KS_eigenvector, KS_eigenvalue, occupations):
    """
    Write a FHIaims restart file in binary format.

    Parameters
    ----------
    filename : str
        The filename for the restart file.
    periodic : bool
        If True, calculation is with periodic boundary conditions
    KS_eigenvector : np.array
        The KS_eigenvector of the system
    KS_eigenvalue : np.array
        The KS_eigenvalues of the system
    occupations : np.array
        The occupations of the system
    """
    n_states = KS_eigenvector.shape[1]
    n_basis = KS_eigenvector.shape[0]
    n_spin = KS_eigenvector.shape[2]
    n_k = 1
    if periodic:
        # TODO careful, hack to work with gamma periodic, need to
        #     change that to proper n_k loop for Martin ;-)
        KS_eigenvector = KS_eigenvector.reshape(n_basis, n_states, n_spin, n_k)
        filename = filename
        restartutils.write_restart_periodic(filename,
                                            KS_eigenvector,
                                            KS_eigenvalue,
                                            occupations)
    else:
        KS_eigenvector = KS_eigenvector.reshape(n_basis, n_states, n_spin)
        restartutils.write_restart_cluster(filename,
                                           KS_eigenvector,
                                           KS_eigenvalue,
                                           occupations)


def combine_wavefunctions(KS_eigenvectors, KS_eigenvalue, occupations):
    """
    Combine any number of FHIaims wavefunctions / molecules to a
    single wavefunction.

    Parameters
    ----------
    KS_eigenvectors : list of np.array
        List of arrays of KS_eigenvectors for each molecule.
    KS_eigenvalue : list of np.array
        List of arrays of KS_eigenvalues for each molecule.
    occupations : list of np.array
        List of arrays of occupations for each molecule.

    Returns
    -------
    super_ks_ev : np.array
        New KS_eigenvector from superimposed single molecule KS_ev.
    super_ks_e : np.array
        New KS_eigenvalue vector.
    super_occ : np.array
        New occupations vector.
    """
    n_spin = KS_eigenvalue[0].shape[1]
    super_ev = [None]*n_spin
    super_e = [None]*n_spin
    super_o = [None]*n_spin

    for i in range(n_spin):
        new_ev = list()
        new_occ = list()
        for mol in KS_eigenvalue:
            new_ev.extend(mol[:,i].tolist())
        for mol in occupations:
            new_occ.extend(mol[:,i].tolist())

        new_ev = np.array(new_ev)
        new_occ = np.array(new_occ)

        old2new = np.argsort(new_ev)
        super_ks_e = new_ev[old2new]
        super_occ = new_occ[old2new]

        # For each other set of KS_ev we need a block of zeros
        # Example system with 3 molecules:
        # x 0 0
        # 0 x 0
        # 0 0 x
        #
        # The order of the rows (=basis) must be identical with the position
        # of the molecules in the geometry.in file.
        KS_ev_work = [ev[:,:,i] for ev in KS_eigenvectors]
        super_ks_ev = scipy.linalg.block_diag(*KS_ev_work)

        # The order of the columns (=states) must be identical with the ordering
        # of the eigenvalues.
        # Degenerated states should not matter.
        super_ks_ev = super_ks_ev[:, old2new]

        super_ks_ev = np.expand_dims(super_ks_ev, axis=2)
        super_ks_e = np.expand_dims(super_ks_e, axis=1)
        super_occ = np.expand_dims(super_occ, axis=1)

        super_ev[i] = np.copy(super_ks_ev)
        super_e[i] = np.copy(super_ks_e)
        super_o[i] = np.copy(super_occ)

    super_ev = np.concatenate(super_ev, axis=2)
    super_e = np.concatenate(super_e, axis=1)
    super_o = np.concatenate(super_o, axis=1)

    return super_ev, super_e, super_o


def combine_geometries(geometries, cell=None):
    """
    Construct a supercell from all molecules in geometries.

    Parameters
    ----------
    geometries : list
        List of ase.atoms.Atoms-objects.
    cell : list, optional
        Lattice vectors for periodic supercell.

    Returns
    -------
    supercell : ase.atoms.Atoms
        A ase.atoms.Atoms object with all molecules.
    """
    supercell = geometries[0]
    for molecule in geometries[1:]:
        supercell += molecule

    if all(cell):
        supercell.set_pbc(True)
        supercell.set_cell(cell)

    return supercell
