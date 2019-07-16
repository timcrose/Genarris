from aimsutils.rotate.common import C, C_a, lsplit
from aimsutils import parser, restarts, ylms
import quaternion
import numbers
import numpy as np
import os
from copy import deepcopy
import re
os.environ["NUMBA_DISABLE_JIT"] = "1"
import spherical_functions as sf
from ase import io
print('io.__file__', io.__file__, flush=True)

class Rotator(object):
    def __init__(self, lmax, *args):
        """
        Setup for a single rotation for given maximum l value. Rotation can be
        defined via Euler angles (ZYZ convention) or quaternion
        (using np.quaternion).

        Parameters
        ----------
        lmax : int
            Highest l value to be considered.
        alpha, beta, gamma : float
            Rotation in Euler angles in degree
        quaternion : np.quaternion
            Rotation as quaternion
        """
        if lmax > 9:
            raise NotImplementedError('Only l values of up to "l = 9" \
are currently implemented (aka hardcoded)!')

        self.R = _identify_rotation(*args)
        self.C = C[:lmax+1]
        self.C_a = C_a[:lmax+1]
        self.lmax = lmax
        self.D = []
        self.calculate_wigner_D()

    def calculate_wigner_D(self):
        """
        Calculate the Wigner D matrices for the given rotation.
        """
        # Get Wigner_D matrices
        D = np.split(sf.Wigner_D_matrices(self.R, 0, self.lmax),
                     lsplit[:self.lmax])
        self.D = D
        Ds = []
        Deltas = []
        for l, x in enumerate(D):
            ell = 2*l+1
            DD = x.reshape([ell, ell])
            Ds.append(DD)
            # could check sum of imag via sum(x.imag.flatten())
            Deltas.append(
                (self.C_a[l].conjugate().dot(DD).dot(self.C_a[l].T)).real)
        self.Deltas = Deltas
        self.pmat = (self.C[1].conjugate().dot(Ds[1]).dot(self.C[1].T)).real

    def rotate_eigenvector(self, KS_ev, lvec):
        """
        Do the rotation using the given angles for a set of eigenvectors.

        Parameters
        ----------
        KS_ev: np.array
            The Kohn-Sham eigenvector from the restart file.
        lvec: list
            The full list of basis function quantum numbers.

        Returns
        -------
        KS_ev_rotated: np.array
            The rotated KS eigenvector
        """
        # n_basis = KS_ev.shape[0]
        n_states = KS_ev.shape[1]
        # n_spin = KS_ev.shape[2]
        KS_ev_work = ylms.group_array(KS_ev, lvec)
        KS_ev_rotated = []
        for basis in KS_ev_work:
            rot = np.array(np.zeros(basis.shape))

            l = int((basis.shape[0]-1)/2)
            Delta = self.Deltas[l]
            for i_state in range(n_states):
                rot[:, i_state] = Delta.dot(basis[:, i_state])
            KS_ev_rotated.append(rot)

        KS_flat = ylms.flatten_array(np.array(KS_ev_rotated))
        return KS_flat

    def rotate_structure(self, structure, center='COM', translation=[0, 0, 0]):
        """
        Parameters
        ----------
        structure : ASE atoms object or str
            An ASE atoms object or a path to a ASE-readable file.
        center : 3-tuple or str
            The point to rotate about or a string to select either
            "COM" (center of mass) or "COP" (center of positions).
        translation : [x, y, z], optional
            The vector for a translation of the rotated molecule.

        Returns
        -------
        atoms : ase.atoms.Atoms
            An atoms object with the new coordinates.
        """
        atoms = deepcopy(structure)

        if isinstance(center, str):
            if center.lower() == 'com':
                center = atoms.get_center_of_mass()
            elif center.lower() == 'cop':
                center = atoms.get_positions().mean(axis=0)
            elif center.lower() == 'cou':
                center = atoms.get_cell().sum(axis=0) / 2
            else:
                raise ValueError('Cannot interpret center')
        else:
            center = np.array(center)

        pos = atoms.get_positions() - center

        R = np.roll(self.pmat, 3)
        for num, row in enumerate(R):
            R[num] = np.roll(row, 1)
            pmat = R

        new_pos = (np.dot(pmat, pos.T)).T + center + np.array(translation)
        atoms.set_positions(new_pos)

        return atoms


class AimsCalculation(object):
    def __init__(self, aimsout, basis_out):
        """
        Load and store data from a finished AIMS calculation.
        """
        self.ks_eigenvector, self.ks_eigenvalues,\
            self.occupations = restarts.read_restart(aimsout)
        self.basis = parser.parse_basis(basis_out)
        self.lvec = ylms.get_l_vector(self.basis)
        self.lmax = max(self.lvec)
        self.structure = io.read(os.path.join(os.path.dirname(aimsout),
                                              "geometry.in"), parallel=False)


class Rotations(object):
    def __init__(self, title="Rotation", cell=[None, None, None]):
        """Handle arbitary rotations for FHIaims wave functions.

        This class is used to handle any number of rotations for any
        number of different molecules.

        Rotations can be added individually or by parsing a file
        with instructions.

        """
        self.title = title
        self.cell = cell
        self.unique_calculations = {}

        self.rot_mol = []
        self.rot_ks_ev = []
        self.rot_ks_e = []
        self.rot_ks_o = []

    def __reset__(self):
        """Reset all data."""
        self.unique_calculations = {}

        self.rot_mol = []
        self.rot_ks_ev = []
        self.rot_ks_e = []
        self.rot_ks_o = []

    def __repr__(self):
        n = min(3, len(self.rot_mol))
        rep = "Rotations"
        if n == 0:
            rep += "()"
        else:
            rep += "(\n"
            for x in [str(r) for r in self.rot_mol[:n]]:
                rep += str(x) + "\n"
            rep += "...)"
        return rep

    def add_rotation(self, parent_mol, rotation, translation=[0, 0, 0]):
        """Add a single rotation.

        Parameters
        ----------
        parent_mol:
            Foldername for parent molecule to be rotated
        rotation:
            The rotation to be done (via Euler angle or Quaternion)
        translation:
            The translation of the molecule w.r.t. parent_mol
        """
        if parent_mol not in self.unique_calculations:
            aimsout = _check_master(parent_mol)
            calc = AimsCalculation(aimsout, os.path.join(parent_mol,
                                                         "basis-indices.out"))
            self.unique_calculations[parent_mol] = calc
        else:
            calc = self.unique_calculations[parent_mol]

        if isinstance(rotation, np.quaternion):
            rotation = [rotation]
        rot = Rotator(calc.lmax, _identify_rotation(*rotation))
        rot_ev = rot.rotate_eigenvector(calc.ks_eigenvector, calc.lvec)
        rot_mol = rot.rotate_structure(calc.structure, translation=translation)
        self.rot_ks_ev.append(rot_ev)
        self.rot_mol.append(rot_mol)
        self.rot_ks_e.append(calc.ks_eigenvalues)
        self.rot_ks_o.append(calc.occupations)

    def write_restartfile(self, folder=None, filename="restart.combined", write_geometry=True):
        """
        Write out files for rotated calculation (geometry.in, restart.combined).

        Parameters
        ----------
        folder : str, optional
            The folder for the rotated calculation.
        filename : str, optional
            The filename for the restart file.
            Default is 'restart.combined'
        """
        if folder is None:
            folder = self.title

        if folder!='.' and folder!='./':
            try:
                os.mkdir(folder)
            except OSError:
                raise

        super_ks_ev, super_ks_e, super_occ = \
            restarts.combine_wavefunctions(self.rot_ks_ev,
                                           self.rot_ks_e,
                                           self.rot_ks_o)

        supercell = restarts.combine_geometries(self.rot_mol, self.cell)
        if all(supercell.get_pbc()):
            periodic = True
        else:
            periodic = False

        restarts.write_restart(os.path.join(folder, filename),
                               periodic,
                               super_ks_ev,
                               super_ks_e,
                               super_occ)

        if write_geometry:
            io.write(os.path.join(folder, "geometry.in"),
                 supercell,
                 format="aims")

    def load_rotations(self, path="rotations.in"):
        """Read a set of rotations from a file.

        Rotations can be defined via Euler angles (alpha, beta, gamma) or
        Quaternions (x0, x1, x2, x3). Mixed definitions are not allowed.

        Example rotations.in:
        .. code-block::
            :linenos:

            # This is system ABC
            # Project X
            title project_x
            lattice_vector 11.2 0.00 0.00
            lattice_vector 0.00 10.5 0.00
            lattice_vector 0.00 0.00 10.7
            # now all rotated molecules
            # alpha beta gamma delta_x delta_y delta_z molecule
            33.4 12.4 167.0 10. 0. 0. h2_ref
            0.   45.  0.    3.5 2.1 5.4 h2_ref
            130. 90. 23.    0.  2.6 12. o2_ref
            # Note that each instruction line is assigned to a base
            # molecule (via m1 or m2)
        """
        print('inside load_rotations', flush=True)
        self.__reset__()
        rotations = []
        cell = []
        with open(path) as f:
            data = f.readlines()
        for line in data:
            if not bool(re.search(r"\s*#|^\s*$", line)):
                if "lattice_vector" in line:
                    cell.append([float(x) for x in line.split()[1:]])
                elif "title" in line:
                    title = line.split()[-1].strip()
                else:
                    if len(line.split()) == 7:
                        # euler case
                        rotation = [float(x) for x in line.split()[0:3]]
                        translation = np.array([float(x) for x
                                                in line.split()[3:6]])
                    elif len(line.split()) == 8:
                        # quaternion case
                        rotation = [float(x) for x in line.split()[0:4]]
                        translation = np.array([float(x) for x
                                                in line.split()[4:7]])

                    molecule = [x for x in line.split()][-1]

                    rotations.append([x for x in
                                     [rotation, translation, molecule]
                                     if x is not None])
            else:
                continue

        if cell == []:
            cell = [None, None, None]

        self.cell = cell
        self.title = title
        for rot, trans, mol in rotations:
            self.add_rotation(mol, rot, trans)
        print('leaving load_rotations', flush=True)


def _identify_rotation(*args):
    if isinstance(args[0], np.quaternion):
        R = args[0]
    elif isinstance(args[0], numbers.Number) \
        and isinstance(args[1], numbers.Number) \
            and isinstance(args[2], numbers.Number):
        R = quaternion.from_euler_angles(np.radians(-args[0]),
                                         np.radians(args[1]),
                                         np.radians(-args[2]))
    else:
        raise ValueError("Can't understand input rotation")
    return R


def _check_master(master):
    """
    Check if the given path is a valid aims output folder with all
    necessary files for rotation.

    Parameters
    ----------
    master : str
        Path to check.

    Returns
    -------
    master : str
        Path to the AIMS output file.
    """
    if os.path.isdir(master):
        p = "Have a nice day"
        fnames = [os.path.join(master, fn) for fn in next(os.walk(master))[2]]

        for file in fnames:
            try:
                with open(file, "r") as f:
                    lines = f.readlines()
                for line in lines:
                    if p in line:
                        aims_out = file
                        break
            except:
                pass

        if not aims_out:
            raise RuntimeError("""Could not find a converged FHIaims
                               calculation.""")
        if not os.path.isfile(os.path.join(master, "basis-indices.out")):
            raise RuntimeError("""Could not find a valid basis-indices.out file
in folder. Please make sure to use the flag
'output h_s_matrices' in your control.in""")

    return aims_out
