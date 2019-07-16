'''
Created on June 17, 2015
Author: Patrick Kilecdi
'''
import sys

import random, copy
random.seed()
import numpy as np
np.random.seed()
from Genarris.core import structure_handling
from Genarris.core.structure import Structure
from Genarris.generation import sgroup
from Genarris.utilities import misc, write_log
from Genarris.utilities.misc import half_gaussian_sampling_upper, \
        half_gaussian_sampling_lower, random_rotation_matrix, \
        frac_coor_adjust, output_structure, \
        retrieve_integer_and_subtract
from Genarris.utilities.check_type import *
from Genarris.utilities.write_log import print_time_log
import time

def get_structure_generator_from_inst(inst, sname):
    ucv_target = inst.get_eval(sname, "ucv_target")
    nmpc = inst.get_eval(sname, "NMPC")
    is_chiral = inst.get_boolean(sname, "is_chiral")
    is_racemic = inst.get_boolean(sname, "is_racemic")
    sname_list = [sname, 'relax_single_molecule', 'estimate_unit_cell_volume', 'harris_single_molecule_prep', 'pygenarris_structure_generation', 'structure_generation_batch', 'harris_approximation_batch']
    molecule_path = inst.get_inferred(sname, sname_list, ['molecule_path'] * 7, type_='file')
    molecule_name = inst.get_with_default(
            sname, "molecule_name", "molecule_original", eval=False)
    enantiomer_name = inst.get_with_default(
            sname, "enantiomer_name", "molecule_enantiomer", eval=False)
    space_groups_allowed = inst.get_or_none(
            sname, "space_groups_allowed", eval=True)
    wyckoff_list = inst.get_with_default(
            sname, "wyckoff_list", [0], eval=True)
    ucv_ratio_range = inst.get_with_default(
            sname, "ucv_ratio_range", [1,1], eval=True)
    ucv_std = inst.get_or_none(sname, "ucv_std", eval=True)
    p_tolerance = inst.get_with_default(
            sname, "p_tolerance", 0.25, eval=True)
    angle_range = inst.get_with_default(
            sname, "angle_range", [30,150], eval=True)
    ax_variance = inst.get_with_default(
            sname, "ax_variance", [0.5,2], eval=True)
    by_variance = inst.get_with_default(
            sname, "by_variance", [0.5,2], eval=True)
    cz_variance = inst.get_with_default(
            sname, "cz_variance", [0.5,2], eval=True)
    com_dist = inst.get_or_none(sname, "com_dist", eval=True)
    atom_dist = inst.get_or_none(sname, "atom_dist", eval=True)
    specific_radius_proportion = inst.get_or_none(
            sname, "specific_radius_proportion", eval=True)
    specific_radius_std = inst.get_or_none(
            sname, "specific_radius_std", eval=True)
    specific_radius_lower_bound = inst.get_or_none(
            sname, "specific_radius_lower_bound", eval=True)
    custom_radii = inst.get_with_default(
            sname, "custom_radii", {}, eval=True)
    custom_radii_by_pairs = inst.get_with_default(
            sname, "custom_radii_by_pairs", {}, eval=True)
    attempt_limit = inst.get_with_default(
            sname, "attempt_limit", 1024, eval=True)
    attempts_per_space_group_change = inst.get_with_default(
            sname, "attempts_per_space_group_change", 64, eval=True)
    attempts_per_unit_cell_change = inst.get_with_default(
            sname, "attempts_per_unit_cell_change", 16, eval=True)
    attempts_per_wyckoff_list_change = inst.get_with_default(
            sname, "attempts_per_wyckoff_list_change", 4, eval=True)
    attempts_per_closeness_criteria_change = inst.get_with_default(
            sname, "attempts_per_closeness_criteria_change", 1, eval=True)
    fill_molecule_attempt_limit = inst.get_with_default(
            sname, "fill_molecule_attempt_limit", 5, eval=True)
    output_dir = inst.get_or_none(sname, "output_dir")
    output_format = inst.get_with_default(
            sname, "output_format", "json")
    struct_id_scheme = inst.get_with_default(sname, "struct_id_scheme",
            ["struct_index","_","random_index"], eval=True)
    struct_index_length = inst.get_with_default(
            sname, "struct_index_length", 0, eval=True)
    index_tracking_file = inst.get_or_none(sname, "index_tracking_file")
    attempt_tracking_file = inst.get_or_none(sname, "attempt_tracking_file")

    check_float(ucv_target, "ucv_target")
    check_int(nmpc, "nmpc")
    check_type_or_none(
            space_groups_allowed, "space_groups_allowed", (int, list))
    check_list(wyckoff_list, "wyckoff_list")
    check_list(ucv_ratio_range, "ucv_ratio_range")
    check_float_or_none(ucv_std, "ucv_std")
    check_float(p_tolerance, "p_tolerance")
    check_list(angle_range, "angle_range")
    check_list(ax_variance, "ax_variance")
    check_list(by_variance, "by_variance")
    check_list(cz_variance, "cz_variance")
    check_float_or_none(com_dist, "com_dist")
    check_float_or_none(com_dist, "atom_dist")
    check_float_or_none(specific_radius_proportion, "specific_radius_proportion")
    check_float_or_none(specific_radius_std, "specific_radius_std")
    check_float_or_none(
            specific_radius_lower_bound, "specific_radius_lower_bound")
    check_dict(custom_radii, "custom_radii")
    check_dict(custom_radii_by_pairs, "custom_radii_by_pairs")
    check_int(attempt_limit, "attempt_limit")
    check_int(attempts_per_space_group_change,
              "attempts_per_space_group_change")
    check_int(attempts_per_unit_cell_change, "attempts_per_unit_cell_change")
    check_int(attempts_per_wyckoff_list_change,
              "attempts_per_wyckoff_list_change")
    check_int(attempts_per_closeness_criteria_change,
              "attempts_per_closeness_criteria_change")
    check_int(fill_molecule_attempt_limit,
              "fill_molecule_attempt_limit")
    check_list(struct_id_scheme, "struct_id_scheme")
    check_int(struct_index_length, "struct_index_length")
    
    print_time_log(
            "Structure generator initialized from section: " + sname)

    if is_chiral and is_racemic:
        raise ValueError("Specify only one of is_chiral and is_racemic.")

    return StructureGenerator(
            ucv_target, nmpc, is_chiral, is_racemic, 
            molecule_path=molecule_path,
            molecule_name=molecule_name, enantiomer_name=enantiomer_name,
            space_groups_allowed=space_groups_allowed,
            wyckoff_list=wyckoff_list,
            ucv_ratio_range=ucv_ratio_range, ucv_std=ucv_std,
            p_tolerance=p_tolerance, angle_range=angle_range,
            ax_variance=ax_variance, by_variance=by_variance,
            cz_variance=cz_variance, com_dist=com_dist,
            atom_dist=atom_dist,
            specific_radius_proportion=specific_radius_proportion,
            specific_radius_std=specific_radius_std,
            specific_radius_lower_bound=specific_radius_lower_bound,
            custom_radii=custom_radii,
            custom_radii_by_pairs=custom_radii_by_pairs,
            attempt_limit=attempt_limit,
            attempts_per_space_group_change=attempts_per_space_group_change,
            attempts_per_unit_cell_change=attempts_per_unit_cell_change,
            attempts_per_wyckoff_list_change=attempts_per_wyckoff_list_change,
            attempts_per_closeness_criteria_change=
            attempts_per_closeness_criteria_change,
            fill_molecule_attempt_limit=fill_molecule_attempt_limit,
            output_dir=output_dir, output_format=output_format,
            struct_id_scheme=struct_id_scheme,
            struct_index_length=struct_index_length,
            index_tracking_file=index_tracking_file,
            attempt_tracking_file=attempt_tracking_file)

class StructureGenerator():
    '''
    This is the master class behind structure generation
    '''
    def __init__(
        self, ucv_target, nmpc, is_chiral, is_racemic, 
        molecule=None, molecule_path=None,
        molecule_name="molecule_original", enantiomer_name="molecule_enantiomer",
        space_groups_allowed=None, wyckoff_list=[0],
        ucv_ratio_range=[1,1], ucv_std=None,
        p_tolerance=0.25, angle_range=[30,150],
        ax_variance=[0.5,2], by_variance=[0.5,2], cz_variance=[0.5,2],
        com_dist=None, atom_dist=None,
        specific_radius_proportion=None,
        specific_radius_std=None, specific_radius_lower_bound=None,
        custom_radii={}, custom_radii_by_pairs={},
        attempt_limit=1024,
        attempts_per_space_group_change=64,
        attempts_per_unit_cell_change=16,
        attempts_per_wyckoff_list_change=4,
        attempts_per_closeness_criteria_change=1,
        fill_molecule_attempt_limit=5,
        output_dir=None, output_format="json",
        struct_id_scheme=["struct_index","_","random_index"],
        struct_index_length=0, index_tracking_file=None,
        attempt_tracking_file=None):
        '''
        ucv: unit cell volume target; actual ucv can be varied based on
            ucv_ratio_range and ucv_std
        nmpc: number of molecule per cell
        is_chiral: whether the generated structure is chiral (or racemic)
        molecule: Structure() object of molecule
        molecule_path: path to molecule file
        molecule_format: format of molecule file
        moleulce_name: molecule name recorded for Harris Approximationi
        enantiomer_name: enantiomer name recorded for Harris Approximation
        space_groups_allowed: A list of allowed space groups; will be trimmed to
            only have valid space groups according to nmpc and is_chiral; empty
            list indicates that all valid space groups are allowed
        wyckoff_list: If not None, forces the Wyckoff positions chosen;
            Default of [0] forces the molecule to always be placed on the
            general Wyckoff position.
        ucv_ratio_range: allowed unit cell ratio range; list of two floats
        ucv_std: If not None, enables half-normal distribution (right) to bias
            towards lower volume
        p_tolerance: Legacy threshold to determine how skewed a unit cell
            can be
        angle_range: List of two int/floats determining range for angles
            that are not fixed as 90 by the Bravais system
        ax_variance: How much the x component of lattice vector a can deviate from
            the cubic root of ucv in ratio
        by_variance: How much the y component of lattice vector b can deviate from
            the cubic root of ucv in ratio
        cz_variance: How much the z component of lattice vector c can deviate from
            the cubic root of ucv in ratio
        com_dist: If not None, enables COM distance check; sets the COM distance
            limit in Angstrom
        atom_dist: If not None, enables interatomic distance check;
            sets the interatomic distance limit in Angstrom
        specific_radius_proportion: If not None, enables specific radius check;
            sets the sr value of the check
        specific_radius_std: If not None, enables fuzzy sr using half-normal
            distribution (left) to biase towards higher sr
        specific_radius_lower_bound: If not None, enforces a lower bound on
            sr value
        custom_radii: A dict mapping atom species to custom specific radii
        custom_radii_by_pairs: A dict mapping pairs of atom species to custom
            specific radii
        attempt_limit: Number of max attempts per structure generation try
        attempts_per_space_group_change: Number of attempts before the space group
            is refreshed
        attempts_per_unit_cell_change: Number of attempts before the 
            unit cell is refreshed
        attempts_per_wyckoff_list_change: Number of attempts before the wyckoff
            list is refreshed
        attempts_per_closeness_criteria_change: Number of attempts before the
            the closeness criteria is refreshed
        fill_molecule_attempt_limit: Number of attempts to fill molecule into
            an existing structure
        '''
        self._ucv_target = ucv_target
        self._nmpc = nmpc
        self._is_chiral = is_chiral
        self._is_racemic = is_racemic
        self._molecule = molecule
        self._molecule_path = molecule_path
        self._molecule_name = molecule_name
        self._enantiomer_name = enantiomer_name
        self._space_groups_allowed = space_groups_allowed
        self._wyckoff_list = wyckoff_list
        self._ucv_ratio_range = ucv_ratio_range
        self._ucv_std = ucv_std
        self._p_tolerance = p_tolerance
        self._angle_range = angle_range
        self._ax_variance = ax_variance
        self._by_variance = by_variance
        self._cz_variance = cz_variance
        self._com_dist = com_dist
        self._atom_dist = atom_dist
        self._specific_radius_proportion = specific_radius_proportion
        self._specific_radius_std = specific_radius_std
        self._specific_radius_lower_bound = specific_radius_lower_bound
        self._custom_radii = custom_radii
        self._custom_radii_by_pairs = custom_radii_by_pairs
        self._attempt_limit = attempt_limit
        self._attempts_per_space_group_change = \
                attempts_per_space_group_change
        self._attempts_per_unit_cell_change = \
                attempts_per_unit_cell_change
        self._attempts_per_wyckoff_list_change = \
                attempts_per_wyckoff_list_change
        self._attempts_per_closeness_criteria_change = \
                attempts_per_closeness_criteria_change
        self._output_dir = output_dir
        self._output_format = output_format
        self._struct_id_scheme = struct_id_scheme
        self._struct_index_length = struct_index_length
        self._index_tracking_file = index_tracking_file
        self._attempt_tracking_file = attempt_tracking_file

        self._fill_molecule_attempt_limit = fill_molecule_attempt_limit

        if type(self._space_groups_allowed) is int:
            self._space_groups_allowed = [self._space_groups_allowed]

        if self._molecule is None and self._molecule_path is None:
            raise ValueError("Must specify one of molecule or molecule_path")
        elif self._molecule is None:
            self._molecule = Structure()
            self._molecule.loads_try_both(self._molecule_path)
        self._molecule = structure_handling.molecule_COM_move(
                self._molecule)
        self._napm = self._molecule.get_n_atoms()

        if (output_format != "json" and output_format != "geometry"
                and output_format != "both"):
            raise ValueError("Unsupported output format: " + output_format +
                    "; select one from json, geometry or both")
        if not self._output_dir is None:
            misc.safe_make_dir(self._output_dir)

        self._space_group_manager = \
            sgroup.SpaceGroupManager(
                nmpc, is_chiral, is_racemic,
                wyckoff_list=wyckoff_list,
                space_groups_allowed=space_groups_allowed)

        self._unit_cell_generator = \
            UnitCellGenerator(
                angle_range=angle_range, ax_variance=ax_variance,
                by_variance=by_variance, cz_variance=cz_variance,
                p_tolerance=p_tolerance)
        
        # Initializes generation settings
        
        number_of_structures = 2500
        processes_limit = 56

        self._update_space_group()
        self._update_closeness_criteria()
        self._current_space_group_attempts = 0
        self._current_unit_cell_attempts = 0
        self._current_wyckoff_list_attempts = 0
        self._current_closeness_criteria_attempts = 0
        self._generation_success = False
        self._structure = None
        #self._structure_index = total_structures / total_processes
        self._structure_index = number_of_structures / processes_limit 

    def generate_structure_until_index_or_attempts_zero(self):
        # First condition makes sure master attempts have not run out
        # Second condition makes sure that structure index has not hit zero
        while self._check_attempts() and self.generate_structure() != True:
            pass

    def generate_structure(self):
        for i in range(self._attempt_limit):
            if self._generate_structure_single_attempt():
                self._update_structure_index()
                self._fill_generated_structure_info()
                self._print_generated_structure_info()
                if self._structure_index < 0:
                    print_time_log("Generation success, but structure index "
                            "found to be < 0; newly generated structure will "
                            "not be outputed; this is expected if generation "
                            "has just completed")
                    return True
            
                print_time_log("Generation success! Outputing structure to "
                        + self._output_dir)
                output_structure(
                        self._structure, self._output_dir,
                        self._output_format)
                self._update_space_group()
                return self._structure

            self._update_generation_counters()

        print_time_log("Generation failed")
        return False

    def get_latest_generated_structure(self):
        return self._structure

    def _generate_structure_single_attempt(self):
        self._generation_success = False
        self._structure = Structure()
        self._current_wyckoff_position = 0
        self._unit_cell_generator.\
            set_latest_generated_unit_cell_for_structure(self._structure)
        return self._fill_molecule_by_wyckoff_list()

    def _fill_molecule_by_wyckoff_list(self):
        original_struct = copy.deepcopy(self._structure)
        for i in range (self._fill_molecule_attempt_limit):
            struct_new = place_molecule_space_group(
                    self._structure,
                    self._molecule,
                    self._space_group.space_group_number,
                    self._wyckoff_list_used[self._current_wyckoff_position],
                    create_duplicate=True)

            if closeness_check(
                    struct_new, self._napm, self._com_dist, self._atom_dist,
                    self._sr, self._custom_radii, self._custom_radii_by_pairs):
                self._structure = struct_new
                if self._current_wyckoff_position == \
                        len(self._wyckoff_list_used) - 1:
                    self._generation_success = True
                    return True
                else: #Move on to the next
                    self._current_wyckoff_position += 1
                    if self._fill_moelcule_by_wyckoff_list():
                        # Recurisvely exit successful generation
                        return True
                    self._current_wyckoff_position -= 1

        return False #Unsuccessful fill
    
    def _print_generated_structure_info(self):
        print_time_log("Structure generation success with space group %i, "
                "wyckoff position list %s, unit cell volume %f, and "
                "sr minimum: %s" 
                % (self._space_group.space_group_number,
                   str(self._wyckoff_list_used),
                   self._ucv, str(self._sr)))

    def _fill_generated_structure_info(self):
        misc.struct_id_assign(self._structure, self._struct_id_scheme,
                str(self._structure_index).zfill(self._struct_index_length))
        self._structure.set_property(
                "space_group", self._space_group.space_group_number)
        self._structure.set_property(
                "wyckoff_list", self._wyckoff_list_used)
        self._structure.set_property("NMPC", self._nmpc)
        self._structure.set_property(
                "molecule_name", self._molecule_name)
        self._structure.set_property(
                "enantiomer_name", self._enantiomer_name)

    def _update_generation_counters(self):
        self._current_space_group_attempts += 1
        self._current_unit_cell_attempts += 1
        self._current_wyckoff_list_attempts += 1
        self._current_closeness_criteria_attempts += 1
        if (self._current_space_group_attempts >=
                self._attempts_per_space_group_change):
            self._update_space_group()
            self._current_space_group_attempts = 0
            self._current_unit_cell_attempts = 0
            self._current_wyckoff_list_attempts = 0
        if (self._current_unit_cell_attempts >=
                self._attempts_per_unit_cell_change):
            self._update_unit_cell()
            self._current_unit_cell_attempts = 0
        if (self._current_wyckoff_list_attempts >=
                self._attempts_per_wyckoff_list_change):
            self._update_wyckoff_list()
            self._current_wyckoff_list_attempts = 0
        if (self._current_closeness_criteria_attempts >=
                self._attempts_per_closeness_criteria_change):
            self._update_closeness_criteria()
            self._current_closeness_criteria_attempts = 0

    def _update_space_group(self):
        # This method updates all the temporary generation settings
        # as space group change requires unit cell and wyckoff list change
        #self._space_group = \
        #    self._space_group_manager.get_space_group_randomly()
        #self._space_group = \
        #     self._space_group_manager.get_space_group_uniformly(self._output_dir)
        self._space_group = \
             self._space_group_manager.get_space_group_uniformly_manny(self._output_dir)
        self._bravais_system = self._space_group.get_bravais_system_type()
        self._update_wyckoff_list()
        self._update_unit_cell()
        print_time_log("Updated space group to: " +
                str(self._space_group.space_group_number))

    def _update_wyckoff_list(self):
        self._wyckoff_list_used = \
            self._space_group_manager.get_wyckoff_list_randomly()

    def _update_unit_cell(self):
        self._update_ucv()
        self._unit_cell_generator.set_ucv(self._ucv)
        self._unit_cell_generator.set_bravais_system(
            self._bravais_system)
        self._unit_cell_generator.generate_unit_cell()

    def _update_ucv(self):
        if self._ucv_std is None:
            # Linear variation mode
            self._ucv = (self._ucv_target *
                random.uniform(self._ucv_ratio_range[0],
                               self._ucv_ratio_range[1]))
        else:
            # Half-normal distribution mode
            self._ucv = \
                half_gaussian_sampling_upper(
                    self._ucv_target * self._ucv_ratio_range[0],
                    self._ucv_std,
                    self._ucv_target * self._ucv_ratio_range[1])
        print_time_log("Updated unit cell volume to: " + str(self._ucv))


    def _update_closeness_criteria(self):
        if self._specific_radius_proportion is None:
            self._sr = None
        elif self._specific_radius_std is None:
            self._sr = self._specific_radius_proportion
        else:
            self._sr = \
                    half_gaussian_sampling_lower(
                            self._specific_radius_proportion,
                            self._specific_radius_std,
                            self._specific_radius_lower_bound)

    def _update_structure_index(self):
        if self._index_tracking_file is None:
            self._structure_index -= 1
        else:
            self._structure_index = retrieve_integer_and_subtract(
                    self._index_tracking_file, subtract=1)
        print_time_log("Structure index updated to %i" % self._structure_index)

    def _check_attempts(self):
        return (self._attempt_tracking_file is None or
                retrieve_integer_and_subtract(
                    self._attempt_tracking_file, subtract=1) > 0)
        

class UnitCellGenerator(object):
    '''
    This module handles the generation of unit cell
    '''
    def __init__(self,
        ucv=None, bravais_system=None, angle_range=[30,150],
        ax_variance=[0.5,2], by_variance=[0.5,2], cz_variance=[0.5,2],
        p_tolerance=0.25):
        '''
        bravais_system: the type of the bravais system to use,
        1-->Triclinic, 2-->monoclinic, 3-->orthorhombic,
        4-->tetragonal, 5-->Cubic
        '''
        if ucv != None:
            self._set_ucv(ucv)
        if bravais_system != None:
            self.set_bravais_system(bravais_system)
        self._ucv = ucv
        self._bravais_system = bravais_system
        self._angle_range = angle_range
        self._ax_variance = ax_variance
        self._by_variance = by_variance
        self._cz_variance = cz_variance
        self._principal_variance = \
            [ax_variance, by_variance, cz_variance]
        self._p_tolerance = p_tolerance
        self._generation_success = False

    def set_ucv(self, ucv):
        self._ucv = ucv
        self._ucv_cube_root = ucv ** (1.0/3)
        self._generation_success = False

    def set_bravais_system(self, bravais_system):
        self._bravais_system = bravais_system
        self._generation_success = False

    def generate_unit_cell(self):
        '''
        Generates and returns a tuple of:
            (lattice_vector_a, lattice_vector_b, lattice_vector_c,
             a, b, c, alpha, beta, gamma, bravais_system)
        '''
        if self._ucv is None:
            raise ValueError("Unit cell volume not initialized;"
                             "call set_ucv before generation attempt")

        if self._bravais_system is None:
            self._bravais_system_used = int(random.uniform(1,6))
        else:
            self._bravais_system_used = self._bravais_system

        for i in range(1000):
            # Generation has a certain rate of failure
            # Uses a for-loop to avoid dead loop
            if self._bravais_system_used == 1:
                result = self._generate_triclinic_unit_cell()
            if self._bravais_system_used == 2:
                result = self._generate_monoclinic_unit_cell()
            if self._bravais_system_used == 3:
                result = self._generate_orthorhombic_unit_cell()
            if self._bravais_system_used == 4:
                result = self._generate_tetragonal_unit_cell()
            if self._bravais_system_used == 5:
                result = self._generate_cubic_unit_cell()

            if result:
                self._generation_success = True
                self._finalize_lattice_vectors()
                print_time_log("Lattice vector generation success; a = %s, "
                        "b = %s, c = %s"
                        % (str(self._lattice_vector_a),
                           str(self._lattice_vector_b),
                           str(self._lattice_vector_c)))
                return self.get_last_generated_unit_cell()

        raise ValueError("Repeated failure to generate lattice vector")

    def set_latest_generated_unit_cell_for_structure(self, struct):
        struct.set_property("lattice_vector_a", self._lattice_vector_a)
        struct.set_property("lattice_vector_b", self._lattice_vector_b)
        struct.set_property("lattice_vector_c", self._lattice_vector_c)

        struct.set_property("a", self._a)
        struct.set_property("b", self._b)
        struct.set_property("c", self._c)

        struct.set_property("alpha", self._alpha)
        struct.set_property("beta", self._beta)
        struct.set_property("gamma", self._gamma)
        
        struct.set_property("bravais_system", self._bravais_system_used)
        struct.set_property("unit_cell_volume", self._ucv)
        struct.set_property("cell_vol", self._ucv)


    def get_last_generated_unit_cell(self):
        return (self._lattice_vector_a, self._lattice_vector_b,
            self._lattice_vector_c, self._a, self._b, self._c,
            self._alpha, self._beta, self._gamma,
            self._bravais_system_used)

    def _generate_triclinic_unit_cell(self):
        if not self._set_triclinic_angles():
            return False
        if not self._set_principal_components_abc():
            return False
        return True

    def _generate_monoclinic_unit_cell(self):
        if not self._set_monoclinic_angles():
            return False
        if not self._set_principal_components_abc():
            return False
        return True

    def _generate_orthorhombic_unit_cell(self):
        self._set_orthogonal_angles()
        if not self._set_principal_components_abc():
            return False
        return True

    def _generate_tetragonal_unit_cell(self):
        self._set_orthogonal_angles()
        self._set_principal_components_aac()
        return True

    def _generate_cubic_unit_cell(self):
        self._set_orthogonal_angles()
        self._set_principal_components_aaa()
        return True

    def _set_triclinic_angles(self):
        self._alpha, self._beta, self._gamma = \
            [random.uniform(self._angle_range[0],
                            self._angle_range[1])
                for x in range(3)]
        return self._check_valid_angles()

    def _set_monoclinic_angles(self):
        self._alpha = self._gamma = 90
        self._beta = random.uniform(self._angle_range[0],
                                    self._angle_range[1])
        return self._check_valid_angles()

    def _set_orthogonal_angles(self):
        self._alpha = self._beta = self._gamma = 90

    def _check_valid_angles(self):
        # Makes sure that generated angles do not result in a
        # unit cell that is (too) flat
        if self._alpha + self._beta + self._gamma > 355:
            return False
        if (self._alpha > self._beta + self._gamma - 5 or
            self._beta > self._gamma + self._alpha - 5 or
            self._gamma > self._alpha + self._beta - 5):
            return False
        if calc_p(self._alpha, self._beta, self._gamma) < self._p_tolerance:
            return False
        return True

    def _set_principal_components_abc(self):
        # Generates three different principal components

        # In order for equivalence between the three lattice vectors,
        # a number is randomized, which is then interpreted by this value into
        # actual sequence in which the lattice vectors are generated
        sequence_interp= \
            {0:[0,1,2], 1:[0,2,1], 2:[1,0,2],
             3:[1,2,0], 4:[2,0,1], 5:[2,1,0]}
        sequence = int(random.uniform(0,6))
        v1, v2, v3 = sequence_interp[sequence]

        result = [0, 0, 0]
        result[v1] = self._ucv_cube_root * \
            random.uniform(
                self._principal_variance[v1][0],
                self._principal_variance[v1][1])

        # v2 has to leave a certain range for v3 to fall into the required range
        vol_remain = self._ucv / result[v1]
        lower = max(self._principal_variance[v2][0] * self._ucv_cube_root,
                    vol_remain / (self._ucv_cube_root * 
                                  self._principal_variance[v3][1]))
        upper = min(self._principal_variance[v2][1] * self._ucv_cube_root,
                    vol_remain / (self._ucv_cube_root * 
                                  self._principal_variance[v3][0]))
        if lower > upper:
            return False

        result[v2] = random.uniform(lower, upper)
        result[v3] = vol_remain / result[v2]

        self._ax = result[0]
        self._by = result[1]
        self._cz = result[2]
        return True

    def _set_principal_components_aac(self):
        # Generates principal components with a=a!=c
        lower = max(
            self._ucv_cube_root * self._ax_variance[0],
            self._ucv_cube_root * self._by_variance[0],
            (self._ucv / (self._ucv_cube_root * self._cz_variance[1]))**0.5)
        
        upper = min(
            self._ucv_cube_root * self._ax_variance[1],
            self._ucv_cube_root * self._by_variance[1],
            (self._ucv / (self._ucv_cube_root * self._cz_variance[0]))**0.5)

        if lower > upper:
            # raise because repeated try will not fix
            raise ValueError(
                "Principal component variance unsuitable for tetragonal "
                "structure generation")
        
        self._ax = self._by = random.uniform(lower, upper)
        self._cz = self._ucv / self._ax**2

    def _set_principal_components_aaa(self):
        # Generate three principal components as the same
        self._ax = self._by = self._cz = self._ucv_cube_root

    def _finalize_lattice_vectors(self):
        alpha = np.deg2rad(self._alpha)
        beta = np.deg2rad(self._beta)
        gamma = np.deg2rad(self._gamma)
        self._lattice_vector_a = [0, 0, 0]
        self._lattice_vector_b = [0, 0, 0]
        self._lattice_vector_c = [0, 0, 0]

        self._a = self._ax
        self._lattice_vector_a[0] = self._ax
        self._b = self._by / np.sin(gamma)
        self._lattice_vector_b[0] = self._b * np.cos(gamma)
        self._lattice_vector_b[1] = self._by
        self._c = \
            (self._cz * self._by / (
                np.sin(beta)**2 * self._by**2 - self._b *
                np.cos(alpha) - self._lattice_vector_b[0] * np.cos(beta)) ** 2
            ) ** 0.5
        self._lattice_vector_c[0] = self._c * np.cos(beta)
        self._lattice_vector_c[1] = (
            (self._c * self._b * np.cos(alpha) -
             self._lattice_vector_c[0] * self._lattice_vector_b[0]) / self._by)
        self._lattice_vector_c[2] = self._cz
                
def place_molecule_space_group(struct, molecule, sgn, wycn, create_duplicate=True):
    '''
    Randomly places a molecule into the cell according to the symmetry requirement of the Wyckoff position
    Both structure and molecule will should be of structure.Structure() class
    '''
    if create_duplicate:
        struct = copy.deepcopy(struct)
    sg = sgroup.Sgroup(sgn)
    lattice_vectors = np.transpose([struct.properties["lattice_vector_a"],
                       struct.properties["lattice_vector_b"],
                       struct.properties["lattice_vector_c"]])
    latinv = np.linalg.inv(lattice_vectors)

    # Allow 10 times of generation to avoid accidental
    # overlapping of molecules leading to wrong special
    # Wyckoff position determination
    for i in range (10):
        # Gets the fractional coordinates the fits the corresponding
        # Wyckoff number 
        fst_frac = sg.wycgen(wycn) 

        fst_orien = random_rotation_matrix()
    
        list_of_fracs = [] #List of translation vectors
        list_of_orien = [] #List of potential rotational matrices
    
        # Loops over Bravais centering
        for ll in range (0, len(sg.btran)):
            # Loops over point group operations
            for l in range (0, len(sg.op)): 
                new_frac = np.dot(sg.op[l],fst_frac)
                #Add the Translation required by the operation
                new_frac = np.add(new_frac,sg.trans[l])
                #Add the translation from the Bravais centering
                new_frac = np.add(new_frac,sg.btran[ll]) 
                frac_coor_adjust(new_frac, False)

                # Fractional rotation to absolute rotation
                space_group_rotation = \
                      np.dot(lattice_vectors,
                             np.dot(sg.op[l], latinv))
                new_orien = np.dot(space_group_rotation,fst_orien)
            
                c = True
                for i in range (len(list_of_fracs)): 
                # If a special Wyckoff selection, 
                # after the symmetry operation is applied,
                # Two or more molecules will sit on the same site.
                # Here groups them into sites for later selection
                    diff = np.linalg.norm(
                          np.subtract(list_of_fracs[i], new_frac))

                    if diff < 0.00001: #Same site
                        c = False
                        list_of_orien[i].append(new_orien) 
                        break
                if c:
                    list_of_fracs.append(new_frac)
                    list_of_orien.append([new_orien])
        # Breaks if the correct amount of unique molecule positions
        # are generated.
        if len(list_of_fracs) == sg.wmult[wycn]*sg.bmult:
            break

    if len(list_of_fracs) != sg.wmult[wycn]*sg.bmult:
        message = ("Repeated failure to match number of molecules "+
        "to that of Wyckoff position; check space group module for: "
        "space group=%i, wycn=%i" % (sg.space_group_number,wycn))
        raise RuntimeError(message) #!!!

    final_orien = []
    for candidates in list_of_orien: 
    # Picks one of the several sitting on the same site 
    # if a special Wyckoff position is selected
        chosen = int(random.uniform(0,len(candidates)))
        final_orien.append(candidates[chosen])

    final_trans = []
    for frac in list_of_fracs:
        # Converts to absolute coordinates
        final_trans.append(np.dot(lattice_vectors,frac))

    for i in range (len(final_trans)):
        # Now add the molecules into the structure
        k = "molecule_harris_info"
        n = misc.molecule_harris_info_prep(final_trans[i], 
                           final_orien[i])
        if k in struct.properties:
            struct.properties[k].append(n)
        else:
            struct.properties[k] = [n]

        # Get geometry of the new molecule
        new_mole = structure_handling.cell_transform_mat(
              molecule, final_orien[i])
        new_mole = structure_handling.cell_translation(
              new_mole, final_trans[i], False)

        struct.geometry = np.concatenate((struct.geometry,
                             new_mole.geometry))
    
    return struct
    

def closeness_check(
        struct, napm, com_dist=None, atom_dist=None, sr_proportion=None,
        custom_radii={}, custom_radii_by_pairs={}):
    '''
    Calls the respective functions in structure_handling to conduct the
    closeness check on the structure. If the closeness criteria is None,
    the respective closeness check is disabled
    '''
    original_struct = struct
    struct = copy.deepcopy(struct)
    if struct.get_n_atoms() % napm != 0:
        raise ValueError("napm does not divide the number of atoms in struct.")

    struct = structure_handling.cell_modification(
            struct, len(struct.geometry)/napm, napm, create_duplicate=False)
    nmpc = struct.get_n_atoms() / napm

    if not atom_dist is None or not sr_proportion is None:
        lowerbound = (atom_dist, sr_proportion)
        continue_check = structure_handling.full_nearest_image(struct, nmpc=nmpc, lowerbound=lowerbound) 

        if not continue_check:
            return False

        min_atom_dist, min_sr = structure_handling.calc_atom_dist_manny(struct,
                                                nmpc=nmpc, 
                                                lowerbound=lowerbound)
        result = True
        if atom_dist != None:
            if min_atom_dist < atom_dist:
                result = False
        if sr_proportion != None:
            if min_sr < sr_proportion:
                result = False
        if result == True:
            original_struct.properties['sr'] = [min_sr]

    if not result:
        return False
    return True

def calc_p (alpha,beta,gamma):
    '''
    Calculates the adjustment to the cell volume posed by the angles
    volume=p*a*b*c
    '''
    ca=np.cos(np.deg2rad(alpha))
    cb=np.cos(np.deg2rad(beta))
    cc=np.cos(np.deg2rad(gamma))
    return (1-ca*ca-cb*cb-cc*cc+2*ca*cb*cc)**0.5
