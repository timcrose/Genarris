import os, copy
from core.structure import Structure
from utilities.misc import list_subdirectories, list_directory_file, \
        safe_make_dir, output_structure, safe_rmdir, safe_remove_files
from utilities.write_log import print_time_log, print_time_warning, \
        print_time_error
from copy import deepcopy

def get_fhi_aims_extractor(inst, section, omit_output=False):
    aims_output_file = inst.get_with_default(
            section, "aims_output_file", "aims.out")
    original_structure_file = inst.get_with_default(
            section, "original_structure_file", "original_structure.json")
    original_structure_dir = inst.get_or_none(
            section, "original_structure_dir")
    extract_energy = inst.get_boolean(section, "extract_energy")
    energy_property_name = inst.get_with_default(
            section, "energy_property_name", "energy")
    update_geometry = inst.get_boolean(section, "update_geometry")
    check_execution_complete = inst.get_boolean(
            section, "check_execution_complete")
    strict_success = inst.get_boolean(section, "strict_success")
    clean_up = inst.get_boolean(section, "clean_up")
    if not omit_output:
        output_dir = inst.get_or_none(section, "output_dir")
        output_format = inst.get_or_none(section, "output_format")
    else:
        output_dir = None
        output_format = None

    return FhiAimsExtractor(
            aims_output_file=aims_output_file,
            original_structure_file=original_structure_file,
            original_structure_dir=original_structure_dir,
            extract_energy=extract_energy,
            energy_property_name=energy_property_name,
            update_geometry=update_geometry,
            check_execution_complete=check_execution_complete,
            strict_success=strict_success,
            clean_up=clean_up,
            output_dir=output_dir,
            output_format=output_format)

class FhiAimsExtractor(object):
    def __init__(self, aims_output_file="aims.out",
            original_structure_file="original_structure.json",
            original_structure_dir=None,
            extract_energy=False, energy_property_name="energy",
            update_geometry=False,
            check_execution_complete=False,
            strict_success=False,
            output_dir=None, output_format="json",
            clean_up=False):
        '''
        aims_output_file: File to look for in each execution folder as 
            FHI-aims output
        original_structure_file: File to look for in each execution folder as 
            original structure in JSON format.
        original_structure_dir: Directory where additional original structures 
            are stored outside of calculation folder. Structures are matched 
            according to struct_id and calc folder's basename.
        extract_energy: Whether or not to extract energy from aims output
        energy_property_name: Property name to store the structure
        update_geometry: Whether or not to update the geometry from aims output
        check_execution_complete: Whether or not to look for
            _exeuction_complete file that indicates execution complete
        strict_success: Whether or not to look for "Have a nice day" to 
            determine execution success (vs. "Leaving FHI-aims")
        clean_up: Whether or not to clean up successfully extracted folder;
            files to be cleaned up include:
            control.in, geometry.in, geometry.in.next_step,
            aims output, oiginal structure file, _execution_complete
            If the folder is empty after removing these files, removes
            the folder as well
        '''
        self._aims_output_file = aims_output_file
        self._extract_energy = extract_energy
        self._energy_property_name = energy_property_name
        self._update_geometry = update_geometry
        self._check_execution_complete = check_execution_complete
        self._strict_success = strict_success
        self._output_dir = output_dir
        self._output_format = output_format
        self._clean_up = clean_up

        self._original_structure_file = original_structure_file 
        self._original_structure_dir = os.path.abspath(original_structure_dir)
        if not self._original_structure_dir is None:
            self._load_original_structures()
        if not self._output_dir is None:
            safe_make_dir(self._output_dir);

        self._dirs_to_extract = []

    def extract(self, calculation_dir):
        self._dirs_to_extract.append(calculation_dir)
        return self._extract();

    def extract_batch(self, calculation_dir):
        subdirs = list_subdirectories(calculation_dir)
        subdirs.sort()
        self._dirs_to_extract += subdirs
        return self._extract()

    def _extract(self):
        result = []
        for calculation_dir in self._dirs_to_extract:
            struct = self._extract_single(calculation_dir)
            if not struct is False:
                result.append(struct)
        return result

    def _extract_single(self, calculation_dir):
        '''
        Main extraction workflow
        '''
        print_time_log("Extracting FHI aims output from folder %s" 
                % calculation_dir)

        if not self._check_complete(calculation_dir):
            return False

        struct = self._find_original_structure(calculation_dir)
        if struct is False:
            struct = Structure()
            found_geometry = False
        else:
            found_geometry = True
        
        geo_file = os.path.join(calculation_dir, "geometry.in")
        if not found_geometry:
            found_geometry = self._load_geometry_file(struct,
                    os.path.join(calculation_dir, "geometry.in"))
            if found_geometry:
                print_time_log("Loaded original geometry file from "
                        "geometry.in in dir: " + calculation_dir)

        extracted_something = False
        if self._update_geometry:
            success = self._extract_geometry(struct, calculation_dir)
            if success:
                extracted_something = True
                found_geometry = True

        if self._extract_energy:
            success = self._extract_energy_do(struct, calculation_dir)
            if success:
                extracted_something = True

        if not found_geometry:
            print_time_error("Extraction failure as no geometry is found. "
                    "Directory: " + calculation_dir)
            return False

        if not extracted_something:
            print_time_error("Extraction failure as nothing is extracted. "
                    "Directory: " + calculation_dir)
            return False

        print_time_log("Successfully extracted dir: " + calculation_dir)
        if struct.struct_id is None:
            struct.struct_id = os.path.basename(calculation_dir)
        if not self._output_dir is None:
            output_structure(struct, self._output_dir, self._output_format)
        
        self._clean_up_directory(calculation_dir)
        return struct

    def _find_original_structure(self, calculation_dir):
        struct = Structure()
        if not self._original_structure_file is None:
            path = os.path.join(calculation_dir, self._original_structure_file)
            if self._load_original_structure(struct, path):
                print_time_log(
                        "Extracted in-folder original structure file: " + path)
                return struct

        matched_struct = [x for x in self._original_structures
                if x[0].struct_id == os.path.basename(calculation_dir)]
        if len(matched_struct) > 2:
            raise ValueError("Multiple original structure files matched "
                "with subdirectory %s. Please check that structures in "
                "%s have unique struct_id"
                % (calculation_dir, self._original_structure_dir))
        elif len(matched_struct) == 1:
            print_time_log("Matched calculation dir %s with original structure"
                    " file %s" % (calculation_dir, matched_struct[0][1]))
            return copy.deepcopy(matched_struct[0][0])

        print_time_warning("Continuing extraction without original "
                "structure file in JSON")
        return False

    def _load_original_structures(self):
        if not os.path.isdir(self._original_structure_dir):
            raise ValueError("Original structure directory not found: "
                    + self._original_structure_dir)

        print_time_log("Loading original structures from dir: "
                + self._original_structure_dir)

        files = list_directory_file(self._original_structure_dir)
        self._original_structures = []
        for f in files:
            struct = Structure()
            if self._load_original_structure(struct, f):
                self._original_structures.append((struct, f))

    def _load_original_structure(self, struct, original_structure_file):
        if os.path.isfile(original_structure_file):
            try:
                struct.build_geo_from_json_file(original_structure_file)
                return True
            except:
                print_time_warning("Failed to load as JSON "
                        "original structure file: "
                        + original_structure_file)
        else:
            print_time_warning("Failed to find original structure "
                    "file " + original_structure_file)
        return False

    def _load_geometry_file(self, struct, geometry_file):
        if os.path.isfile(geometry_file):
            try:
                struct.build_geo_from_atom_file(geometry_file)
                return True
            except:
                print_time_error("Failed to load geometry "
                        "from file: " + geometry_file)
        else:
            print_time_warning("No geometry file: " + geometry_file)


    def _check_complete(self, calculation_dir):
        if self._check_execution_complete and \
                not os.path.isfile(
                        os.path.join(calculation_dir, "_execution_complete")):
            print_time_warning("_execution_complete file not found; skipping "
                    "extraction of folder: " + calculation_dir)
            return False

        aims_file = os.path.join(calculation_dir, self._aims_output_file)
        if not os.path.isfile(aims_file):
            print_time_error("Extraction failed as aims output file "
                    "does not exist: " + aims_file)
            return False

        if not self._check_success(aims_file):
            print_time_error("Extraction failed as aims execution "
                    "determined as failure: %s. Strict success set to %s."
                    % (aims_file, str(self._strict_success)))
            return False
        return True


    def _check_success(self, aims_file):
        try:
            aims_out = open(aims_file,"r")
        except:
            return False #File missing

        counter = 0
        while True:
            line = aims_out.readline()
            if (not self._strict_success and "Leaving FHI-aims" in line)\
                or "Have a nice day" in line:
                    aims_out.close()
                    return True
            elif line == '': counter += 1
            else: counter = 0

            if counter > 10:
            #Allowing 10 empty lines in a row before determining eof
                break

        aims_out.close()
        return False

    def _extract_geometry(self, struct, calculation_dir):
        next_step_file = os.path.join(
                calculation_dir, "geometry.in.next_step")
        if not os.path.isfile(next_step_file):
            print_time_error("Failed to find geometry.in.next_step "
                    "in dir: " + calculation_dir)
            return False
        test_struct = Structure()
        if not self._load_geometry_file(test_struct, next_step_file):
            print_time_error("Failed to extract geometry.in.next_step "
                    "from dir: " + calculation_dir)
            return False
        struct.clear_geometry()
        self._load_geometry_file(struct, next_step_file)
        print_time_log("Updated geometry from geometry.in.next_step from dir: "
                + calculation_dir)
        return True

    def _extract_energy_do(self, struct, calculation_dir):
        # By this point, we know that calculation_dir/aims_output_file exists
        aims_out = open(os.path.join(
            calculation_dir, self._aims_output_file))

        s = '  | Total energy of the DFT / Hartree-Fock s.c.f. calculation      :'
        while True:
            line = aims_out.readline()
            if not line:
                print_time_error("Failed to extract converged energy from file: "
                    + os.path.join(calculation_dir, self._aims_output_file))
                return False # energy not converged
            if s in line:
                tokens = line.split()
                energy = float(tokens[11])
                struct.set_property(self._energy_property_name, energy)
                aims_out.close()
                return True

        return False

    def _clean_up_directory(self, calculation_dir):
        if not self._clean_up: return

        files = ["control.in", "geometry.in", "geometry.in.next_step",
                "_execution_complete", self._aims_output_file]
        if not self._original_structure_file is None:
            files.append(self._original_structure_file)
        safe_remove_files([os.path.join(calculation_dir, x)
                           for x in files])
        if len(os.listdir(calculation_dir)) == 0:
            safe_rmdir(calculation_dir)

