from utilities import misc, write_log
from core.structure import Structure
import os, multiprocessing


class BatchSingleStructureOperation(object):
    def __init__(self, operation, name="SingleStructureOperation",
            args=(), kwargs={},
            structure_list=None, structure_path=None,
            structure_dir=None, structure_dir_depth=0, structure_suffix="",
            output_dir=None, output_format="json",
            processes_limit=1, disable_preload=False):
        '''
        operation: Function or list of functions to perform on a structure
        args: A tuple or a list of tuple of arguments passed to the function
            Must match length of operation if a list

        structure_list: List of Structure() class object
        structure_path: A single path to a structure
        structure_dir: A path to a directory of structures
        output_dir: A path where the completed structure will be printed with
            file name based on the struct_id
        output_format: Supports "json", "geometry", or "both"
        disable_preload: Disable preloading the structure list into this 
            object and instead use .fin file to mark completion. Useful for
            parallelization realized through external launch of multiple 
            instances of Genarris. Cannot be set to True when processes_limit
            > 1, which indicates that internal multiprocessing parallel
            scheme is used.
        '''
        self._operation = operation
        self._name = name
        self._args = args
        self._kwargs = kwargs
        self._structure_list = structure_list
        self._structure_path = structure_path
        self._structure_dir = structure_dir
        self._structure_dir_depth = structure_dir_depth
        self._structure_suffix = structure_suffix
        self._output_dir = output_dir
        self._output_format = output_format
        self._processes_limit = processes_limit
        self._disable_preload = disable_preload

        self._load_structures()
        if processes_limit > 1 and disable_preload:
            raise ValueError("Cannot disable preload when processes_limit "
                    "is greater than 1.")
        if (output_format != "json" and output_format != "geometry"
                and output_format != "both"):
            raise ValueError("Unsupported output format: " + output_format +
                    "; select one from json, geometry or both")
        misc.safe_make_dir(self._output_dir)

    def run(self):
        if not self._disable_preload:
            write_log.print_time_log(
                    "Executing single structure operation %s on %i structures "
                    "with %i processes" % (self._name,
                        len(self._structure_list), self._processes_limit))
        if self._processes_limit > 1:
            result = self._parallel_run()
        else:
            result = self._serial_run()

        if not self._disable_preload:
            write_log.print_time_log(
                    "Completed single structure operation %s on %i structures "
                    "with %i processes" % (self._name,
                        len(self._structure_list), self._processes_limit))
        return result

    def _parallel_run(self):
        p = multiprocessing.Pool(self._processes_limit)
        op_args_list = \
                [(self._operation, struct, self._args, self._kwargs,
                    self._output_dir, self._output_format)
                    for struct in self._structure_list]
        
        if type(self._operation) is list:
            return p.map(_run_operations_with_arguments, op_args_list)
        else:
            return p.map(_run_operation_with_arguments, op_args_list)

    def _serial_run(self):
        result = []
        if not self._disable_preload:
            for struct in self._structure_list:
                op_args = (self._operation, struct, self._args, self._kwargs,
                        self._output_dir, self._output_format)
                if type(self._operation) is list:
                    result.append(_run_operations_with_arguments(op_args))
                else:
                    result.append(_run_operation_with_arguments(op_args))
        else:
            # TODO: Create serial run for non-preload
            pass

        return result

    def _load_structures(self):
        set_options = 0
        if self._structure_list != None:
            set_options += 1

        if self._structure_path != None:
            set_options += 1

        if self._structure_dir != None:
            set_options += 1

        if set_options != 1:
            raise ValueError("Must specificy exactly one of structure_list, "
                    "structure_path or structure_dir; %i is/are set"
                    % set_options)

        if self._structure_list != None:
            # Structures already loaded
            return
        self._structure_list = []
        if self._structure_path != None:
            self._load_structure(self._structure_path)
        if self._structure_dir != None and not self._disable_preload:
            self._preload_structure_list()

    def _load_structure(self, path):
        struct = Structure()
        try:
            struct.build_geo_from_json_file(path)
        except:
            try:
                struct.build_geo_from_atom_file(path)
                struct_id = os.path.basename(path)[:-self._structure_suffix]
                struct.struct_id = struct_id
            except:
                write_log.print_time_log(
                        "WARNING: Failed to load structure file: "
                        + path)
                return
                
        self._structure_list.append(struct)

    def _preload_structure_list(self):
        # TODO: Merge with general load structure directory util
        for path in misc.list_directory_file(
                self._structure_dir,
                depth=self._structure_dir_depth,
                suffix=self._structure_suffix):
            self._load_structure(path)

    def _add_structure_to_queue(self):
        path = misc.get_and_lock_file(
                self._structure_dir,
                depth=self._structure_dir_depth,
                suffix=self._structure_suffix,
                message="Locked by " + self._name)
        if path != False:
            self._load_structure(path)
            return True
        return False


def load_batch_single_structure_operation_keywords(inst, section):
    kwargs = inst.get_kv_single_section(
            section, [
                "structure_path", "structure_dir",
                "structure_suffix", "output_dir", "output_format"],
            eval=False)
    kwargs.update(inst.get_kv_single_section(section,
        ["structure_dir_depth"], eval=True))
    kwargs["processes_limit"] = inst.get_processes_limit(section)
    return kwargs

def _run_operation_with_arguments(op_args):
    (operation_list, struct, args_list, kwargs_list,
            output_dir, output_format) = op_args
    if len(operation_list) != len(args_list):
        raise ValueError(
                "Length of operation_list and args_list must match.")
    if len(operation_list) != len(kwargs_list):
        raise ValueError(
                "Length of operation_list and kwargs_list must match.")

    returned_result = []
    for operation, args, kwargs in zip(operation_list, args_list, kwargs_list):
        result = _run_operation_with_arguments(
                operation, struct, args, kwargs, None, "")
        if type(result) is tuple:
            _check_return_tuple(result)
            struct = result[0]
            returned_result.append(result[1])
        else:
            _check_return_single(result)
            struct = result

    _output_structure(struct, output_dir, output_format)
    return struct

def _run_operation_with_arguments(op_args):
    operation, struct, args, kwargs, output_dir, output_format = op_args
    result = operation(struct, *args, **kwargs)
    if type(result) is tuple:
        _check_return_tuple(result)
        _output_structure(result[0], output_dir, output_format)
    else:
        _check_return_single(result)
        _output_structure(result, output_dir, output_format)
    return result

def _check_return_single(result):
    if not type(result) is Structure:
        raise ValueError("Single structure operation must return a Structure "
                "class object")

def _check_return_tuple(result):
    if len(result) > 2 or not type(result[0]) is Structure:
        raise ValueError("Single structure operation returning tuple "
            "must contain only two elements: the updated "
                "structure and additional result")

def _output_structure(struct, output_dir, output_format):
    if output_dir == None:
        return
    if output_format == "json" or output_format == "both":
        path = os.path.join(output_dir, struct.struct_id + ".json")
        f = open(path, "w")
        f.write(struct.dumps())
        f.close()
    if output_format == "geometry" or output_format == "both":
        path = os.path.join(output_dir, struct.struct_id + ".in")
        f = open(path, "w")
        f.write(struct.get_geometry_atom_format())
        f.close()
