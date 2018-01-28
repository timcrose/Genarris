from utilities import misc, write_log
from core.structure import Structure
import os, multiprocessing, shutil, time
from utilities.misc import output_structure, list_subdirectories, \
       retrieve_float
from write_log import print_time_log


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
        if not self._output_dir is None:
            misc.safe_make_dir(self._output_dir)

    def run(self):
        if not self._disable_preload:
            print_time_log(
                    "Executing single structure operation %s on %i structures "
                    "with %i processes" % (self._name,
                        len(self._structure_list), self._processes_limit))
        if self._processes_limit > 1:
            result = self._parallel_run()
        else:
            result = self._serial_run()

        if not self._disable_preload:
            print_time_log(
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
                print_time_log(
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

    output_structure(struct, output_dir, output_format)
    return struct

def _run_operation_with_arguments(op_args):
    operation, struct, args, kwargs, output_dir, output_format = op_args
    result = operation(struct, *args, **kwargs)
    if type(result) is tuple:
        _check_return_tuple(result)
        output_structure(result[0], output_dir, output_format)
    else:
        _check_return_single(result)
        output_structure(result, output_dir, output_format)
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

class BatchSingleStrucutreExternalEvaluation(object):
    def __init__(self, execute_commands,
            setup_method, restart_method, execute_method, extract_method,
            setup_method_args=(), restart_method_args=(),
            execute_method_args=(), extract_method_args=(),
            structure_list=None, structure_path=None,
            structure_dir=None, structure_dir_depth=0, structure_suffix="",
            tmp_dir="./tmp",
            output_dir=None, output_format="json",
            update_interval=0.1):
        '''
        setup_method: Method that takes a structure, a sub tmp dir, and 
            setup_method_args to set up the structure specific tmp dir 
            for external evaluation to take place. Returns whether the 
            preparation is successful.
        restart_method: Method that takes a sub tmp dir and 
            restart_method_args to prepare an existing, executed sub tmp dir
            (marked by _evaluation_lock file with timestamp older than 
            60 s) for execute_method. Returns whether the preparation is 
            successful.
        execute_method: Method that takes an execute command, a sub tmp dir,
            and execute_method_args and returns a subprocess.Popen object
        extract_method: Method that takes a sub tmp dir where execute_method
            has completed and returns a single structure or (struct, result) 
            tuple
        '''
        # TODO: Implement retry that will prepare a subdir and place it
        # at the end of self._subdirs_to_evaluate list.
        self._execute_commands = execute_commands
        self._setup_method = setup_method
        self._execute_method = execute_method
        self._restart_method = restart_method
        self._extract_method = extract_method
        self._setup_method_args = setup_method_args

        self._structure_list = structure_list
        self._structure_path = structure_path
        self._structure_dir = structure_dir
        self._structure_dir_depth = structure_dir_depth
        self._structure_suffix = structure_suffix
        self._output_dir = output_dir
        self._output_format = output_format

        self._returned_result = []
        # stores tuple of (process, execute_command, subdir)
        self._executing_processes = []
        self._available_executing_commands = self._execute_commands[:]

    def setup(self):
        self._load_structures()
        struct_id_list = [x.struct_id for x in self._structure_list]
        struct_id_set = set(struct_id_list)
        if len(struct_id_list) != len(struct_id_set):
            raise RuntimeError("BatchSingleStrucutreExternalEvaluation "
                    "expects Structure IDs to be unique.")

        existing_subdirs = list_subdirectories(self._tmp_dir)
        _subdirs_to_evaluateself._subdirs_to_evaluate = []
        for subdir in [x for x in existing_subdirs if self._is_stale(x)]:
            if self._restart_method(subdir, *restart_method_args):
                print_time_log("Successfully prepares existing, executed "
                        "subdirectory %s for restart." % subdir)
                self._subdirs_to_evaluate.append(subdir)
            elif os.path.basename(subdir) in struct_id_set:
                print_time_log("Failed to prepare subdirectory %s for "
                        "restart. Removing this directory as newly loaded "
                        "structures contain the same basename as struct_id."
                        % subdir)
                shutil.rmtree(subdir)
            else:
                print_time_log("Failed to prepare subdirectory %s for "
                        "restart." % subdir)

        # Does not include subdirs that cannot be restarted
        # If any new structure has existing subdir that cannot be restarted,
        # the subdir would have been deleted at this point
        existing_subdir_struct_id = [os.path.basename(x)
                for x in self._subdirs_to_evaluate]
        for struct in self._structure_list:
            if struct.struct_id in existing_subdir_struct_id:
                print_time_log("Subdirectory already exists for structure %s; "
                        "skipping evaluation setup" % struct.struct_id)
                continue
            new_subdir = os.path.join(self._tmp_dir, struct.struct_id)
            if self._setup_method(struct, new_subdir, *setup_method_args):
                print_time_log("Successfully set up subdirectory for structure "
                        "%s at %s." % (self.struct_id, new_subdir))
                self._subdirs_to_evaluate.append(new_subdir)
            else:
                print_time_log("Failed to set up subdirectory for structure " +
                        self.struct_id)

    def execute(self):
        while (len(self._subdirs_to_evaluate) != 0 or
                len(self_executing_processes) != 0):
            self._check_processes()
            self._attempt_launch()
            time.sleep(update_interval)
 
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
        if self._structure_dir != None:
            self._load_structure_list()

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
                print_time_log(
                        "WARNING: Failed to load structure file: "
                        + path)
                return
                
        self._structure_list.append(struct)

    def _load_structure_list(self):
        # TODO: Merge with general load structure directory util
        for path in misc.list_directory_file(
                self._structure_dir,
                depth=self._structure_dir_depth,
                suffix=self._structure_suffix):
            self._load_structure(path)

    def _check_processes(self):
        not_completed = []
        for i in range(len(self._executing_processes)):
            process, execute_command, subdir = self._executing_processes[i]
            if process.poll() is None:
                self._mark_latest_execution_timestamp(subdir)
                # Process has not completed
                not_completed.append(i)
                continue
            print_time_log("Evaluation of subdir %s with execute command %s has "
                    "completed with exit code %s."
                    % (subdir, execute_command, str(process.returncode)))
            self._extract_subdir(subdir)
            self._mark_complete(subdir)
            self._available_execute_commands.append(execute_command)

        if len(not_completed) < len(self._executing_processes):
            self._executing_processes = [
                    self._executing_processes[x] for x in not_completed]

    def _extract_subdir(self, subdir):
        result = self._extract_method(subdir, *extract_method_args)
        if result is False:
            print_time_log("Extraction of execution-complete subdir %s failed."
                    % subdir)
            return
        if type(result) is tuple:
            _check_return_tuple(result)
            output_structure(result[0], self._output_dir, self._output_format)
        else:
            _check_return_single(result)
            output_structure(result, self._output_dir, self._output_format)

        self._returned_result.append(result)

    def _attempt_launch(self):
        if (len(self._available_execute_commands) == 0 or
                len(self._subdirs_to_evaluate) == 0):
            return
        execute_command = self._available_execute_commands.pop()
        subdir = self._subdirs_to_evaluate.pop()
        process = self._execute_method(execute_command, subdir,
                *execute_command_args)
        if process is False:
            print_time_log("Failed to launch evaluation in subdir %s."
                    % subdir)
            # Place execute_command at the front of queue to avoid
            # immediate reuse
            self._available_execute_commands = [execute_command] + \
                    self._available_execute_commands
            return
        print_time_log("Launched execution process to evaluate subdirectory "
                "%s with command %s" % (subdir, execute_command))
        self._mark_latest_execution_timestamp(subdir)
        self._executing_processes.append(process, execute_command, subdir)

    def _is_stale(self, subdir):
        if not os.path.isfile(os.path.join(subdir, "_evaluation_lock")):
            # Execution has not started
            return False
        last_updated_timestamp = retrieve_float(
                os.path.join(subdir, "_evaluation_lock"))
        if time.time() - last_updated_timestamp > 60:
            return True
        return False

    def _is_complete(self, subdir):
        return os.path.isfile(os.path.join(subdir, "_evaluation_complete"))

    def _mark_latest_execution_timestamp(self, subdir):
        f = open(os.path.join(subdir, "_evaluation_lock"), "w")
        f.write(time.time())
        f.close()

    def _mark_complete(self, subdir):
        # subdir can be removed by extract_method
        if os.path.isdir(subdir):
            safe_remove_files([os.path.join(subdir, "_evaluation_lock")])
            f = open(os.path.join(subdir, "_evaluation_complete"), "w")
            f.write("Evaluation complete at " + str(time.time()))
            f.close()
