import os, sys
from glob import glob
from core.structure import Structure
from utilities import file_utils

def check_type(var, desired_type):
    if desired_type == 'path':
        if (type(var) is not str) or not os.path.exists(var):
            raise TypeError('var', var, 'was not a found path')
        else:
            return
    if type(var) is not desired_type:
        raise TypeError('var', var, 'is not type', desired_type, 'as desired')

def run_aims(aims_lib_dir, comm):
        '''
        comm: mpi4py.MPI object
            MPI communicator to pass into aims
        aims_lib_dir: str
            File path to the f2py-compiled aims_w.so library file. This file included
            libaims.XXXX.scalapack.mpi.so in its compilation.

        Purpose:
            Run a single instance of aims.
        '''
        if aims_lib_dir not in sys.path:
            sys.path.append(aims_lib_dir)
        import aims_w

        rank = comm.rank
        commf = comm.py2f()
        print('comm.rank', comm.rank, 'entering Barrier', flush=True)
        comm.Barrier()
        print('comm.rank', comm.rank, 'about to run aims', flush=True)
        try:
            aims_w.aims_w(commf)
            print('comm.rank', comm.rank, 'ran aims', flush=True)
        except:
            print('comm.rank', comm.rank, 'could not run aims', flush=True)
        sys.stdout.flush()


def setup_aims_dirs(aims_output_dir, structure_dir, control_path):
    '''
    aims_output_dir: str
        Path to folder that will contain a folder for every structure you
        want to run FHI-aims on
    structure_dir: str
        Path to folder that contains every structure you want to run FHI-aims on
        in json format.
    control_path: str
        path to FHI-aims control.in file

    Return: list
        Return task list which has a 0 for incomplete structures and 1 for completed
        calculations for that structure.

    Purpose: FHI-aims requires a new directory for every run. So, we need to
        create a directory for each run we anticipate doing: one per structure
        in structure_dir.
    '''
    
    if os.path.isdir(aims_output_dir):
        #Get ready for restarting calculations
        task_list = []
        struct_folds = glob(os.path.join(aims_output_dir, '*/'))
        for struct_fold in struct_folds:
            aims_out = os.path.join(struct_fold, 'aims.out')
            geometry_in = os.path.join(struct_fold, 'geometry.in')
            geometry_in_next_step = os.path.join(struct_fold, 'geometry.in.next_step')
            if os.path.isfile(geometry_in_next_step):
                file_utils.cp(geometry_in_next_step, geometry_in)
            if os.path.isfile(aims_out):
                # If still running or didn't converge then don't rerun it.
                task_list.append(1)
            else:
                task_list.append(0)
        print('Ready for restart', flush=True)
        return task_list

    print('Setting working directories up', flush=True)
    # Setup aims working dirs
    os.makedirs(aims_output_dir)
    structure_files = glob(os.path.join(structure_dir, '*.json')) + glob(os.path.join(structure_dir, '*.in'))
    for structure_file in structure_files:
        name = file_utils.fname_from_fpath(structure_file)
        struct_fold = os.path.join(aims_output_dir, name)
        file_utils.mkdir_if_DNE(struct_fold)
        struct = Structure()
        if '.json' in structure_file:
            struct.build_geo_from_json_file(structure_file)
        elif '.in' in structure_file:
            struct.build_geo_from_atom_file(structure_file)
            file_utils.write_to_file(os.path.join(struct_fold, name + '.json'), struct.dumps(), mode='w')
            #struct_dct = file_utils.get_dct_from_json(structure_file)
            #file_utils.write_dct_to_json(os.path.join(struct_fold, name + '.json'), struct_dct)
        
        file_utils.cp(structure_file, struct_fold)
        geometry_in = os.path.join(struct_fold, 'geometry.in')
        file_utils.write_to_file(geometry_in, str(struct.get_geometry_atom_format()))
        file_utils.cp(control_path, struct_fold, dest_fname='control.in')
    print('Ready for start', flush=True)
    return [0] * len(structure_files)


def extract_energy(aims_out):
    search_str = '| Total energy of the DFT'
    grep_results = file_utils.grep(search_str, aims_out)
    if len(grep_results) > 1 or len(grep_results) == 0:
        print('wanted one and only one match to "| Total energy of the DFT" but didnt get that for ' + aims_out, flush=True)
        return 'none'
    energy = float(grep_results[0].split()[11])
    return energy


def run_fhi_aims_batch(comm, world_comm, MPI_ANY_SOURCE, num_replicas, inst=None, sname=None, structure_dir=None, aims_output_dir=None, output_dir=None, aims_lib_dir=None, control_path=None, energy_name='energy', verbose=False):
    if aims_output_dir is None:
        aims_output_dir = os.getcwd()
    if inst is not None:
        verbose = inst.get_boolean(sname, 'verbose')
        energy_name = inst.get(sname, 'energy_name')
        output_dir = inst.get_or_none(sname, 'output_dir')
        aims_output_dir = inst.get(sname, 'aims_output_dir')
        sname_list = ['harris_single_molecule_prep', 'harris_approximation_batch', 'run_fhi_aims_batch']
        aims_lib_dir = inst.get_inferred(sname, sname_list, ['aims_lib_dir'] * 3)

        sname_list = ['harris_single_molecule_prep', 'estimate_unit_cell_volume', 'pygenarris_structure_generation', 'structure_generation_batch', 'harris_approximation_batch']
        molecule_path = inst.get_inferred(sname, sname_list, ['molecule_path'] * 5)

        if sname == 'harris_single_molecule_prep':
            if inst.has_option(sname, 'structure_dir'):
                structure_dir = inst.get(sname, 'structure_dir')
            else:
                structure_dir = 'structure_dir_for_harris_single_molecule_prep'
            if not os.path.isfile(molecule_path):
                raise Exception("molecule_path", molecule_path, 'DNE')
            if not os.path.isfile(os.path.join(structure_dir, os.path.basename(molecule_path))) and world_comm.rank == 0:
                file_utils.cp(molecule_path, structure_dir)
        else:
            sname_list = [sname, 'affinity_propagation_fixed_clusters', 'affinity_propagation_fixed_clusters', 'rcd_calculation', 'harris_approximation_batch', 'pygenarris_structure_generation', 'structure_generation_batch'] * 2
            structure_dir = inst.get_inferred(sname, sname_list,
                                                    ['structure_dir', 'exemplars_output_dir_2', 'exemplars_output_dir', 'output_dir', 'output_dir', 'output_dir', 'output_dir'] + (len(sname_list) // 2) * ['structure_dir'], required=False)

        control_path = inst.get(sname, 'control_path')

    if world_comm.rank == 0:
        # creates the working dirs if DNE. Gets them ready to restart otherwise.
        task_list = setup_aims_dirs(aims_output_dir, structure_dir, control_path)
        if verbose:
            print('world_comm.size', world_comm.size, flush=True)

    if world_comm.rank == 0:
        done_ranks = []
        while len(done_ranks) != num_replicas:
            if verbose:
                print('world rank 0 waiting for ready task', flush=True)
            ready_rank, task_id = world_comm.recv(source=MPI_ANY_SOURCE, tag=0)
            if verbose:
                print('world rank 0 received task_id', task_id, ' from ready rank', ready_rank, flush=True)
            if task_id != 'starting':
                task_list[task_id] = 1
            if 0 in task_list:
                message = task_list.index(0)
                task_list[message] = 1
            else: #Done!
                message = 'done'
                done_ranks.append(ready_rank)
            if verbose:
                print('world rank 0 sending message', message, 'to rank' , ready_rank, flush=True)
            world_comm.send(message, dest=ready_rank, tag=1)
            if verbose:
                print('world rank 0 sent message', message, 'to rank' , ready_rank, flush=True)
    else:
        if comm.rank == 0:
            if verbose:
                print('According to comm.rank 0, comm.size is ', comm.size, flush=True)
                print('According to comm.rank 0, world_comm.size is ', world_comm.size, flush=True)
            world_comm.send([world_comm.rank, 'starting'], dest=0, tag=0)
            if verbose:
                print('comm.rank 0 sent ',[world_comm.rank, 'starting'], ' to world_comm.rank 0', flush=True)
            message = world_comm.recv(source=0, tag=1)
            if verbose:
                print('comm.rank 0 received message', message, flush=True)
        else:
            message = None
        message = comm.bcast(message, root=0)
        if verbose:
            print('world_comm.rank', world_comm.rank, 'has message', message, flush=True)
        struct_folds = glob(os.path.join(aims_output_dir, '*/'))
        while message != 'done':
            struct_fold = struct_folds[message]
            check_type(aims_lib_dir, 'path')
            check_type(struct_fold, 'path')
            aims_out = os.path.abspath(os.path.join(struct_fold, 'aims.out'))
            if not os.path.isfile(aims_out):
                os.chdir(struct_fold)
                run_aims(aims_lib_dir, comm)
                if comm.rank == 0:
                    # Update struct json with energy
                    structure_files = glob('*.json')
                    if len(structure_files) == 0:
                        struct = Structure()
                        struct.build_geo_from_atom_file('geometry.in')
                        name = file_utils.fname_from_fpath(molecule_path)
                        structure_file = name + '.json'
                        file_utils.write_to_file(structure_file, struct.dumps(), mode='w')
                    else:
                        structure_file = structure_files[0]
                    struct_dct = file_utils.get_dct_from_json(structure_file)
                    struct_dct['properties'][energy_name] = extract_energy(aims_out)
                    file_utils.write_dct_to_json(structure_file, struct_dct)

                os.chdir('../..')

            if comm.rank == 0:
                world_comm.send([world_comm.rank, message], dest=0, tag=0)
                message = world_comm.recv(source=0, tag=1)
            else:
                message = None
            message = comm.bcast(message, root=0)

    if world_comm.rank == 0:
        # Output (copy) jsons to output_dir
        if output_dir is not None:
            json_flist = file_utils.glob(os.path.join(aims_output_dir, '**', '*.json'), recursive=True)
            file_utils.cp(json_flist, output_dir)



if __name__ == '__main__':
    from mpi4py import MPI
    aims_lib_dir = "/home/trose/aims_python_lib/cmake_190527"
    comm = MPI.COMM_WORLD
    structure_dir = 'structure_dir'
    control_path = 'control.in'
    sname = 'run_fhi_aims_batch'
    verbose = False
    aims_output_dir = os.getcwd()
    output_dir = 'output_dir'

    run_fhi_aims_batch(comm, comm, MPI.ANY_SOURCE, num_replicas=1, sname=sname, structure_dir=structure_dir, aims_output_dir=aims_output_dir, output_dir=output_dir, aims_lib_dir=aims_lib_dir, control_path=control_path, energy_name='energy', verbose=verbose)
