import os, sys
from glob import glob
from Genarris.core.structure import Structure
from Genarris.utilities import file_utils
from Genarris.core.instruct import get_last_active_procedure_name, get_molecule_path
from Genarris.core.structure_handling import cell_niggli_reduction
from ibslib.io import read,write
import ase, spglib
import ase.io
from ase.spacegroup import Spacegroup
print('ase.__file__', ase.__file__, flush=True)
print('ase.io.__file__', ase.io.__file__, flush=True)
print('dir(ase.io)', dir(ase.io), flush=True)
import numpy as np
from copy import deepcopy


def update_lattice_lengths_in_struct(struct):
    if type(struct) is bool:
        return struct
    directions = ['a', 'b', 'c']
    for d in directions:
        if not 'lattice_vector_' + d in struct.properties:
            break
        struct.properties[d] = np.linalg.norm(
                    np.array(struct.properties['lattice_vector_' + d]))

    return struct

def check_type(var, desired_type):
    if desired_type == 'path':
        if (type(var) is not str) or not os.path.exists(var):
            raise TypeError('var', var, 'was not a found path')
        else:
            return
    if type(var) is not desired_type:
        raise TypeError('var', var, 'is not type', desired_type, 'as desired')

def run_aims(aims_lib_dir, comm, verbose):
    '''
    Run a single instance of aims.
    
    Arguments
    ---------
    comm: mpi4py.MPI object
        MPI communicator to pass into aims
    aims_lib_dir: str
        File path to the f2py-compiled aims_w.so library file. This file included
        libaims.XXXX.scalapack.mpi.so in its compilation.
    verbose: bool
        True: Print debugging output
        False: Do not print debugging output

    Returns
    -------
    None
    '''
    if aims_lib_dir not in sys.path:
        sys.path.append(aims_lib_dir)
    import aims_w

    rank = comm.rank
    commf = comm.py2f()
    if verbose:
        print('comm.rank', comm.rank, 'entering Barrier', flush=True)
    comm.Barrier()
    if verbose:
        print('comm.rank', comm.rank, 'about to run aims', flush=True)
    try:
        aims_w.aims_w(commf)
        if verbose:
            print('comm.rank', comm.rank, 'ran aims', flush=True)
    except:
        if verbose:
            print('comm.rank', comm.rank, 'could not run aims', flush=True)
    sys.stdout.flush()


def setup_aims_dirs(aims_output_dir, structure_dir, control_path):
    '''
    FHI-aims requires a new directory for every run. So, we need to
    create a directory for each run we anticipate doing: one per structure
    in structure_dir.
    
    Arguments
    ---------
    aims_output_dir: str
        Path to folder that will contain a folder for every structure you
        want to run FHI-aims on
    structure_dir: str
        Path to folder that contains every structure you want to run FHI-aims 
        on in json format.
    control_path: str
        path to FHI-aims control.in file

    
    Returns
    -------
    list
        Return task list which has a 0 for incomplete structures and 1 for 
        completed calculations for that structure.

    '''
    if os.path.isdir(aims_output_dir):
        #Get ready for restarting calculations
        task_list = []
        struct_folds = glob(os.path.join(aims_output_dir, '*/'))
        for struct_fold in struct_folds:
            aims_out = os.path.join(struct_fold, 'aims.out')
            geometry_in_next_step = os.path.join(struct_fold, 'geometry.in.next_step')
            if os.path.isfile(geometry_in_next_step):
                file_utils.cp(geometry_in_next_step, struct_fold, dest_fname='geometry.in')
            if len(file_utils.grep('Inconsistency of forces', aims_out)) > 0:
                file_utils.rm(aims_out)
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
    structure_files = glob(os.path.join(structure_dir, '*.json')) + glob(os.path.join(structure_dir, '*.in')) + glob(os.path.join(structure_dir, '*.in.next_step'))
    for structure_file in structure_files:
        name = file_utils.fname_from_fpath(structure_file)
        struct_fold = os.path.join(aims_output_dir, name)
        file_utils.mkdir_if_DNE(struct_fold)
        struct = Structure()
        if '.json' in structure_file:
            struct.build_geo_from_json_file(structure_file)
        elif '.in' in structure_file or '.in.next_step' in structure_file:
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
        print('wanted one and only one match to "| Total energy of the DFT"'+
              'but didnt get that for ' + aims_out, flush=True)
        return 'none'
    energy = float(grep_results[0].split()[11])
    return energy


def run_fhi_aims_batch(comm, world_comm, MPI_ANY_SOURCE, num_replicas, inst=None, sname=None, structure_dir=None, aims_output_dir=None, output_dir=None, aims_lib_dir=None, control_path=None, energy_name='energy', verbose=False):
    """ Performs multiple FHI calculations
    
    Arguments
    ---------
    comm: mpi4py.MPI object
        MPI communicator to pass into aims
    world_comm: mpi4py.MPI object
        World MPI communicator
    MPI_ANY_SOURCE: mpi4py.MPI.ANY_SOURCE
        Any source object for communication.
    num_replicas: int
        Number of replicas to perform calculation.
    inst: genarris.core.instruct.Instruct
        Config Parser object which contains all the configuration file sections
        and options for calculation.
    sname: str
        Section name which called run_fhi_aims_batch
    struct_dir: str
        Path to directory of structures to perform calculation.
    aims_output_dir: str
        Path to the directory where FHI-aims calculations should take place.
    output_dir: str
        Path to the directory where the Structure files should be saved.
    aims_lib_dir: str
        Path to the directory containing the FHI-aims library file.
    control_path: str
        Path to the directory containing the control file to use.
    energy_name : str
        Property name which the calculated energy will be stored with in the 
        Structure file.
    verbose : bool
        Controls verbosity of output.
        
    
    Configuration File Options
    --------------------------
    verbose : bool
        Controls verbosity of output.
    energy_name : str
        Property name which the calculated energy will be stored with in the 
        Structure file.
    output_dir : str
        Path to the directory where the output structure file will be saved.
    aims_output_dir : str
        Path where the aims calculation will take place.
    aims_lib_dir : str
        Path to the location of the directory containing the FHI-aims library 
        file. 
    molecule_path : str
        Path to the geometry.in file of the molecule to be calculated if 
        called using harris_single_molecule_prep or relax_single_molecule.
    structure_dir : str
        Path to the directory of structures to be calculated if calculation
        was called not using harris_single_molecule_prep or 
        relax_single_molecule.
        
        
    Returns
    -------
    None : None
        
    
    """
    if aims_output_dir is None:
        aims_output_dir = os.getcwd()
    if inst is not None:
        verbose = inst.get_boolean(sname, 'verbose')
        energy_name = inst.get_with_default(sname, 'energy_name', 'energy')
        output_dir = inst.get_or_none(sname, 'output_dir')
        aims_output_dir = inst.get(sname, 'aims_output_dir')
        sname_list = [sname, 'relax_single_molecule', 'fhi_aims_energy_evaluation', 'harris_single_molecule_prep', 'harris_approximation_batch', 'run_fhi_aims_batch']
        aims_lib_dir = inst.get_inferred(sname, sname_list, ['aims_lib_dir'] * 6, type_='dir', required=True)

        molecule_path = get_molecule_path(inst, sname)

        if sname == 'harris_single_molecule_prep' or sname == 'relax_single_molecule':
            if inst.has_option(sname, 'structure_dir'):
                structure_dir = inst.get(sname, 'structure_dir')
            elif sname == 'harris_single_molecule_prep':
                structure_dir = 'structure_dir_for_harris_single_molecule_prep'
            else:
                structure_dir = 'structure_dir_for_relaxing_single_molecule'
            if not os.path.isfile(molecule_path):
                raise Exception("molecule_path", molecule_path, 'DNE')
            if not os.path.isfile(os.path.join(structure_dir, os.path.basename(molecule_path))) and world_comm.rank == 0:
                file_utils.cp(molecule_path, structure_dir)
        else:
            last_section = get_last_active_procedure_name(inst, sname)
            sname_list = [sname, last_section, last_section, last_section, 'affinity_propagation_fixed_clusters', 'affinity_propagation_fixed_clusters', 'rcd_calculation', 'harris_approximation_batch', 'pygenarris_structure_generation', 'structure_generation_batch'] * 2
            structure_dir = inst.get_inferred(sname, sname_list,
                                                    ['structure_dir', 'exemplars_output_dir_2', 'exemplars_output_dir', 'output_dir', 'exemplars_output_dir_2', 'exemplars_output_dir', 'output_dir', 'output_dir', 'output_dir', 'output_dir'] + (len(sname_list) // 2) * ['structure_dir'], type_='dir', required=False)

        control_path = inst.get(sname, 'control_path')
        Z = int(inst.get_inferred(sname, [sname, 'pygenarris_structure_generation', 'estimate_unit_cell_volume'],
                                        ['Z']*3))

    if world_comm.rank == 0:
        print('run_fhi_aims_batch using molecule_path', molecule_path, flush=True)
        # creates the working dirs if DNE. Gets them ready to restart otherwise.
        task_list = setup_aims_dirs(aims_output_dir, structure_dir, control_path)
        if verbose:
            print('task_list',task_list, flush=True)
            print('aims_output_dir', aims_output_dir, flush=True)
            print('structure_dir', structure_dir,flush=True)
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
                run_aims(aims_lib_dir, comm, verbose)
                if comm.rank == 0:
                    # Update struct json with energy
                    structure_files = glob('*.json')
                    if len(structure_files) == 0:
                        struct = read('geometry.in', file_format='aims')
                        name = file_utils.fname_from_fpath(molecule_path)
                        structure_file = name + '.json'
                        struct.struct_id = name
                        write(structure_file, struct)
                    else:
                        structure_file = structure_files[0]
                    struct = read(structure_file)
                    struct_is_periodic = 'lattice_vector_a' in struct.properties
                    # update json with new geometry if geometry.in.next_step exists
                    geometry_in_next_step = 'geometry.in.next_step'
                    if verbose:
                        print('checking if', geometry_in_next_step, 'exists', flush=True)
                    if os.path.exists(geometry_in_next_step):
                        if verbose:
                            print(geometry_in_next_step, 'exists', flush=True)
                        file_utils.cp(geometry_in_next_step, '.', dest_fname='geometry.in')
                        geometry_in_next_step = 'geometry.in'
                        relaxed_struct = read(geometry_in_next_step, file_format='aims')
                        relaxed_struct.struct_id = struct.struct_id
                        if verbose:
                            print('struct_id', relaxed_struct.struct_id, flush=True)
                            print('Z', Z, flush=True)
                        if not struct_is_periodic:
                            napm = None
                        else:
                            napm = len(relaxed_struct.geometry) / Z
                        if verbose:
                            print('napm', napm, 'Z', Z, flush=True)
                        if struct_is_periodic:
                            relaxed_struct = cell_niggli_reduction(relaxed_struct, napm, create_duplicate=False)
                            relaxed_struct = update_lattice_lengths_in_struct(relaxed_struct)
                            if verbose:
                                print('updated lattice lengths', flush=True)
                        for struct_property in struct.properties:
                            if struct_property not in relaxed_struct.properties:
                                relaxed_struct.properties[struct_property] = struct.properties[struct_property]
                        if verbose:
                            print('updated struct properties', flush=True)
                        write(geometry_in_next_step, relaxed_struct, file_format='aims', overwrite=True)
                        write('geometry.in.next_step', relaxed_struct, file_format='aims', overwrite=True)

                        if verbose:
                            print('relaxed_struct.geometry',relaxed_struct.geometry, flush=True)
                            print('updated', geometry_in_next_step, flush=True)
                        if struct_is_periodic:
                            if verbose:
                                print('structure is periodic', flush=True)
                            if not os.path.exists(geometry_in_next_step):
                                raise Exception('path', geometry_in_next_step, 'DNE')
                            else:
                                print(geometry_in_next_step, 'exists', flush=True)
                            try: ase.io.read
                            except: print('ase.io.read not defined')
                            ase_struct = ase.io.read(geometry_in_next_step, format='aims', parallel=False)
                            if verbose:
                                print('read ase_struct', ase_struct, flush=True)
                            lattice = ase_struct.get_cell()
                            if verbose:
                                print('ase lattice', lattice, flush=True)
                            positions = ase_struct.get_scaled_positions()
                            if verbose:
                                print('ase positions', positions, flush=True)
                            numbers = ase_struct.get_atomic_numbers()
                            if verbose:
                                print('ase numbers', numbers, flush=True)
                            cell = (lattice, positions, numbers)
                            if verbose:
                                print('cell', cell, flush=True)
                            spglib_data = spglib.get_symmetry_dataset(cell, symprec=1e-1)
                            spglib_spg = int(spglib_data['number'])
                            if verbose:
                                print('spglib_spg', spglib_spg, flush=True)
                            multiplicity = Spacegroup(spglib_spg)
                            if multiplicity.nsymop == Z:
                                # general position
                                site_symmetry = 1
                            else:
                                # special position
                                site_symmetry = 0
                            if verbose:
                                print('site_symmetry', site_symmetry, flush=True)
                            if 'spg' in relaxed_struct.properties:
                                relaxed_struct.properties['spg_before_relaxation'] = relaxed_struct.properties['spg']
                            relaxed_struct.properties['spg'] = spglib_spg
                            if 'cell_vol' in relaxed_struct.properties:
                                relaxed_struct.properties['unit_cell_volume'] = relaxed_struct.properties['cell_vol']
                        relaxed_struct.struct_id = struct.struct_id
                        if verbose:
                            print('relaxed struct_id', relaxed_struct.struct_id, flush=True)

                        if struct_is_periodic and 'site_symmetry_group' in relaxed_struct.properties and 'site_symmetry_group_before_relaxation' not in relaxed_struct.properties:
                            relaxed_struct.properties['site_symmetry_group_before_relaxation'] = relaxed_struct.properties['site_symmetry_group']
                        if struct_is_periodic:
                            relaxed_struct.properties['site_symmetry_group'] = site_symmetry
                            if verbose:
                                print('updated symmetry properties', flush=True)
                        relaxed_struct.properties[energy_name] = extract_energy(aims_out)
                        if verbose:
                            print('updated energy', flush=True)
                        write(structure_file, relaxed_struct, overwrite=True)
                        
                    else:
                        struct.properties[energy_name] = extract_energy(aims_out)
                        write(structure_file, struct, overwrite=True)

                os.chdir('../..')

            if comm.rank == 0:
                world_comm.send([world_comm.rank, message], dest=0, tag=0)
                message = world_comm.recv(source=0, tag=1)
            else:
                message = None
            message = comm.bcast(message, root=0)

    if world_comm.rank == 0:
        # Output (copy) jsons to output_dir
        if verbose:
            print('output_dir', output_dir, flush=True)
        if output_dir is not None:
            json_flist = file_utils.glob(os.path.join(aims_output_dir, '**', '*.json'), recursive=True)
            if verbose:
                print('aims_output_dir', aims_output_dir, flush=True)
                print('json_flist', json_flist, flush=True)
            # Only save those structures into jsons in output_dir which successfully gave an energy according to the desired control file
            for json_fpath in json_flist:
                name = file_utils.fname_from_fpath(json_fpath)
                aims_out = os.path.join(aims_output_dir, name, 'aims.out')
                if verbose:
                    print('aims_out', aims_out, flush=True)
                    print('sname', sname, flush=True)
                struct = read(json_fpath)
                if energy_name in struct.properties and struct.properties[energy_name] != 'none' and \
                    (sname != 'harris_approximation_batch' or len(file_utils.grep('Reading periodic restart information from file', aims_out)) > 0):
                    
                    if verbose:
                        print('file_utils.grep(Reading periodic restart information from file, aims_out)', file_utils.grep('Reading periodic restart information from file', aims_out), flush=True)
                    file_utils.cp(json_fpath, output_dir)
                elif verbose:
                    print('Not copying', json_fpath, 'to output_dir because energy was not obtained by aims', flush=True)
        plot_histograms = inst.get_boolean(sname, 'plot_histograms')
        if plot_histograms:
            prop = inst.get_with_default(sname, 'prop', 'unit_cell_volume')
            prop_figname = inst.get_with_default(sname, 'prop_figname', 'raw_pool_volume_histogram.pdf')
            prop_xlabel = inst.get_with_default(sname, 'prop_xlabel', 'Structure Volume, $\AA^3$')
            prop_ylabel = inst.get_with_default(sname, 'prop_ylabel', 'Counts')
            prop_figure_size = inst.get_with_default(sname, 'prop_figure_size', (12,8), eval=True)
            prop_label_size = inst.get_with_default(sname, 'prop_label_size', 24, eval=True)
            prop_tick_size = inst.get_with_default(sname, 'prop_tick_size', 18, eval=True)
            prop_tick_width = inst.get_with_default(sname, 'prop_tick_width', 3, eval=True)
            prop_GAtor_IP = inst.get_boolean(sname, 'prop_GAtor_IP')

            pygenarris_outfile = inst.get_with_default(sname, 'pygenarris_outfile', 'outfile')
            spg_bar_chart_fname = inst.get_with_default(sname, 'spg_bar_chart_fname', 'raw_pool_spg_bar_chart.pdf')
            spg_bar_width = inst.get_with_default(sname, 'spg_bar_width', 0.5, eval=True)
            spg_bar_xlabel = inst.get_with_default(sname, 'spg_bar_xlabel', 'Allowed space groups')
            spg_bar_ylabel = inst.get_with_default(sname, 'spg_bar_ylabel', 'Count')
            spg_bar_title = inst.get_or_none(sname, 'spg_bar_title')
            spg_bar_tick_rotation = inst.get_with_default(sname, 'spg_bar_tick_rotation', 'vertical')

            plot_property(output_dir, prop=prop, nmpc=Z, figname=prop_figname, 
                                xlabel=prop_xlabel, ylabel=prop_ylabel, figure_size=prop_figure_size,
                                label_size=prop_label_size, tick_size=prop_tick_size, tick_width=prop_tick_width, GAtor_IP=prop_GAtor_IP)

            plot_spg_bar_chart(output_dir, pygenarris_outfile=pygenarris_outfile, spg_bar_chart_fname=spg_bar_chart_fname,
                                width=spg_bar_width, ylabel=spg_bar_ylabel, xlabel=spg_bar_xlabel,
                                title=spg_bar_title, tick_rotation=spg_bar_tick_rotation)


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
