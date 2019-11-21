import sys, os, random
import numpy as np
import _pygenarris_mpi as pygenarris_mpi
import _pygenarris as pygenarris
import math
from Genarris.utilities import file_utils, list_utils, time_utils
from Genarris.core import structure
from Genarris.core import file_handler
from Genarris.core.instruct import get_last_active_procedure_name, get_molecule_path
from Genarris.utilities.find_bonding import MoleculeBonding
from Genarris.core.structure_handling import cell_niggli_reduction
from ibslib.io import read,write

def write_json_files( output_dir, structs_list, Z, output_format):
    """
    function to create jsons from pygenarris geometry.out

    """
    if output_format == 'json' or output_format == 'both':
        print('writing json files...', flush=True)
        file_utils.mkdir_if_DNE(output_dir)

    geometry_out_fpath = "geometry.out"
    outfile_basename = os.path.basename(geometry_out_fpath)
    outfile_dirname = os.path.dirname(os.path.abspath(geometry_out_fpath))
    tmp_geo_fname = 'tmp_file_for_geo.in' # must have .in extension
    if not type(Z) is int:
        Z = 'none'
    
    for i, struct_lines in enumerate(structs_list):
        number_of_atoms_in_molecule, unit_cell_volume, attempted_spg, attempted_wyckoff_position, site_symmetry_group, spg, first_geo_line = \
            'none', 'none', 'none', 'none', 'none', 'none', 'none'
        
        for j,line in enumerate(struct_lines):
            if '#Z =' in line:
                Z = int(line.split('=')[-1].split()[0])
            elif '#number_of_atoms_in_molecule =' in line:
                number_of_atoms_in_molecule = int(line.split('=')[-1].split()[0])
            elif '#unit_cell_volume =' in line:
                unit_cell_volume = float(line.split('=')[-1].split()[0])
            elif '#attempted_spacegroup =' in line:
                attempted_spg = int(line.split('=')[-1].split()[0])
            elif '#attempted_wyckoff_position =' in line:
                attempted_wyckoff_position = line.split('=')[-1].split()[0]
            elif '#site_symmetry_group =' in line:
                site_symmetry_group = line.split('=')[-1].split()[0]
            elif '#SPGLIB_detected_spacegroup =' in line:
                spg = int(line.split('=')[-1].split()[0])
            elif first_geo_line == 'none' and 'lattice_vector' == line.split()[0]:
                first_geo_line = j
            elif first_geo_line == 'none' and 'atom' == line.split()[0]:
                first_geo_line = j

        #write to one structure to file "tmp_file_for_geo.in" and read from it.    
        file_utils.write_lines_to_file(tmp_geo_fname, struct_lines, mode='w')
        struct = read(tmp_geo_fname, file_format='aims')
        struct.properties['Z'] = Z
        struct.properties['number_of_atoms_in_molecule'] = number_of_atoms_in_molecule
        struct.properties['unit_cell_volume'] = unit_cell_volume
        struct.properties['attempted_spg'] = attempted_spg
        struct.properties['attempted_wyckoff_position'] = attempted_wyckoff_position
        struct.properties['site_symmetry_group'] = site_symmetry_group
        struct.properties['spg'] = spg
        random_str = file_handler.get_random_index()
        struct_id = file_utils.fname_from_fpath(geometry_out_fpath) + '_' + random_str
        struct.struct_id = struct_id
        outfile_path = os.path.join(outfile_dirname, output_dir, struct_id + '.json')
        if 'lattice_vector_a' in struct.properties:
            napm = len(struct.geometry) / Z
            struct = cell_niggli_reduction(struct, napm, create_duplicate=False)
            write(tmp_geo_fname, struct, file_format='aims', overwrite=True)
            geo_lines = file_utils.get_lines_of_file(tmp_geo_fname)
            structs_list[i] = structs_list[i][:first_geo_line] + geo_lines
        if output_format == 'json' or output_format == 'both':
            write(outfile_path, struct, file_format='json', overwrite=True)
        structs_list[i][1] = '#structure_number = ' + str(i) + '\n'
        structs_list[i] = structs_list[i][:2] + ['#struct_id = ' + struct_id + '\n'] + structs_list[i][2:]
        first_geo_line = 'none'
    os.remove(tmp_geo_fname)


def set_up(molecule_path):
    # pygenarris requires that "geometry.in" be in the current working directory and be the single molecule geometry in aims format.
    if not os.path.isfile(molecule_path):
        raise Exception('molecule_path path DNE. molecule_path = ', molecule_path)
    
    file_utils.cp(molecule_path, '.', dest_fname='geometry.in')


def format_output(output_format, output_dir,  num_structures, Z):
    '''
    After all structures have been output by various cores into a file, collect all the 
    structures into a usable format: output a single file which is fhi-aims 
    geometry format if "geometry" or "both" is output_format or a folder of 
    json files (one for each structure) if "json" or "both" is output_format.

    
    Arguments
    ---------        
    output_format: str
        "json" for outputting a folder of json files - each containing a structure
        "geometry" for outputting a single geometry file where each structure's fhi-aims-formatted structure
            is concatenated.
        "both" for outputting according to "json" and "geometry"

    output_dir: str
        path of folder to dump the json files into


    num_structures: int or None
        Number of structures in the final raw pool or if None, output them all


    Z: int
        Number of molecules per unit cell

    Returns
    -------
    None : None
    
    '''
    if not (output_format == 'json' or output_format == 'geometry' or output_format == 'both'):
        raise Exception('Only "json", "geometry", or "both" are supported output formats. output_format = ', output_format)

    if os.path.isdir(output_dir):
        file_utils.rm(output_dir)

    lines_list = file_utils.get_lines_of_file("geometry.out")
    
    indices_list = file_utils.grep('BEGIN STRUCTURE', "geometry.out", return_line_nums=True)[1] + [len(lines_list)]
    structs_list = [
                    lines_list[
                                indices_list[i]
                                :
                                indices_list[i + 1]
                                ]
                    for i in range(len(indices_list) - 1)
                    ]
    
    structs_list = write_json_files(output_dir, structs_list, Z, output_format)
    

def clean_up(output_format):

    if output_format == 'json':
        os.remove("geometry.out")

def pygenarris_structure_generation(inst=None, comm=None, 
        num_structures=None, Z=None, volume_mean=None, 
        volume_std=None, sr=None, tol=None, max_attempts_per_spg=None, 
        molecule_path=None):
    """
    Passes in all necessary arguments for generation to Pygenarris.
    """
    # Currently does not support multiple instances running simultaneously. 
    # If you want more simultaneous processes, submit more MPI ranks.
    if comm is not None:
        comm.barrier()
    if inst is not None:
        
        sname = 'pygenarris_structure_generation'
        num_structures = inst.get_or_none(sname, 'num_structures', eval=True)
        Z = inst.get_eval(sname, 'Z')
        volume_mean = inst.get_eval(sname, 'volume_mean')
        volume_std = inst.get_eval(sname, 'volume_std')
        sr = inst.get_eval(sname, 'sr')
        tol = inst.get_eval(sname, 'tol')
        max_attempts_per_spg = inst.get_eval(sname, 'max_attempts_per_spg')
        molecule_path = get_molecule_path(inst, sname)
        output_format = inst.get_with_default(sname, 'output_format', 'json') #options are json, geometry, both
        output_dir = inst.get_with_default(sname, 'output_dir', '.')

    if comm is None:
        set_up(molecule_path)
    else:
        if comm.rank == 0:
            print('molecule path in Pygenarris Structure Generation:', molecule_path, flush=True)
            set_up(molecule_path)
        comm.barrier()

        
    molecule_struct = read(molecule_path)
    mb = MoleculeBonding(molecule_struct)
    cutoff_matrix = mb.get_crystal_cutoff_matrix(Z, vdw_mult=sr)
    cutoff_matrix = np.array(cutoff_matrix, dtype='float32')

    if num_structures is None:
        num_structures_per_allowed_SG = 1
    else:
        num_compatible_spg = pygenarris.num_compatible_spacegroups(Z, tol)
        num_structures_per_allowed_SG = math.ceil(num_structures / float(num_compatible_spg))
    if comm is not None:
        if comm.rank==0:
             print('number of structures per allowed space group:', num_structures_per_allowed_SG, flush=True)
             print('Z:', Z, flush=True)
             print('volume mean:', volume_mean, flush=True)
             print('volume standard deviation:', volume_std, flush=True)
             print('tol:', tol, flush=True)
             print('maximum attempts per space group:', max_attempts_per_spg, flush=True)
             print("************iterations for each space group is listed in detail in standard out file*************")
             np.savetxt('cutoff matrix:', cutoff_matrix)
    
    start_pygenarris_time = time_utils.gtime()
    if comm is not None:
        comm.barrier()
    pygenarris_mpi.mpi_generate_molecular_crystals_with_vdw_cutoff_matrix(cutoff_matrix,
                   num_structures_per_allowed_SG, Z, volume_mean, volume_std, tol,
                   max_attempts_per_spg, comm)
    '''
    if comm is not None:
        print('Time for just pygenarris =', time_utils.gtime() - start_pygenarris_time, 'comm.rank', comm.rank,  flush=True)
    '''
    if comm is not None:
        comm.barrier()
        if comm.rank == 0:
            format_output(output_format, output_dir, num_structures,  Z)
            clean_up(output_format)

