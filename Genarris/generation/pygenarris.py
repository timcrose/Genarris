import sys, os
import numpy as np
import pygenarris
from Genarris.utilities import file_utils, list_utils
from Genarris.core import structure
from Genarris.core import file_handler
from Genarris.utilities.find_bonding import MoleculeBonding
from ibslib.io import read

def write_json_files(geometry_out_fpath, output_dir, structs_list):
    file_utils.mkdir_if_DNE(output_dir)

    outfile_basename = os.path.basename(geometry_out_fpath)
    outfile_dirname = os.path.dirname(os.path.abspath(geometry_out_fpath))
    tmp_geo_fname = 'tmp_file_for_geo'
    for i, struct_lines in enumerate(structs_list):
        properties_list = struct_lines[1].split()
        spg = int(properties_list[properties_list.index('spg') + 2])
        file_utils.write_lines_to_file(tmp_geo_fname, struct_lines, mode='w')
        struct = structure.Structure()
        struct.build_geo_from_atom_file(tmp_geo_fname)
        struct.properties['spg'] = spg
        random_str = file_handler.get_random_index()
        struct_id = file_utils.fname_from_fpath(geometry_out_fpath) + '_' + random_str
        struct.struct_id = struct_id
        outfile_path = os.path.join(outfile_dirname, output_dir, struct_id + '.json')
        file_utils.write_to_file(outfile_path, struct.dumps(), mode='w')
        structs_list[i][0] = '#structure number = ' + str(i) + '\n'
        structs_list[i] = structs_list[i][:2] + ['#struct_id = ' + struct_id + '\n'] + structs_list[i][2:]
    os.remove(tmp_geo_fname)
    return structs_list


def set_up(molecule_path):
    # pygenarris requires that "geometry.in" be in the current working directory and be the single molecule geometry in aims format.
    if not os.path.isfile(molecule_path):
        raise Exception('molecule_path path DNE. molecule_path = ', molecule_path)
    
    file_utils.cp(molecule_path, '.', dest_fname='geometry.in')


def format_output(output_format, output_dir, final_filename, num_structures):
    '''
    output_format: str
        "json" for outputting a folder of json files - each containing a structure
        "geometry" for outputting a single geometry file where each structure's fhi-aims-formatted structure
            is concatenated.
        "both" for outputting according to "json" and "geometry"

    output_dir: str
        path of folder to dump the json files into
    
    final_filename: str
        File name used to derive the file name for each core. It is the file name that
        the concatenated file (if output_format != json) will have.

    num_structures: int
        Number of structures in the final raw pool.

    Return: None

    Purpose: After all structures have been output by various cores into a file (with a filename being a variant of
        final_filename), collect all the structures into a usable format: output a single file which is fhi-aims geometry 
        format if "geometry" or "both" is output_format or a folder of json files (one for each structure) if "json" or "both" is output_format.
        Pick a random subset of structures from the ones generated with length num_structures and remove the rest.
    '''
    if not (output_format == 'json' or output_format == 'geometry' or output_format == 'both'):
        raise Exception('Only "json", "geometry", or "both" are supported output formats. output_format = ', output_format)

    if os.path.isdir(output_dir):
        file_utils.rm(output_dir)
    print('concatenating', len(outfiles), 'output files into a list of structures', flush=True)
    outfiles = file_utils.glob(file_utils.fname_from_fpath(final_filename) + '*')
    lines_list = file_utils.concatenate_files(outfiles, final_filename, return_lines=True)
    
    indices_list = file_utils.grep('structure number', final_filename, return_line_nums=True)[1] + [len(lines_list)]
    structs_list = [
                    lines_list[
                                indices_list[i]
                                :
                                indices_list[i + 1]
                                ]
                    for i in range(len(indices_list) - 1)
                    ]
    if num_structures >= len(indices_list):
        selected_structs_idx_list = np.arange(len(structs_list))
    else:
        selected_structs_idx_list = np.random.choice(len(structs_list), num_structures, replace=False)
    
    selected_structs_list = list(np.array(structs_list)[selected_structs_idx_list])
    selected_structs_list = list(map(list, selected_structs_list))
    print('done concatenating', flush=True)
    if output_format == 'json' or output_format == "both":
        print('writing json files', flush=True)
        selected_structs_list = write_json_files(final_filename, output_dir, selected_structs_list)
    else:
        for i, struct_lines in enumerate(selected_structs_list):
            random_str = file_handler.get_random_index()
            struct_id = file_utils.fname_from_fpath(final_filename) + '_' + random_str
            selected_structs_list[i][0] = '#structure number = ' + str(i) + '\n'
            selected_structs_list[i] = selected_structs_list[i][:2] + ['#struct_id = ' + struct_id + '\n'] + selected_structs_list[i][2:]
    
    if output_format != 'json':
        print('writing geometry output file', flush=True)
        # Update the current geometry.out file with struct ids and chronological structure numbers
        for i, struct_lines in enumerate(selected_structs_list):
            if i == 0:
                mode = 'w'
            else:
                mode = 'a'
            file_utils.write_lines_to_file(final_filename, struct_lines, mode=mode)


def clean_up(final_filename, comm_size, output_format):
    for i in range(comm_size):
        os.remove(file_utils.fname_from_fpath(final_filename) + str(i) + '.' + final_filename.split('.')[-1])

    if output_format == 'json':
        os.remove(final_filename)

def pygenarris_structure_generation(inst=None, comm=None, filename=None, num_structures_per_allowed_SG_per_rank=None, Z=None, volume_mean=None, volume_std=None, sr=None, tol=None, max_attempts=None, molecule_path=None):
    # Currently does not support multiple instances running simultaneously. If you want more simultaneous processes, submit more MPI ranks.
    # Might support restart? (i.e. each rank pickup up from where it left off)
    
    if inst is not None:
        sname = 'pygenarris_structure_generation'
        num_structures = inst.get_eval(sname, 'num_structures')
        num_structures_per_allowed_SG_per_rank = inst.get_eval(sname, 'num_structures_per_allowed_SG_per_rank')
        Z = inst.get_eval(sname, 'Z')
        volume_mean = inst.get_eval(sname, 'volume_mean')
        volume_std = inst.get_eval(sname, 'volume_std')
        sr = inst.get_eval(sname, 'sr')
        tol = inst.get_eval(sname, 'tol')
        max_attempts = inst.get_eval(sname, 'max_attempts')
        final_filename = inst.get_with_default(sname, 'geometry_out_filename', 'geometry.out')
        sname_list = [sname, 'relax_single_molecule', 'estimate_unit_cell_volume', 'harris_single_molecule_prep', 'pygenarris_structure_generation', 'structure_generation_batch', 'harris_approximation_batch']
        molecule_path = inst.get_inferred(sname, sname_list, ['molecule_path'] * 7, type_='file')
        output_format = inst.get_with_default(sname, 'output_format', 'json') #options are json, geometry, both
        output_dir = inst.get_with_default(sname, 'output_dir', '.')
    else:
        final_filename = filename
    
    if comm is None:
        filename = final_filename
        set_up(molecule_path)
    else:
        if comm.rank == 0:
            set_up(molecule_path)
        filename = file_utils.fname_from_fpath(final_filename) + str(comm.rank) + '.' + final_filename.split('.')[-1]


    comm.barrier()
    molecule_struct = read(molecule_path)
    mb = MoleculeBonding(molecule_struct)
    cutoff_matrix = mb.get_crystal_cutoff_matrix(Z, vdw_mult=sr)
    cutoff_matrix = np.array(cutoff_matrix, dtype='float32')
    pygenarris.generate_molecular_crystals_with_vdw_cutoff_matrix(filename, cutoff_matrix, num_structures_per_allowed_SG_per_rank, Z, volume_mean, volume_std, tol, max_attempts)
    if comm is not None:
        comm.barrier()
        if comm.rank == 0:
            format_output(output_format, output_dir, final_filename, num_structures)
            clean_up(final_filename, comm.size, output_format)
    
    
