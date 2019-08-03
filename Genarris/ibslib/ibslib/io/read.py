# -*- coding: utf-8 -*-


import os

import ase
from ase.io.formats import all_formats as ase_all_formats

from pymatgen.io.cif import CifParser

from ibslib import Structure 

implemented_file_formats = ["geo", "geometry", "aims", "json", "cif", 
                            "ase"]

def read(struct_path, file_format=''):
    """ 
    Wrapper which will call loading a directory or loading a file depending 
    on the input type. 
    """
    # If file format is used, check if it's an accepted value
    if len(file_format) > 0:
        if file_format not in implemented_file_formats:
            raise Exception("File format {} is not an implemented "
                    .format(file_format) +
                    "format type. Try one of {}."
                    .format(implemented_file_formats))
    
    if os.path.isdir(struct_path):
        return read_dir(struct_path, file_format)
    elif os.path.isfile(struct_path):
        return read_file(struct_path, file_format)
    else:
        raise Exception("Input {} ".format(struct_path) +
                "was not recoginized as a file or a directory")


def read_dir(struct_dir, file_format=''):
    """
    Import any type of structure from directory to the structure class and 
    returns all structure objects as a dictionary, the best python data 
    structure due to hashing for O(1) lookup time.
    """
    file_list = os.listdir(struct_dir)
    struct_dict = {}
    for file_name in file_list:
        file_path = os.path.join(struct_dir,file_name)
        struct = read_file(file_path, file_format)
        if struct == None:
            continue
        struct_dict[struct.struct_id] = struct
            
    return struct_dict


def read_file(file_path, file_format=''):
    """
    Imports a single file to a structure object.
    """
    if len(file_format) > 0:
        if file_format in ["geometry","geo", "aims"]:
            struct = import_geo(file_path)
        elif file_format == "json":
            struct = import_json
        elif file_format == "cif":
            struct = import_cif(file_path)
        elif file_format == "ase":
            struct = import_ase(file_format)
            
    elif '.json' == file_path[-5:]:
            struct = import_json(file_path)
    elif '.cif' == file_path[-4:]:
        struct = import_cif(file_path)
    elif '.in' == file_path[-3:]:
        struct = import_geo(file_path)
    else:
        try: struct = import_ase(file_path)
        except: 
            print("Could not load file {}".format(file_path))
            return None
        
    return struct

    
def import_structures(struct_dir):
    """ Legacy function
    """
    print("WARNING: import_structures is changed to read_dir. "+
          "Please change source code to reflect this.")
    return read_dir(struct_dir)
    

def import_json(file_path):
    struct = Structure()
    struct.build_struct_from_json_path(file_path)
    if struct.struct_id == None:
        file_name = os.path.basename(file_path)
        struct.struct_id = file_name.replace('.json','')
    return struct


def import_geo_dir(struct_dir):
    '''    
    Import all geometry files in a directory to the structure class and 
    returns all the structure objects as a dictionary. Dictionary keys
    are the struct_id of each structure which is set to the filename
    of the geometry file without '.in'. 
    '''
    file_list = os.listdir(struct_dir)
    
    struct_dict = {}
    for file_name in file_list:
        file_path = os.path.join(struct_dir, file_name)
        struct = import_geo(file_path)
        struct_dict[struct.struct_id] = struct
    
    return struct_dict


def import_geo(file_path, struct_id=''):
    """
    Import single geometry file. 
    
    Arguments
    ---------
    file_path: str
        Path to the geometry file to be loaded
    struct_id: str
        Option to specify a struct_id. The default is an empty string which 
        specifies the behavior to use the file name as the struct_id. 
    """
    struct = Structure()
    struct.build_geo_from_atom_file(file_path)
    if len(struct_id) == 0:
        file_name = os.path.basename(file_path)
        struct.struct_id = file_name.replace('.in','')
    else:
        struct.struct_id = struct_id
    return struct


def import_cif(file_path,occupancy_tolerance=100):
    """
    Import cif with pymatgen. This is the current default.
    
    Aruments
    --------
    occupancy_tolerance: int
        Occupancy tolerance for pymatgen.io.cif.CifParser. Files from the 
        CSD can contain multiple atoms which sit on the same site when 
        PBC are taken into account. Pymatgen correctly recognizes these as the
        same atom, however, it has a tolerance for the number of atoms which it
        will assume are the same. This is the occupancy tolerance. I don't 
        believe this value should be limited and thus it's set at 100.
    """
    pycif = CifParser(file_path, occupancy_tolerance=occupancy_tolerance)
    pstruct_list = pycif.get_structures()
    struct = Structure()
    struct.from_pymatgen(pstruct_list[0])
    file_name = os.path.basename(file_path)
    struct.struct_id = file_name.replace('.cif','')
    return struct


def import_cif_ase(file_path):
    """
    Import a single cif file using ase.io
    """
    atoms = ase.io.read(file_path,format='cif')
    struct = Structure()
    struct.from_ase(atoms)
    file_name = os.path.basename(file_path)
    struct.struct_id = file_name.replace('.cif','')
    return struct


def import_ase(file_path):
    """
    Import general file using ase.io
    """
    atoms = ase.io.read(file_path)
    struct = Structure()
    struct.from_ase(atoms)
    file_name = os.path.basename(file_path)
    # Struct ID is filename before decimal 
    struct.struct_id = file_name.split('.')[0]
    return struct




