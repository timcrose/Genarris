"""
If any part of this module is used for a publication please cite:

X. Li, F. Curtis, T. Rose, C. Schober, A. Vazquez-Mayagoitia, K. Reuter,
H. Oberhofer, and N. Marom "Genarris: Random Generation of Molecular Crystal 
Structures and Fast Screening with a Harris Approximation, ",
J. Chem. Phys., DOI: 10.1063/1.5014038; arXiv 1803.02145 (2018)
"""
import multiprocessing, copy
import numpy as np
from utilities import write_log
from core import structure_handling
from evaluation.evaluation_util import BatchSingleStructureOperation, \
        load_batch_single_structure_operation_keywords


__author__ = "Xiayue Li, Timothy Rose, Christoph Schober, and Farren Curtis"
__copyright__ = "Copyright 2018, Carnegie Mellon University and "+\
                "Fritz-Haber-Institut der Max-Planck-Gessellschaft"
__credits__ = ["Xiayue Li", "Luca Ghiringhelli", "Farren Curtis", "Tim Rose",
               "Christoph Schober", "Alvaro Vazquez-Mayagoita",
               "Karsten Reuter", "Harald Oberhofer", "Noa Marom"]
__license__ = "BSD-3"
__version__ = "1.0"
__maintainer__ = "Timothy Rose"
__email__ = "trose@andrew.cmu.edu"
__url__ = "http://www.noamarom.com"

def main(inst):
    '''
    Returns a list of interatomic distances between pairs of certain species

    '''
    sname = "interatomic_distance_evaluation" #Module name
    kwargs = load_batch_single_structure_operation_keywords(inst, sname)

    napm = inst.get_eval(sname, "NAPM")
    atomic_pairs = inst.get_eval(sname, "atomic_pairs")
    default_extension = [[x,y,z] for z in range (-1,2)
                                 for y in range (-1,2) 
                                 for x in range (-1,2)
                                 if x!=0 or y!=0 or z!=0]
    structure_extension = inst.get_with_default(sname, "structure_extension", \
                        default_extension, eval=True)
    allow_all_pairs = inst.get_boolean(sname, "allow_all_pairs")
    block_same_pair = inst.get_boolean(sname, "block_same_pair")
    property_name = inst.get_with_default(sname,
            "property_name", "interatomic_distances")
    op_args = (napm, atomic_pairs)
    op_kwargs = {"structure_extension" : structure_extension,
                 "allow_all_pairs"     : allow_all_pairs,
                 "block_same_pair"     : block_same_pair,
                 "property_name"       : property_name}

    op = BatchSingleStructureOperation(
            _interatomic_distance_evaluation,
            name=sname, args=op_args,
            kwargs=op_kwargs, **kwargs)
    result = op.run()

    if inst.has_option(sname,"output_info_file"):
        output_file = inst.get(sname,"output_info_file")
        _interatomic_distance_evaluation_callback(
                result, output_file)

def _interatomic_distance_evaluation(struct, napm, atomic_pairs,
        structure_extension=[[x, y, z] for z in range(-1, 2)
                                       for y in range(-1, 2)
                                       for x in range(-1, 2)
                                       if x!=0 or y!=0 or z!=0],
        allow_all_pairs=False, block_same_pair=False,
        property_name="interatomic_distance"):

    # Makes deepcopy to avoid numpy's resizing error
    struct = copy.deepcopy(struct)
    structure_handling.cell_modification(
            struct, struct.get_n_atoms()//napm,
            napm, create_duplicate=False)
    structure_handling.cell_extension(struct, extension=structure_extension,
            create_duplicate=False)
    
    result = []

    for pair in atomic_pairs:
        l = specific_interatomic_distances(struct, pair[0], pair[1], napm, \
                allow_all_pairs=allow_all_pairs,
                block_same_pair=block_same_pair)

        if len(l) < pair[2]:
            # Not enough interatomic distance returned to select
            # the pair[2]th closest
            result.append([-1,-1,-1])
            continue
        
        for i in range(pair[2]): #Now select the pair[2]th closest
            r = _min_distance(l)
            l.remove(r)
        result.append(r)

    return (struct, result)

def _interatomic_distance_evaluation_callback(result, output_file):
    result.sort(key=lambda x : x[0].struct_id)
    message = ""
    for case in result:
        struct = case[0]
        case_result = case[1]
        message += struct.struct_id
        for r in case_result:
            message += " " + " ".join(map(str, r))
        message += "\n"
    write_log.write_log(output_file, message, time_stamp=False)

# TODO: Refactor all interatomic distance calculation into an
# InteratomicDistanceCalculation class
def specific_interatomic_distances (struct,
        a1, a2, napm, allow_all_pairs=False, block_same_pair=False):
    '''
    Returns a list of interatomic distances
    
    If a1 is a string, then it specifies the species of the first atom
    If a1 is a integer, then it specifies the index of the first atom 
    (This overrules the setting of allow_all_pairs)
    
    if allow_all_pairs set to True, then a1 does not need to come from the first molecule
    if block_same_pair is set to True, then only the closest a1-a2 contact will be recorded for each molecule pair

    '''
    # TODO: Write unit test specifically testing out this module
    if not len(struct.geometry) % napm == 0:
        raise ValueError("napm does not divide structure geometry length")

    nmpc = len(struct.geometry)/napm
    napc = struct.get_n_atoms()

    if isinstance(a1, int) and isinstance(a2, int):
        if a1 // napm == a2 // napm:
            # Does not return result if from the same molecule
            return []
        return [[a1 , a2, np.linalg.norm([\
            struct.geometry[a1]["x"]-struct.geometry[a2]["x"],\
            struct.geometry[a1]["y"]-struct.geometry[a2]["y"],\
            struct.geometry[a1]["z"]-struct.geometry[a2]["z"]])]]

    elif isinstance(a1, int) and isinstance(a2, str):
        a1_element = struct.geometry[a1]["element"]
        if not block_same_pair:
            return all_interatomic_distances(
                    struct, a1_element, a2,
                    a1_range=[a1], a2_range=range(0, napc),
                    napm=napm)
        else:
            result = []
            # Iterate through all the molecules except the one a1 is in
            for i in range(0, a1//napm) + range(a1//napm+1, nmpc):
                d = all_interatomic_distances(
                        struct, a1_element, a2,
                        a1_range=[a1], a2_range=range(i*napm, (i+1)*napm),
                        napm=napm)
                result.append(_min_distance(d))
            return [x for x in result if x != None]

    elif isinstance(a1, str) and isinstance(a2, int):
        a2_element = struct.geometry[a2]["element"]
        if not block_same_pair:
            if not allow_all_pairs:
                if a2 < napm:
                    # Immediately returns because atoms must be from different
                    # molecules. However, if a1 must be selected from the first
                    # molecule since allow_all_pairs=False, and a2 is from the
                    # first molecule as well, then no result can be returned
                    return []
                return all_interatomic_distances(
                        struct, a1, a2_element,
                        a1_range=range(0, napm), a2_range=[a2],
                        napm=napm)
            else:
                # Avoids a1_range from including any atoms of a2 molecule
                a1_range = range(0, napm * (a2 // napm)) + \
                        range(napm * (a2 // napm + 1), napc)
                return all_interatomic_distances(
                        struct, a1, a2_element,
                        a1_range=a1_range, a2_range=[a2],
                        napm=napm)
        else:
            if not allow_all_pairs:
                if a2 < napm:
                    # Immediately returns because atoms must be from different
                    # molecules. However, if a1 must be selected from the first
                    # molecule since allow_all_pairs=False, and a2 is from the
                    # first molecule as well, then no result can be returned    
                    return []
                else:
                    d = all_interatomic_distances(
                            struct, a1, a2_element,
                            a1_range=range(0, napm), a2_range=[a2],
                            napm=napm)
                    return [_min_distance(d)]
            else:
                result = []
                # Iterates through all molecules in the struct
                # avoiding the molecule of a2
                for i in range (0, a2 // napm) + range(a2 // napm + 1, nmpc):
                    d = all_interatomic_distances(
                            struct, a1, a2_element,
                            a1_range=range(i*napm, (i+1)*napm),
                            a2_range=[a2], napm=napm)
                    result.append(_min_distance(d))
                return [x for x in result if x != None]

    elif isinstance(a1, str) and isinstance(a2, str):
        if not block_same_pair and not allow_all_pairs:
            return all_interatomic_distances(
                    struct, a1, a2, a1_range=range(0, napm),
                    a2_range=range(napm, napc), napm=napm)

        elif not block_same_pair and allow_all_pairs:
            return all_interatomic_distances(
                    struct, a1, a2, a1_range=range(0, napc),
                    a2_range=range(0, napc), napm=napm)

        elif block_same_pair and not allow_all_pairs:
            result = []
            for i in range (1, nmpc):
                d = all_interatomic_distances(struct, a1, a2,
                        a1_range=range(0, napm),
                        a2_range=range(i*napm,(i+1)*napm))
                result.append(_min_distance(d))
            return [x for x in result if x != None]
        # From this point on, both block_same_pair and
        # allow_all_pairs are true
        if a1 == a2:
            result = []
            for i in range (0, nmpc-1):
                for j in range (i, nmpc):
                    d = all_interatomic_distances(
                            struct, a1, a2,
                            a1_range=range(i*napm, (i+1)*napm),
                            a2_range=range(j*napm, (j+1)*napm),
                            napm=napm)
                    result.append(_min_distance(d))
            return [x for x in result if x != None]
        else:
            # block_same_pair does not block the simultaneous selection
            # of (X from M1, Y from M2) and (X from M2, Y from M1) if X != Y
            result = []
            for i in range (0, nmpc):
                for j in range(0, i) + range(i+1, nmpc):
                    d = all_interatomic_distances(
                            struct, a1, a2,
                            a1_range=range(i*napm,(i+1)*napm),
                            a2_range=range(j*napm,(j+1)*napm),
                            napm=napm)
                    result.append(_min_distance(d))
            return [x for x in result if x!=None]
    else:
        raise ValueError("Invalid input of a1 and a2")

def all_interatomic_distances (
        struct, a1, a2, a1_range=None, a2_range=None, napm=1):
    '''
    Returns a list of interatomic distances 
    a1, a2 specifies the species of the two atoms
    Returns [[a1_index_1, a2_index_1, distance_1],
             [a1_index_2, a2_index_2, distance_2] . . .]

    Requires a1_index/napm to be different from a2_index/napm
    '''
    result = []
    if a1_range == None:
        a1_range = range(len(struct.geometry))
    if a2_range == None:
        a2_range = range(len(struct.geometry))

    if a1 == a2: #Need to avoid double counting
        for i in a1_range:
            for j in a2_range:
                if (i // napm != j // napm) and \
                not ((i in a2_range) and (j in a1_range) and (i>j)) \
                and (a1=="X" or struct.geometry[i]["element"] == a1) \
                and (a2=="X" or struct.geometry[j]["element"] == a2):
                    result.append([i,j,np.linalg.norm([\
                    struct.geometry[i]["x"]-struct.geometry[j]["x"],\
                    struct.geometry[i]["y"]-struct.geometry[j]["y"],\
                    struct.geometry[i]["z"]-struct.geometry[j]["z"]])])
    
    if a1 != a2:
        for i in a1_range:
            for j in a2_range:
                if i // napm != j // napm \
                and (a1=="X" or struct.geometry[i]["element"] == a1) \
                and (a2=="X" or struct.geometry[j]["element"] == a2):
                    result.append([i,j,np.linalg.norm([\
                    struct.geometry[i]["x"]-struct.geometry[j]["x"],\
                    struct.geometry[i]["y"]-struct.geometry[j]["y"],\
                    struct.geometry[i]["z"]-struct.geometry[j]["z"]])])
    return result

def _min_distance(distance_list):
    '''
    Returns the pair of atoms with the largest distance from the list.
    '''
    if len(distance_list)==0:
        return None
    m = 0
    for i in range (1, len(distance_list)):
        if distance_list[i][2] < distance_list[m][2]:
            m = i
    return distance_list[m]

