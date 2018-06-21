"""
If any part of this module is used for a publication please cite:

X. Li, F. Curtis, T. Rose, C. Schober, A. Vazquez-Mayagoitia, K. Reuter,
H. Oberhofer, and N. Marom "Genarris: Random Generation of Molecular Crystal 
Structures and Fast Screening with a Harris Approximation, ",
J. Chem. Phys., DOI: 10.1063/1.5014038; arXiv 1803.02145 (2018)
"""
import copy, math
import numpy as np
import numpy # To allow for explicit user input
from utilities import write_log, misc
from core import structure_handling
from evaluation.evaluation_util import BatchSingleStructureOperation, \
        load_batch_single_structure_operation_keywords
from evaluation.interatomic_distance_evaluation import all_interatomic_distances


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

def rdf_descriptor_by_bin(inst):
    '''
    Returns a list of interatomic distances between pairs of certain species
    '''
    sname = "rdf_descriptor_by_bin" #Module name
    kwargs = load_batch_single_structure_operation_keywords(inst, sname)

    atomic_pairs = inst.get_eval(sname, "atomic_pairs")
    if isinstance(atomic_pairs[0],str):
        # Operation function expecting list of list(s)
        atomic_pairs = [atomic_pairs]

    distance_range = inst.get_eval(sname, "distance_range")
    bins = inst.get_with_default(sname, "bins", 10, eval=True)
    normalized = inst.get_boolean(sname, "normalized")
    property_name = inst.get_with_default(sname, "property_name",
                                          "rdf_descriptor_by_bin")
    output_info_file = inst.get_or_none(sname, "output_info_file")
    
    op_args = (atomic_pairs, distance_range)
    op_kwargs = {"bins"          : bins,
                 "normalized"    : normalized,
                 "property_name" : property_name}

    op = BatchSingleStructureOperation(
            _rdf_descriptor_by_bin,
            name=sname, args=op_args,
            kwargs=op_kwargs, **kwargs)
    result = op.run()

    if output_info_file != None:
        _radial_distribution_function_callback(
                result, output_info_file)

    return result

def _rdf_descriptor_by_bin(struct, atomic_pairs, distance_range,
        bins=10, normalized=False, property_name="rdf_descriptor_by_bin"):
    original_struct = struct
    struct = copy.deepcopy(original_struct)

    # First figures out the cell extension needed
    # Uses napm=1 because of no need to keep atoms
    # from the same molecule together
    structure_handling.cell_modification(struct,
            struct.get_n_atoms(), 1, False)

    structure_handling.cell_lower_triangular(struct, create_duplicate=False)

    a_ext = int(math.ceil(distance_range[1] /
        struct.properties["lattice_vector_a"][0]))
    b_ext = int(math.ceil(distance_range[1] /
        struct.properties["lattice_vector_b"][1]))
    c_ext = int(math.ceil(distance_range[1] /
        struct.properties["lattice_vector_c"][2]))

    structure_extension = [[x,y,z]
                           for x in range (-a_ext,a_ext+1) \
                           for y in range (-b_ext,b_ext+1) \
                           for z in range (-c_ext,c_ext+1)]
    structure_handling.cell_extension(struct,
            extension=structure_extension, create_duplicate=False)

    rdf_descriptor = []
    for pair in atomic_pairs:
        num_of_reference_atoms, distances = \
                _get_relevant_interatomic_distances(
                        struct, original_struct, pair)
        vector, bins_boundaries = np.histogram(distances, bins=bins,
                range=tuple(distance_range))
        vector = list(vector)
        bins_boundaries = list(bins_boundaries)

        # Always normalizes according to number of reference atoms (pair[0])
        vector = [x / num_of_reference_atoms for x in vector]

        if normalized:
            # Normalizes each bin by the volume of the sphere
            for i in range (len(vector)):
                vector[i] = (vector[i] /
                        (4 / 3.0 * math.pi *
                            (bins_boundaries[i+1] ** 3 -
                             bins_boundaries[i] ** 3)))

        rdf_descriptor += vector[:]

    original_struct.set_property(property_name, rdf_descriptor)
    return (original_struct, rdf_descriptor)
       
def rdf_descriptor_by_point(inst):
    '''
    Returns a list of interatomic distances between pairs of certain species

    ''' 
    sname = "rdf_descriptor_by_point" #Module name
    kwargs = load_batch_single_structure_operation_keywords(inst, sname)

    atomic_pairs = inst.get_eval(sname,"atomic_pairs")
    if isinstance(atomic_pairs[0],str):
        # Operation function expecting list of list(s)
        atomic_pairs = [atomic_pairs]

    sampled_distances = eval(inst.get(sname,"sampled_distances"))
    smoothing_parameter = inst.get_with_default(
            sname, "smoothing_parameter", 1, eval=True)
    normalized = inst.get_boolean(sname,"normalized")
    property_name = inst.get_with_default(
            sname, "property_name", "rdf_descriptor_by_point")
    output_info_file = inst.get_or_none(sname, "output_info_file")
    
    op_args = (atomic_pairs, sampled_distances)
    op_kwargs = {"smoothing_parameter" : smoothing_parameter,
                 "normalized"          : normalized,
                 "property_name"       : property_name}

    op = BatchSingleStructureOperation(
            _rdf_descriptor_by_point,
            name=sname, args=op_args,
            kwargs=op_kwargs, **kwargs)
    result = op.run()

    if output_info_file != None:
        _radial_distribution_function_callback(
                result, output_info_file)

    return result

def _rdf_descriptor_by_point(struct, atomic_pairs, sampled_distances,
       smoothing_parameter=1, normalized=False,
       property_name="rdf_descriptor_by_point"):
    original_struct = struct
    struct = copy.deepcopy(original_struct)

    # First figures out the cell extension needed
    # Uses napm=1 because of no need to keep atoms
    # from the same molecule together
    structure_handling.cell_modification(struct,
            struct.get_n_atoms(), 1, False)

    structure_handling.cell_lower_triangular(struct, create_duplicate=False)

    a_ext = int(math.ceil(max(sampled_distances) /
        struct.properties["lattice_vector_a"][0]))
    b_ext = int(math.ceil(max(sampled_distances) /
        struct.properties["lattice_vector_b"][1]))
    c_ext = int(math.ceil(max(sampled_distances) /
        struct.properties["lattice_vector_c"][2]))
    structure_extension = [[x,y,z]
                           for x in range (-a_ext-1,a_ext+2)\
                           for y in range (-b_ext-1,b_ext+2)\
                           for z in range (-c_ext-1,c_ext+2)]
    structure_handling.cell_extension(
            struct,extension=structure_extension,create_duplicate=False)

    if normalized:
        a = smoothing_parameter
        # normalizing factor is the integral of the Gaussian function
        n_factor = [4*math.pi*(0.25*(math.pi/a**3)**0.5+\
        0.5*r**2*(math.pi/a)**0.5+3*r/a*math.exp(-a*r**2)+\
        (0.25*(math.pi/a**3)**0.5+0.5*r**2*(math.pi/a)**0.5)*math.erf(r*a**0.5)) \
        for r in sampled_distances]

    rdf_descriptor = []
    for pair in atomic_pairs:
        num_of_reference_atoms, distances = \
                _get_relevant_interatomic_distances(
                        struct, original_struct, pair)
        g = [sum([math.exp(-smoothing_parameter*(r-r_ij)**2)
                  for r_ij in distances]) / num_of_reference_atoms
            for r in sampled_distances]
                  
        if normalized:
            g = [g[i]/n_factor[i] for i in range (len(sampled_distances))]

        rdf_descriptor += g[:]

    original_struct.set_property(property_name, rdf_descriptor)
    return (original_struct, rdf_descriptor)

def _get_relevant_interatomic_distances(struct, original_struct, atomic_pair):
    reference_atom = atomic_pair[0]
    target_atom_type = atomic_pair[1]
    if isinstance(reference_atom, int):
        a1_range = [reference_atom]
    elif isinstance(reference_atom, list):
        a1_range = reference_atom
    elif reference_atom == "X":
        a1_range = range(struct.get_n_atoms())
    elif isinstance(reference_atom, str):
        a1_range = [i for i in range (original_struct.get_n_atoms())\
        if original_struct.geometry[i]["element"]==reference_atom]
    else:
        raise ValueError("Unsupported format of atomic pairs input; ")

    if not isinstance(target_atom_type, str):
        raise ValueError("Second atom in the atomic pair must be "
                         "a letter indicating a species")

    rs = all_interatomic_distances(
            struct, "X", target_atom_type, a1_range=a1_range)
    return len(a1_range), [case[2] for case in rs]

def _radial_distribution_function_callback(result, output_info_file):
    result.sort(key=lambda x : x[0].struct_id)
    message = ""
    for case in result:
        struct = case[0]
        rdf_descriptor = case[1]
        message += struct.struct_id + " " + " ".join(map(str, rdf_descriptor))
        message += "\n"
    write_log.write_log(output_info_file, message, time_stamp=False)
