# -*- coding: utf-8 -*-

import math # need for error function
import copy
from copy import deepcopy
import numpy as np
import time

from ibslib.structures import Structure
from ibslib.structures.structure_handling import *
from ibslib.descriptor.distance_matrix import calc_euclidean_dist_vectorized,calc_atom_dist
from ibslib.structures.return_geo_array import get_geo_array
from ibslib.analysis.pltlib import standard_formatting

lat_interp = {0:'lattice_vector_a',1:'lattice_vector_b',2:'lattice_vector_c'}
molar_mass = {
        'H': 1.01,
        'N': 14.01,
        'C': 12.01,
        'O': 16.00,
        'Br': 79.90
        }

def eval_dict(struct_dict, parallel=1, **kwargs):
    '''
    Purpose:
        Evaluates RDF of an entire struct_dict
    '''
    total = len(struct_dict)
    for i,struct_id in enumerate(struct_dict):
        print('############# {} Structures Left for {} #############'
              .format(total-i,kwargs['atomic_pairs_list']))
        struct = struct_dict[struct_id]
        # Old version is slightly incorrect but much faster
        compute_RDF_vector_neighbor(struct, **kwargs)

def eval_list(struct_list,**kwargs):
    for struct in struct_list:
        compute_RDF_vector_neighbor(struct, **kwargs)

def compute_RDF_vector(original_struct, atomic_pairs_list=['C','C'], 
                       atomic_distance_range=[1,8],
                       atomic_distance_increment=1, 
                       smoothing_parameter=1):
    
    ###########################################################################
    start_time = time.time()
    ###########################################################################
    struct = deepcopy(original_struct)
    
    atomic_distance_increment = float(atomic_distance_increment)
    dist = np.arange(float(atomic_distance_range[0]), 
                     # Add atomic distance increment to make dist inclusive 
                     #  of the final distance value
                     float(atomic_distance_range[1])+atomic_distance_increment,
                     float(atomic_distance_increment))
    
    unit_cell_density = len(original_struct.geometry)/ \
                        (original_struct.get_property('cell_vol'))
    inv_density = 1/unit_cell_density
    #Computer Normalization Factors
    a = smoothing_parameter
    n_factor = []
    
    # Will perform shell volume calculations, for example for [0,10] by 1:
    # [[1-0], [2-1], [3-2], [4-3], ..., [20-19]]
    bin_value_list = []
    for i,value in enumerate(dist[1:]):
        shell_vol = (4.0/3.0)*np.pi*(value**3 - (value-atomic_distance_increment)**3)
        inv_shell_vol = 1/shell_vol
        n_factor.append(inv_shell_vol*inv_density)
        bin_value_list.append(value - atomic_distance_increment/2)
    
    A = struct.properties["lattice_vector_a"]
    B = struct.properties["lattice_vector_b"]
    C = struct.properties["lattice_vector_c"]
 
    #First figure out the cell extension needed
    cell_modification(struct,struct.get_n_atoms(), False)
    cell_lower_triangular(struct, False)
    a_ext = int(math.ceil(max(dist)/A[0]))
    b_ext = int(math.ceil(max(dist)/B[1]))
    c_ext = int(math.ceil(max(dist)/C[2]))
    structure_extension = [[x,y,z] for x in range (-a_ext-1,a_ext+2)\
                           for y in range (-b_ext-1,b_ext+2)\
                           for z in range (-c_ext-1,c_ext+2)]
    
    # THIS WASN'T HERE PREVIOUSLY
    structure_extension.remove([0,0,0])
    
    ex_struct = cell_extension(struct, extension=structure_extension, 
                               create_duplicate=True)

    # Compute interatomic distances and build g vector
    vector_all = [struct.struct_id]
    atomic_pairs = [atomic_pairs_list[i:i+2] for i in range(0, len(atomic_pairs_list), 2)]
                
    ###########################################################################
    before_for_loop = time.time()
    print('Time taken before for loop: {}'.format(before_for_loop-start_time))
    ###########################################################################
    napm = int(len(original_struct.geometry)/nmpc)
    
    for pair in atomic_pairs:
        #print pair
        ref_atom_type = pair[0]
        target_atom_type = pair[1]

        a1_range = [i for i in range(original_struct.get_n_atoms()) 
                        if original_struct.geometry[i]["element"] == ref_atom_type]
        a2_range = [i for i in range(original_struct.get_n_atoms()) 
                        if original_struct.geometry[i]["element"] == target_atom_type]
        
        interatomic_start = time.time()
        
        rs_manny = all_interatomic_distances_manny(ex_struct, napm,
                                       ref_atom_type, target_atom_type, 
                                       a1_range)
        distances = rs_manny
        
        interatomic_end = time.time()
        print('Time for interatomic distances: {}'
              .format(interatomic_end - interatomic_start))
        
        print('MANNY len: ',len(distances))
        print('MANNY min: ',np.min(distances))
        
        #######################################################################
        g_start_time = time.time()
        #######################################################################
        g = []
        
        counting_array = [0 for x in range(len(n_factor))]
        
        # Binning distance values
        for r in distances:
            bin_num = int(r/atomic_distance_increment)
            if bin_num > len(n_factor)-1:
                continue
            counting_array[bin_num] += 1
            bin_value = bin_value_list[bin_num]
#            print(bin_num,bin_value,r)

            
        #######################################################################
        first_g_end_time = time.time()
        print('first g compute: {}'.format(first_g_end_time - g_start_time))
        print('len g: {}'.format(len(g)))
        #######################################################################
        
        g = np.divide(counting_array,np.sum(counting_array))
        g = np.multiply(g,n_factor)
        g = list(g)
        vector_all += g[:]

        #######################################################################
        g_end_time = time.time()
        print('g compute final steps: {}'.format(g_end_time - first_g_end_time))
        #######################################################################
        
    RDF_vec = vector_all[1:]
    original_struct.set_property("RDF_vector", RDF_vec)
    
    ###########################################################################
    after_for_loop = time.time()
    print('Time from for loop to the end: {}'.format(after_for_loop-
                                                     before_for_loop))
    print('Total time: {}'.format(after_for_loop-start_time))
    print('RDF total: {}'.format(np.sum(RDF_vec)))
    ###########################################################################
    
    return original_struct

def compute_RDF_vector_newest(original_struct, atomic_pairs_list=['C','C'], 
                       atomic_distance_range=[1,8],
                       atomic_distance_increment=1, 
                       smoothing_parameter=1,
                       descriptor_key='rdf'):
    
    ###########################################################################
    start_time = time.time()
    ###########################################################################
    struct = deepcopy(original_struct)
    
    
    atomic_distance_increment = float(atomic_distance_increment)
    dist = np.arange(float(atomic_distance_range[0]), 
                     # Add atomic distance increment to make dist inclusive 
                     #  of the final distance value
                     float(atomic_distance_range[1])+atomic_distance_increment,
                     float(atomic_distance_increment))
    
    unit_cell_density = len(original_struct.geometry)/ \
                        (original_struct.get_property('cell_vol'))
    inv_density = 1/unit_cell_density
    
    #Computer Normalization Factors
    a = smoothing_parameter
    n_factor = []
    
    # Will perform shell volume calculations, for example for [0,10] by 1:
    # [[1-0], [2-1], [3-2], [4-3], ..., [20-19]]
    bin_value_list = [0]
    for i,value in enumerate(dist[1:]):
        shell_vol = (4.0/3.0)*np.pi*(value**3 - (value-atomic_distance_increment)**3)
        inv_shell_vol = 1/shell_vol
        n_factor.append(inv_shell_vol*inv_density)
        bin_value_list.append(value - atomic_distance_increment/2)
    
    A = struct.properties["lattice_vector_a"]
    B = struct.properties["lattice_vector_b"]
    C = struct.properties["lattice_vector_c"]
 
    #First figure out the cell extension needed
    cell_modification(struct,struct.get_n_atoms(), False)
    cell_lower_triangular(struct, False)
    a_ext = int(math.ceil(max(dist)/A[0]))
    b_ext = int(math.ceil(max(dist)/B[1]))
    c_ext = int(math.ceil(max(dist)/C[2]))
    structure_extension = [[x,y,z] for x in range (-a_ext-1,a_ext+2)\
                           for y in range (-b_ext-1,b_ext+2)\
                           for z in range (-c_ext-1,c_ext+2)]
    
    # THIS WASN'T HERE PREVIOUSLY
    structure_extension.remove([0,0,0])
    
    ex_struct = cell_extension(struct, extension=structure_extension, 
                               create_duplicate=True)

    # Compute interatomic distances and build g vector
    atomic_pairs = [atomic_pairs_list[i:i+2] for i in range(0, len(atomic_pairs_list), 2)]
                
    ###########################################################################
    before_for_loop = time.time()
    print('Time taken before for loop: {}'.format(before_for_loop-start_time))
    ###########################################################################
    
    g = [0 for x in range(len(dist))]
    for pair in atomic_pairs:

        ref_atom_type = pair[0]
        target_atom_type = pair[1]

        interatomic_start = time.time()
        
        distances = all_interatomic_distances_manny_newest(ex_struct,
                                       ref_atom_type, target_atom_type)
        
        interatomic_end = time.time()
        print('Time for interatomic distances newest: {}'
              .format(interatomic_end - interatomic_start))
        
        print('MANNY newest num distances: ',len(distances))
        print('MANNY newest min: ',np.min(distances))
        print('MANNY newest max: ',np.max(distances))
        
        #######################################################################
        g_start_time = time.time()
        #######################################################################
        
        index = np.where(distances < atomic_distance_range[1])
        distances = distances[index]
        
        # Binning distance values
        hist,bin_edges = np.histogram(distances, bins=dist, density=False,
                                      range=(0,atomic_distance_range[1]))
        total = np.sum(hist)
        # Solve for the number of particles using biniomial coefficient
        num_particles = 0.5*(np.sqrt(8*total)+1)
        print(total,num_particles)
        
        for i,bin_value in enumerate(hist):
            norm = bin_value / num_particles
            norm = norm / ((4/3)*np.pi*(bin_edges[i]**3-bin_edges[i-1]**3))
            norm = norm*inv_density
            g[i] += norm
        
        #######################################################################
        first_g_end_time = time.time()
        print('first g compute: {}'.format(first_g_end_time - g_start_time))
        print('len g: {}'.format(len(g)))
        #######################################################################
        
    original_struct.set_property(descriptor_key, g)
    
    ###########################################################################
    after_for_loop = time.time()
    print('Time from for loop to the end: {}'.format(after_for_loop-
                                                     before_for_loop))
    print('Newest total time: {}'.format(after_for_loop-start_time))
    ###########################################################################
    
    return original_struct

def compute_RDF_vector_neighbor(original_struct, atomic_pairs_list=['C','C'], 
                       atomic_distance_range=[1,8],
                       atomic_distance_increment=1, 
                       smoothing_parameter=1,
                       descriptor_key='rdf'):
    '''
    - This version does not compute all interatomic distances but builds a 
        a neighborhood within Rc of the relevant atoms in the unit cell. 
    '''
    
    ###########################################################################
    start_time = time.time()
    ###########################################################################
    struct = deepcopy(original_struct)
    
    N = len(np.where(struct.geometry['element'] == atomic_pairs_list[0])[0])
    print('N: ', N)
    
    atomic_distance_increment = float(atomic_distance_increment)
    dist = np.arange(float(atomic_distance_range[0]), 
                     # Add atomic distance increment to make dist inclusive 
                     #  of the final distance value
                     float(atomic_distance_range[1])+atomic_distance_increment,
                     float(atomic_distance_increment))
    unit_cell_density = original_struct.geometry.shape[0]/ \
                        (original_struct.get_property('cell_vol'))
                        
    inv_density = 1/unit_cell_density

    # Compute interatomic distances and build g vector
    atomic_pairs = [atomic_pairs_list[i:i+2] for i in 
                    range(0, len(atomic_pairs_list), 2)]
                
    ###########################################################################
    before_for_loop = time.time()
    print('Time taken before for loop: {}'.format(before_for_loop-start_time))
    ###########################################################################
    
    g = [0 for x in range(len(dist)-1)]
    for pair in atomic_pairs:
        a1 = pair[0]
        a2 = pair[1]
        # Cutoff radius
        Rc = atomic_distance_range[1]

        interatomic_start = time.time()
        
        distances = all_interatomic_distances_neighbor(struct,a1,a2,Rc)
        
        #######################################################################
        interatomic_end = time.time()
        print('Time for interatomic distances newest: {}'
              .format(interatomic_end - interatomic_start))
        
        print('MANNY neighbor num distances: ',len(distances))
        print('MANNY neighbor min: ',np.min(distances))
        print('MANNY neighbor max: ',np.max(distances))
        g_start_time = time.time()
        #######################################################################
        
        # Binning distance values
        hist,bin_edges = np.histogram(distances, bins=len(dist)-1, density=False,
                                      range=(0,atomic_distance_range[1]))
        
        for i,bin_value in enumerate(hist):
            norm = bin_value / N #/ num_particles
            norm = norm / ((4/3)*np.pi*(bin_edges[i]**3-bin_edges[i-1]**3))
            norm = norm*inv_density
            g[i] += norm
        
        #######################################################################
        first_g_end_time = time.time()
        print('first g compute: {}'.format(first_g_end_time - g_start_time))
        print('len g: {}'.format(len(g)))
        #######################################################################
        
    original_struct.set_property(descriptor_key, g)
    
    ###########################################################################
    after_for_loop = time.time()
    print('Neighbor total time: {}'.format(after_for_loop-start_time))
    ###########################################################################
    
    return original_struct

def all_interatomic_distances_neighbor(struct, a1, a2, Rc, nmpc=1):
    '''
    1. Makes a matrix of only the relevant atoms
    2. Computes neighborhood of relevant atom pairs
    3. Distance matrix is set up such that parsing it is trivially easy. This 
         is usually the rate limiting step so improvements here will make a 
         big difference.
    Arguments:
        struct: Structure object
        a1: str of atom element 1
        a2: str of atom element 2
        nmpc: If nmpc is greater than 1, then only intermolecular distances 
                will be computed. 
    '''
    geo = struct.geometry
    
    element_list = geo['element']
    
    a1_index = np.where(element_list == a1)
    
    # Cell extension to compute neighborhood with PBC
    A = struct.properties["lattice_vector_a"]
    B = struct.properties["lattice_vector_b"]
    C = struct.properties["lattice_vector_c"]
 
    #First figure out the cell extension needed
    cell_modification(struct,struct.get_n_atoms(), False)
    cell_lower_triangular(struct, False)
    a_ext = int(math.ceil(Rc/A[0]))
    b_ext = int(math.ceil(Rc/B[1]))
    c_ext = int(math.ceil(Rc/C[2]))
    structure_extension = [[x,y,z] for x in range (-a_ext-1,a_ext+2)\
                           for y in range (-b_ext-1,b_ext+2)\
                           for z in range (-c_ext-1,c_ext+2)]
    try: structure_extension.remove([0,0,0])
    except: pass
    
    ex_struct = cell_extension(struct, extension=structure_extension, 
                               create_duplicate=True)
    
    a2_index = np.where(ex_struct.geometry['element'] == a2)
    geo_array = get_geo_array(ex_struct)
    
    # a1_array from unit cell
    a1_array = geo_array[a1_index[0],:]
    # a2_array from supercell
    a2_array = geo_array[a2_index[0],:]
    
    important_array = np.concatenate([a1_array,a2_array])
    num_a1 = a1_array.shape[0]
    
    dist_matrix = calc_atom_dist(important_array)
    important_dist = dist_matrix[0:num_a1,num_a1:]
    
    # Get upper triangle indices
    if a1 == a2:
        # If the elements are the same then the diagonal of the important_dist 
        #   matrix will be zero. These are not relevant distances if the elements
        #   are the same. If the elements are different, then these are relevant
        #   distances and must be none zero.
        upper = np.triu_indices(important_dist.shape[0], 1, 
                                important_dist.shape[1])
        final_dist = important_dist[upper]
    else:
        # if a1 != a2 then the entire matrix is relevant
        final_dist = important_dist
    
    # Keep only distances within the cutoff radius
    index = np.where(final_dist < Rc)
    final_dist = final_dist[index]
    
    return final_dist

###############################################################################
# Use this in newest rdf evaluation to smooth the function similar to the old #
#  version                                                                    #
###############################################################################
#        dist = np.arange(float(atomic_distance_range[0]), 
#                     float(atomic_distance_range[1]),
#                     float(atomic_distance_increment))
#        a = smoothing_parameter
#        n_factor = [4 * np.pi * (0.25 * (np.pi/a**3)**0.5 + \
#            0.5 *r**2 * (np.pi/a)**0.5 + 3*r/a * math.exp(-a*r**2) + \
#            (0.25 * (np.pi/a**3)**0.5 + 0.5*r**2*(np.pi/a)**0.5)*\
#            math.erf(r*a**0.5)) \
#            for r in dist]
#        
#        g = []    
#        for r in dist:       
#            g.append(np.sum(np.exp(
#                            -smoothing_parameter*
#                            np.square(r-distances))))
#        g = np.divide(g,n_factor)
#        g = list(g)
#        vector_all += g[:]
###############################################################################
    
def compute_RDF_vector_old(original_struct, atomic_pairs_list=['C','C'], 
                       atomic_distance_range=[1,8],
                       atomic_distance_increment=1, 
                       smoothing_parameter=1,
                       descriptor_key='rdf'):
    
    old_start = time.time()

    dist = np.arange(float(atomic_distance_range[0]), 
                     float(atomic_distance_range[1]),
                     float(atomic_distance_increment))
    #Computer Normalization Factors
    a = smoothing_parameter
    n_factor = [4 * np.pi * (0.25 * (np.pi/a**3)**0.5 + \
                0.5 *r**2 * (np.pi/a)**0.5 + 3*r/a * math.exp(-a*r**2) + \
                (0.25 * (np.pi/a**3)**0.5 + 0.5*r**2*(np.pi/a)**0.5)*\
                math.erf(r*a**0.5)) \
                for r in dist]

    A = struct.properties["lattice_vector_a"]
    B = struct.properties["lattice_vector_b"]
    C = struct.properties["lattice_vector_c"]
 
    #First figure out the cell extension needed
    cell_modification(struct,struct.get_n_atoms(), False)
    cell_lower_triangular(struct, False)
    a_ext = int(math.ceil(max(dist)/A[0]))
    b_ext = int(math.ceil(max(dist)/B[1]))
    c_ext = int(math.ceil(max(dist)/C[2]))
    structure_extension = [[x,y,z] for x in range (-a_ext-1,a_ext+2)\
                           for y in range (-b_ext-1,b_ext+2)\
                           for z in range (-c_ext-1,c_ext+2)]
    # THIS WASN'T HERE PREVIOUSLY
    structure_extension.remove([0,0,0])
    ex_struct = cell_extension(struct, extension=structure_extension, create_duplicate=True)

    # Compute interatomic distances and build g vector
    vector_all = [struct.struct_id]
    atomic_pairs = [atomic_pairs_list[i:i+2] for i in range(0, len(atomic_pairs_list), 2)]
    for pair in atomic_pairs:
        #print pair
        ref_atom_type = pair[0]
        target_atom_type = pair[1]
        
        # This only considers the RDF between the orignial unit cell and 
        a1_range = [i for i in range(original_struct.get_n_atoms()) 
                        if original_struct.geometry[i]["element"] == ref_atom_type]
        a2_range = [i for i in range(original_struct.get_n_atoms()) 
                        if original_struct.geometry[i]["element"] == target_atom_type]
        rs = all_interatomic_distances(ex_struct, ref_atom_type, target_atom_type, a1_range)
        distances = [i[2] for i in rs]
        print('Old num distances: {}'.format(len(distances)))
        print('Old min distance: {}'.format(min(distances)))
        print('Old max distance: {}'.format(max(distances)))
        g = [sum([math.exp(-smoothing_parameter*(r-r_ij)**2) 
            for r_ij in distances])/len(a1_range) for r in dist]
        g = [g[i]/n_factor[i] for i in range (len(dist))]
        vector_all += g[:]
    RDF_vec = vector_all[1:]
    struct.set_property(descriptor_key, RDF_vec)
    
    old_end = time.time()
    print('Total old time: {}'.format(old_end - old_start))
    
    return struct

def all_interatomic_distances (struct, a1, a2, a1_range=None, a2_range=None):
    '''
    Returns a list of interatomic distances 
    a1, a2 specifies the species of the two atoms
    Returns [[a1_index_1, a2_index_1, distance_1], [a1_index_2, a2_index_2, distance_2] . . .]
    Requires a1_index/napm to be different from a2_index/napm
    '''
    result = []
    geo = struct.geometry
#    napm=int(len(geo)/nmpc)
    napm = 1
    print('OLD NAPM: {}'.format(napm))
    if a1_range == None:
        a1_range = list(range(len(geo)))
    if a2_range == None:
        a2_range = list(range(len(geo)))
    if a1 == a2: #Need to avoid double counting
        for i in a1_range:
            for j in a2_range:
                if (i//napm!=j//napm) and not ((i in a2_range) and (j in a1_range) and (i>j))\
                and (a1=="X" or geo[i]["element"] == a1) \
                and (a2=="X" or geo[j]["element"] == a2):
                    result.append([i,j,np.linalg.norm([\
                    geo[i]["x"] - geo[j]["x"],\
                    geo[i]["y"] - geo[j]["y"],\
                    geo[i]["z"] - geo[j]["z"]])])
    
#    print('a1_range: ',a1_range)
#    print('a2_range: ', a2_range)
    
    if a1 != a2:
        for i in a1_range:
            for j in a2_range:
                if i//napm!=j//napm\
                and (a1=="X" or geo[i]["element"] == a1) \
                and (a2=="X" or geo[j]["element"] == a2):
                    result.append([i,j,np.linalg.norm([\
                    geo[i]["x"] - geo[j]["x"],\
                    geo[i]["y"] - geo[j]["y"],\
                    geo[i]["z"] - geo[j]["z"]])])
    return result

def all_interatomic_distances_manny_newest(struct, a1, a2):
    '''
    1. Makes a matrix of only the relevant atoms
    2. Computes all the distances
    3. Distance matrix is set up such that parsing it is trivially easy. This 
         is usually the rate limiting step so improvements here will make a 
         big difference.s
    '''
    geo = struct.geometry
    
    element_list = [x[3] for x in geo]   
    element_list = np.array(element_list)
    
    a1_index = np.where(element_list == a1)
    a2_index = np.where(element_list == a2)
    
    geo_array = get_geo_array(struct)
    
    a1_array = geo_array[a1_index[0],:]
    a2_array = geo_array[a2_index[0],:]
    important_array = np.concatenate([a1_array,a2_array])
    num_a1 = a1_array.shape[0]
    num_a2 = a2_array.shape[0]
    
    dist_matrix = calc_atom_dist(important_array)
    important_dist = dist_matrix[0:num_a1,num_a1:]
    
    # Get upper triangle indices
    if a1 == a2:
        # If the elements are the same then the diagonal of the important_dist 
        #   matrix will be zero. These are not relevant distances if the elements
        #   are the same. If the elements are different, then these are relevant
        #   distances and must be none zero.
        upper = np.triu_indices(important_dist.shape[0], 1, 
                                important_dist.shape[1])
        final_dist = important_dist[upper]
    else:
        # if a1 != a2 then the entire matrix is relevant
        final_dist = important_dist
        
    
    
    return final_dist

def all_interatomic_distances_manny(struct, napm, 
                               a1, a2, a1_range=None, a2_range=None):
    '''
    Returns a list of interatomic distances 
    a1, a2 specifies the species of the two atoms
    Returns [[a1_index_1, a2_index_1, distance_1], [a1_index_2, a2_index_2, distance_2] . . .]
    Requires a1_index/napm to be different from a2_index/napm
    '''
    ###########################################################################
    start_time = time.time()
    ###########################################################################
    
    geo = struct.geometry
    num_molecules = int(len(geo)/napm)
    print('Number molecules: ',len(geo)/napm)
    
    element_list = [x[3] for x in geo[0:napm]]   
    element_list = np.array(element_list)
    
    a1_index = np.where(element_list == a1)
    a2_index = np.where(element_list == a2)            

    a1_index = [x for x in a1_index[0]]
    a2_index = [x for x in a2_index[0]]
        
   
    # Building index list that will be uesd to parse the distance matrix
    a1_index_list = []
    a2_index_list = []
    a1_factor = 0
    a2_factor = 0
    # Count from first molecule to the second to last molecule for a1 list
    for num in range(0,num_molecules-1):
        a1_factor = num*napm
        a2_factor = (num+1)*napm
        a1_index_list.append([x+a1_factor for x in a1_index])
        # Adding extra num*napm so that a2_list is built from separate 
        #  molecules thanm the a1 value such that it refers to 
        #  intermolecular information
        a2_index_list.append([x+a2_factor for x in a2_index])
    
    final_dist_index = []
    for i,a1_list in enumerate(a1_index_list):
        for a1 in a1_list:
           for a2_list in a2_index_list[i+1:]:
               for a2 in a2_list:
                   final_dist_index.append((a1,a2))
                      
    atom_coords = []
    for atom_info in geo:
        atom_coords.append([atom_info[0],atom_info[1], atom_info[2]])
    
    atom_coords = np.array(atom_coords)
    all_distances = calc_euclidean_dist_vectorized(atom_coords)
    
    final_distances = []
    for index in final_dist_index:
        final_distances.append(all_distances[index[0],index[1]])
        
    ###########################################################################
    end_time = time.time()
    print('Time in this dist function: {}'.format(end_time-start_time))
    ###########################################################################
    
    return final_distances
    
def cell_lower_triangular(struct,create_duplicate=True):
    '''
    Sets the cell back to lower triangular form
    Returns a boolean to indicate whether or not the reset was required
    '''
    if create_duplicate:
        struct=copy.deepcopy(struct)

    if (abs(struct.properties["lattice_vector_a"][1])<0.001 and 
            abs(struct.properties["lattice_vector_a"][2])<0.001 and 
	    abs(struct.properties["lattice_vector_b"][2])<0.001):
        return struct

    struct.properties.update(lattice_parameters(struct)) 
	#Add in case not calculated

    new_lattice=lattice_lower_triangular(struct)
    old_lattice=struct.get_lattice_vectors()
    rots = np.dot(np.transpose(new_lattice),np.linalg.inv(np.transpose(old_lattice)))
    cell_transform_mat(struct,rots,create_duplicate=False)
    struct.reset_lattice_vectors(new_lattice)
    return struct

def cell_modification(struct,napm=None,create_duplicate=True):
    '''
    Cell modification using Niggli reduction
    '''
    if create_duplicate:
        struct = copy.deepcopy(struct)
    lats = struct.get_lattice_vectors()
    from spglib import niggli_reduce
    reduced_lats =  niggli_reduce(lats)
    if reduced_lats is None:
        return False
    del(struct.properties["lattice_vector_a"])
    del(struct.properties["lattice_vector_b"])
    del(struct.properties["lattice_vector_c"])
    struct.set_lattice_vectors(reduced_lats)
    nmpc = int(len(struct.geometry)/napm)
    cell_lower_triangular(struct,False)
    move_molecule_in(struct,nmpc,False)
    #struct.set_lattice_angles()
    return struct

def cell_extension(struct,extension=[[0,0,1],[0,1,0],[0,1,1],
                                     [1,0,0],[1,0,1],[1,1,0],
                                     [1,1,1]],
                   create_duplicate=True):
    '''
    Extends the structure by the specified extensions
    '''

    if create_duplicate:
        struct=deepcopy(struct)
    napc = len(struct.geometry)
    for extend in extension:
        A = struct.properties["lattice_vector_a"]
        B = struct.properties["lattice_vector_b"]
        C = struct.properties["lattice_vector_c"]
        trans_vector = np.dot(extend,[A, B, C])

        #Create a new copy of each of the atom in the original geometry
        for i in range (napc):
            atom = copy.deepcopy(struct.geometry[i])
            struct.build_geo_by_atom(atom['x']+trans_vector[0],
                                     atom['y']+trans_vector[1],
                                     atom['z']+trans_vector[2],
                                     atom['element'],
                                     atom['spin'],
                                     atom['charge'],atom['fixed'])
    return struct

def move_molecule_in(struct,nmpc=None, create_duplicate=True):
	'''
	Translate the molecules by the cell vector such that their center of mass lies within the cell
	'''
	if create_duplicate:
		struct=deepcopy(struct)
	napm=int(len(struct.geometry)/nmpc)
	lattice=[struct.properties["lattice_vector_a"],
          struct.properties["lattice_vector_b"],
          struct.properties["lattice_vector_c"]]
	lattice=np.transpose(lattice)
	latinv=np.linalg.inv(lattice)    
	for i in range(int(nmpc)):
		cm=cm_calculation(struct,list(range(i*napm,i*napm+napm)))
		frac=np.dot(latinv,cm)
		for j in range (0,3):
			lat=lat_interp[j]
			vec=struct.properties[lat]
			if (frac[j]<-0.0001):
				kk=int(-frac[j]+1)
				for k in range(i*napm,i*napm+napm):
					for l in range (3):
						struct.geometry[k][l]+=kk*vec[l]
			elif (frac[j]>0.99999):
				kk=int(frac[j]+0.00001)
				for k in range(i*napm,i*napm+napm):
					for l in range (3):
						struct.geometry[k][l]-=kk*vec[l]
	return struct

def cm_calculation (struct,atom_list):
    '''
    Reads in a list of atom
    Find the center of mass
    '''
    cm = [0,0,0]
    tm = 0

    for i in range(len(atom_list)):
        tm += molar_mass[struct.geometry[atom_list[i]][3]]
        for j in range (3):
            cm[j] += molar_mass[struct.geometry[atom_list[i]][3]]*struct.geometry[atom_list[i]][j]
    for j in range (3):
        cm[j]/=tm
    return cm

# HARD CODED IN FOR NOW
nmpc = 4
if __name__ == '__main__':
#    test = all_interatomic_distances_manny(harringbone_list[0], 13, a1='C', a2='C')
#    test_old = all_interatomic_distances(harringbone_list[0], 13, a1='C', a2='C')
#    test_old = [x[2] for x in test_old]
    
    import matplotlib.pyplot as plt
    
    
    ###########################################################################
    # Very simple FCC test for RDF calculation                                #
    ###########################################################################
    struct = Structure()
    struct.set_lattice_vectors([[5,0,0],[0,5,0],[0,0,5]])
    struct.set_property('cell_vol',125)
    
    # FCC Structure
    struct.build_geo_by_atom_array(0, 0, 0, 'N')
    struct.build_geo_by_atom_array(2.5, 2.5, 0, 'N')
    struct.build_geo_by_atom_array(2.5, 0, 2.5, 'N')
    struct.build_geo_by_atom_array(0, 2.5, 2.5, 'N')
    
    # BCC Structure
#    struct.build_geo_by_atom_array(0, 0, 0, 'N')
#    struct.build_geo_by_atom_array(2.5, 2.5, 2.5, 'N')
    
    atomic_pairs_list=['N','N']
    atomic_distance_range=[0,15]
    atomic_distance_increment=0.5
    smoothing_parameter=1
    descriptor_key = 'rdf'
    
#    key_test = [x for x in harringbone_dict.keys()]
#    struct = harringbone_dict[key_test[0]]
    test_list = []
#    test_list.append(compute_RDF_vector(struct,atomic_pairs_list,
#                                             atomic_distance_range,
#                                             atomic_distance_increment,
#                                             smoothing_parameter))
    starting_struct_1 = deepcopy(struct)
    starting_struct_2 = deepcopy(struct)
    
    test_list.append(compute_RDF_vector_old(starting_struct_2,atomic_pairs_list,
                                     atomic_distance_range,
                                     atomic_distance_increment,
                                     smoothing_parameter,
                                     descriptor_key))
    
#    test_list.append(compute_RDF_vector_newest(starting_struct_1,atomic_pairs_list,
#                                             atomic_distance_range,
#                                             atomic_distance_increment,
#                                             smoothing_parameter))

    test_list.append(compute_RDF_vector_neighbor(starting_struct_1,atomic_pairs_list,
                                             atomic_distance_range,
                                             atomic_distance_increment,
                                             smoothing_parameter,
                                             descriptor_key)) 
        
    dist = np.arange(float(atomic_distance_range[0]), 
                     # Add atomic distance increment to make dist inclusive 
                     #  of the final distance value
                     float(atomic_distance_range[1])+atomic_distance_increment,
                     float(atomic_distance_increment))
    
    X = []
    for i,value in enumerate(dist[1:]):
        X.append(value-atomic_distance_increment/2)
    
        
    colors = ['tab:orange', 'tab:green', 
          'tab:red', 'tab:purple', 
          'tab:pink', 'tab:olive', 
          'tab:cyan']
        
    fig = plt.figure()
    fig.add_subplot(111)
    ax1 = fig.gca()
    plt.plot(X,test_list[0].get_property(descriptor_key),c='b')
    
    ax2= fig.add_subplot(111,sharex=ax1, frameon=False)
    plt.plot(X,test_list[1].get_property(descriptor_key),c='g')
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position('right')
    
    ###########################################################################
    # Very simple FCC test for RDF calculation                                #
    ###########################################################################
    
#    test_list = []
#    for i,struct_name in enumerate(harringbone_dict):
#        struct_temp = harringbone_dict[struct_name]
#        struct_temp = compute_RDF_vector(struct_temp,atomic_pairs_list,
#                                             atomic_distance_range,
#                                             atomic_distance_increment,
#                                             smoothing_parameter)
#        test_list.append(struct_temp)
#        if i == 0:
#            break
    
#    for i,struct_name in enumerate(graphite_dict):
#        struct_temp = graphite_dict[struct_name]
#        struct_temp = compute_RDF_vector(struct_temp,atomic_pairs_list,
#                                             atomic_distance_range,
#                                             atomic_distance_increment,
#                                             smoothing_parameter)
#        test_list.append(struct_temp)
#        if i == 1:
#            break
    
#    dist = np.arange(float(atomic_distance_range[0]), 
#                     # Add atomic distance increment to make dist inclusive 
#                     #  of the final distance value
#                     float(atomic_distance_range[1])+atomic_distance_increment,
#                     float(atomic_distance_increment))
#    X = []
#    for i,value in enumerate(dist[1:]):
#        X.append(value-atomic_distance_increment/2)
#
##    for struct in test_list:
##        fig = plt.figure()
##        plt.plot(X,struct.get_property('RDF_vector'),
##                 label=struct.struct_id)
##        plt.legend()
##        plt.show()
#        
#    fig = plt.figure()
#    fig.add_subplot(111)
#    ax1 = fig.gca()
#    plt.plot(X,test_list[0].get_property('RDF_vector'))
#    
#    ax2= fig.add_subplot(111,sharex=ax1, frameon=False)
#    plt.plot(X,test_list[1].get_property('RDF_vector'),c='g')
#    ax2.yaxis.tick_right()
#    ax2.yaxis.set_label_position('right')