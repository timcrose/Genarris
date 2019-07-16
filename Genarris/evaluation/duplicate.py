"""
If any part of this module is used for a publication please cite:

X. Li, F. Curtis, T. Rose, C. Schober, A. Vazquez-Mayagoitia, K. Reuter,
H. Oberhofer, and N. Marom "Genarris: Random Generation of Molecular Crystal 
Structures and Fast Screening with a Harris Approximation, ",
J. Chem. Phys., DOI: 10.1063/1.5014038; arXiv 1803.02145 (2018)
"""
'''
Created on 12/30/2016 by Patrick Kilecdi
'''


from pymatgen import Lattice as LatticeP
from pymatgen import Structure as StructureP
from pymatgen.analysis.structure_matcher import StructureMatcher,ElementComparator,SpeciesComparator,FrameworkComparator
import multiprocessing, copy
from Genarris.utilities import misc


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

def find_duplicates(inst):
    '''
    Module wrapper
    '''
    sname = "find_duplicates"
    structure_dir = inst.get(sname, "structure_dir")
    structure_dir_depth = inst.get_with_default(
            sname, "structure_dir_depth", 0, eval=True)
    structure_suffix = inst.get_with_default(sname, "structure_suffix", "")
    coll = misc.input_pool(structure_dir,
            structure_dir_depth=structure_dir_depth,
            structure_suffix=structure_suffix)

    ltol = inst.get_with_default(sname,"ltol",0.2,eval=True)
    stol = inst.get_with_default(sname,"stol",0.3,eval=True)
    angle_tol = inst.get_with_default(sname,"angle_tol",3,eval=True)
    scale_vol = inst.get_boolean(sname,"scale_vol")
    output_file = inst.get_with_default(sname,"duplicate_output_file","./duplicate.info")
    processes = inst.get_processes_limit(sname)

    dup_pairs = find_duplicate_pairs(coll,ltol,stol,angle_tol,scale_vol,processes)

    f = open(output_file,"a")
    f.write("Found %i pairs of duplicates:\n" % len(dup_pairs))
    f.write("\n".join(map(str,["%s %s" % x for x in dup_pairs])))
    f.write("\n")

    coll = remove_duplicates(coll,dup_pairs,make_copy=False)
    f.write("%i unique structures after removing duplicates:\n" % len(coll))
    f.write("\n".join(map(str,[x.struct_id for x in coll])))
    f.write("\n") 

    output_dir = inst.get_or_none(sname, "output_dir")
    output_format = inst.get_with_default(sname, "output_format", "both")
    if not output_dir is None:
        misc.output_pool(coll, output_dir, output_format)


def find_duplicates_distance_matrix(inst):
    '''
    Module wrapper
    '''
    sname = "find_duplicates_distance_matrix"
    coll = misc.load_collection_with_inst(inst, sname)
    dist_mat_path = inst.get(sname,"dist_mat_input_file")
    threshold = inst.get_eval(sname,"threshold")
    output_file = inst.get_with_default(sname,"duplicate_output_file","./duplicate.info")

    f = open(dist_mat_path,"r")
    lines = f.read().split("\n")
    while lines[-1]=="":
        lines.pop()
    dist_mat = [[float(x) for x in y.split()] for y in lines]

    dup_pairs = duplicate_distance_matrix(coll, dist_mat, threshold)

    f = open(output_file,"a")
    f.write("Using distance matrix: " + dist_mat_path)
    f.write("Found %i pairs of duplicates:\n" % len(dup_pairs))
    f.write("\n".join(map(str,["%s %s" % x for x in dup_pairs])))
    f.write("\n")

    coll = remove_duplicates(coll,dup_pairs,make_copy=False)
    f.write("%i unique structures after removing duplicates:\n" % len(coll))
    f.write("\n".join(map(str,[x.struct_id for x in coll])))
    f.write("\n") 

    if inst.has_option(sname,"output_folder"):
        misc.dump_collection_with_inst(inst,sname,coll)


def duplicate_distance_matrix(coll, dist_mat, threshold):
    '''
    Given a collection and a distance matrix, returns the pairs of duplicates
    '''
    coll.sort(key=lambda struct: struct.struct_id)

    pairs =[]
    for i in range(len(coll)-1):
        for j in range(i+1, len(coll)):
            if dist_mat[i][j]<threshold:
                pairs.append((coll[i].struct_id, coll[j].struct_id))

    return pairs

def remove_duplicates(coll, dup_pairs,make_copy=True):
    '''
    Given a collection, a list of pairs of duplicate structures' ID
    Remove until no pairs co-exist in the collection
    '''
    if make_copy:
        coll = copy.deepcopy(coll)

    for pair in dup_pairs:
        id1, id2 = pair
        s1 = [x for x in coll if x.struct_id==id1]
        if len(s1)==0: continue;

        s2 = [x for x in coll if x.struct_id==id2]
        if len(s2)==0: continue;
        coll.pop(coll.index(s2[0]))

    return coll

def find_duplicate_pairs(coll, ltol=.2, stol=.3, angle_tol=3, scale_vol=False, processes=1):
    '''
    Give a collection of structures
    Return all pairs of duplicates
    '''
    arglist = [(s1, s2, ltol, stol, angle_tol, scale_vol)
                for s2 in coll
                for s1 in coll if coll.index(s1)<coll.index(s2)]

    if processes ==1:
        result = [compare_struct_multiprocessing_wrapper(x) for x in arglist]
    else:
        pool = multiprocessing.Pool(processes)
        result = pool.map(compare_struct_multiprocessing_wrapper,arglist)
    return [(x[0].struct_id, x[1].struct_id) for x in result if x[2]==True]

def compare_struct(s1, s2, ltol=.2, stol=.3, angle_tol=3, scale_vol=False):
    '''
    Compares the 2 structures using pymatgen 
    '''
    sp1 = get_pymatgen_structure(s1)
    sp2 = get_pymatgen_structure(s2)
    sm = StructureMatcher(ltol=ltol, stol=stol, angle_tol=angle_tol, primitive_cell=True,
                          scale=scale_vol, 
                          attempt_supercell=False, 
                          comparator=SpeciesComparator())
    return sm.fit(sp1, sp2)

def compare_struct_multiprocessing_wrapper(args):
    s1, s2, ltol, stol, angle_tol, scale_vol = args
    return (s1, s2, compare_struct(s1, s2, ltol, stol, angle_tol, scale_vol))

def get_pymatgen_structure(struct):
    '''
    Args: Geometric data from GAtor's Structure() object
    Returns: A pymatgen StructureP() object with the same geometric properties
    '''
    frac_data = struct.get_frac_data()
    coords = frac_data[0] # frac coordinates
    atoms = frac_data[1] # site labels
    lattice = LatticeP.from_parameters(a=frac_data[2], 
                                       b=frac_data[3], 
                                       c=frac_data[4],
                                       alpha=frac_data[5],
                                       beta=frac_data[6], 
                                       gamma=frac_data[7])
    structp = StructureP(lattice, atoms, coords)
    return structp
