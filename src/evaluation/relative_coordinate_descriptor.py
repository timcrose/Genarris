'''
Created on 12/30/2016

By Patrick Kilecdi

Assigns a new descriptor to crystal structures,
and conducts Affinity Propagation for pools
'''

import multiprocessing
from core.structure import Structure
import numpy as np
from copy import deepcopy
from utilities import misc, write_log
import os, random, time, socket

from core.structure_handling import cm_calculation, cell_modification, \
        mole_translation
from evaluation.evaluation_util import BatchSingleStructureOperation, \
        load_batch_single_structure_operation_keywords

def rcd_calculation(inst):
    '''
    Main RCD calculation module
    '''
    sname = "rcd_calculation"
    kwargs = load_batch_single_structure_operation_keywords(inst, sname)

    napm = inst.get_eval(sname, "NAPM")
    axes = inst.get_eval(sname,"axes")
    close_picks = inst.get_with_default(sname, "close_picks", 16, eval=True)
    property_name = inst.get_with_default(
            sname, "property_name", "rcd")
    output_info_file = inst.get_or_none(sname, "output_info_file")

    op_args = (napm, axes)
    op_kwargs = {"close_picks"   : close_picks,
                 "property_name" : property_name}

    op = BatchSingleStructureOperation(
            _rcd_calculation,
            name=sname, args=op_args,
            kwargs=op_kwargs, **kwargs)
    result = op.run()

    if output_info_file != None:
        _rcd_calculation_callback(
                result, output_info_file)

    return result

def _rcd_calculation(struct, napm, axes, close_picks=16, property_name="rcd"):
    '''
    Generates and stores the rcd for the given structure
    '''
    if struct.get_n_atoms() %napm != 0:
        raise ValueError("Specified NAPM does not divide %s's atom list length"
                                 % struct.struct_id)
    nmpc = len(struct.geometry)/napm
 
    ref_struct = _create_ref_struct(struct,nmpc,napm,close_picks)

    axes_list = [
            _calculate_molecule_axes(ref_struct.geometry[x*napm:(x+1)*napm],
                                     axes)
            for x in range(close_picks+1)]

    COM = [cm_calculation(ref_struct,range(x*napm,(x+1)*napm))
           for x in range(close_picks+1)]

    # Calculates Cartesian relative coordinate
    diff = [np.subtract(COM[x], COM[0]) for x in range(1, close_picks+1)]
    
    # Calculates relative coordinate in axes basis of reference molecule
    axes_i = np.linalg.inv(axes_list[0].T)
    diff_t = [list(np.dot(axes_i,x)) for x in diff]
    
    orien = [[float(np.dot(axes_list[x][y], axes_list[0][y]))
              for y in range(3)]
             for x in range(1, close_picks+1)]

    descriptor = [[diff_t[x], orien[x]] for x in range(close_picks)]
    struct.set_property(property_name, descriptor)
    return struct, descriptor

def _create_ref_struct(struct, nmpc, napm, close_picks=16):
    '''
    Using the 1st molecule in the structure as reference
    Select and return a structure with the molecules that are closest to it

    struct: input structure
    nmpc: number of molecupe per cell
    napm: number of atom per cell
    close_picks: number of closest molecules to select

    '''
    ref_struct = cell_modification(struct, nmpc, napm)

    COM = [cm_calculation(ref_struct,
                          range(x*napm,(x+1)*napm))
           for x in range(nmpc)]


    lat_mat = np.transpose(ref_struct.get_lattice_vectors())
    COMt = [(k, x, y, z, 
            np.linalg.norm(np.subtract(np.add(COM[k],np.dot(lat_mat,[x,y,z])),
                                       COM[0]))) 
            for z in range(-2,3)
            for y in range(-2,3)
            for x in range(-2,3)
            for k in range(nmpc)] #Calculate the distance between the COM's
    COMt.sort(key=lambda com: com[4])

    all_geo = deepcopy(ref_struct.geometry)
    ref_struct.geometry = np.delete(ref_struct.geometry,
                                    range(napm,nmpc*napm))

    for i in range(1, close_picks+1): #Skip the 1st original molecule
        k, x, y, z, dist = COMt[i]
        ref_struct.geometry = np.concatenate((ref_struct.geometry,
                                              all_geo[k*napm:(k+1)*napm]))

        mole_translation(
                ref_struct, i, napm, frac=[x,y,z], create_duplicate=False)

    return ref_struct

def _calculate_molecule_axes(geo, axes):
    '''
    Calculates the orientation axes of the molecule
    geo: geometry slice of the molecule
    axes: a list of 2 tuples, each tuple include the indeces of two atoms
    The first 2 axes will be established as the 2 vectors formed by the
    pairs of atoms. The third will be their cross product
    '''
    a1 = [geo[axes[0][0]][x]-geo[axes[0][1]][x] for x in range(3)]
    a2 = [geo[axes[1][0]][x]-geo[axes[1][1]][x] for x in range(3)]
    a3 = list(np.cross(a1,a2))
    X = np.array([a1,a2,a3])
    X = _gram_schmidt(X)
    return X


def _gram_schmidt(X, row_vecs=True, norm=True):
    if not row_vecs:
        X = X.T
    Y = X[0:1,:].copy()
    for i in range(1, X.shape[0]):
        proj = np.diag((X[i,:].dot(Y.T)/np.linalg.norm(Y,axis=1)**2).flat).dot(Y)
        Y = np.vstack((Y, X[i,:] - proj.sum(0)))
    if norm:
        Y = np.diag(1/np.linalg.norm(Y,axis=1)).dot(Y)
    if row_vecs:
        return Y
    else:
        return Y.T

def _rcd_calculation_callback(result, output_info_file):
    result.sort(key=lambda x : x[0].struct_id)
    message = ""
    for struct, rcd in result:
        message += struct.struct_id + " "
        message += " ".join(map(str, _flatten_rcd(rcd)))
        message += "\n"
    write_log.write_log(output_info_file, message, time_stamp=False)

def _flatten_rcd(v):
    return [entry for mol in v for rel in mol for entry in rel]

def rcd_difference_calculation(inst):
    '''
    Main module of rcd difference
    '''
    sname = "rcd_difference_calculation"
    key = inst.get_with_default(sname,"stored_property_key","rcd")
    ratio = inst.get_with_default(sname,"contribution_ratio",1,eval=True)
    pairs = inst.get_with_default(sname,"select_pairs",4,eval=True)
    dis_en = inst.get_boolean(sname,"disable_enantiomer")
    diff_list_output = inst.get_with_default(sname,"diff_list_output",
                                             "./rcd_difference_list.info")
    diff_matrix_output = inst.get_with_default(
            sname, "diff_matrix_output", "./rcd_difference_matrix.info")
    processes = inst.get_processes_limit(sname)

    coll = misc.load_collection_with_inst(inst,sname)
    coll.sort(key=lambda struct: struct.struct_id)

    result, diff_mat = calculate_collection_rcd_difference(coll, key=key,
                                                           ratio=ratio, 
                                                           select_pairs=pairs,
                                                           allow_enantiomer=not dis_en, 
                                                           processes=processes)
  
    f = open(diff_list_output,"a") 

    for k in result:
        f.write("%s %s %f \n" % (coll[k[0]].struct_id, coll[k[1]].struct_id, k[2]))
    f.close()

    if diff_matrix_output!="":
        f = open(diff_matrix_output,"a")
        for k in diff_mat:
#            print k
            f.write(" ".join(map(str,k))+"\n")
        f.close()

    return result, diff_mat

def rcd_difference_compare_single(inst):
    '''
    Main module that compares a single structure with a group of structures
    '''
    sname = "rcd_difference_compare_single"
    key = inst.get_with_default(sname,"stored_property_key","RCD_vector")
    ratio = inst.get_with_default(sname,"contribution_ratio",1,eval=True)
    pairs = inst.get_with_default(sname,"select_pairs",4,eval=True)
    dis_en = inst.get_boolean(sname,"disable_enantiomer")
    diff_list_output = inst.get_with_default(sname,"diff_list_output",
                                             "./rcd_difference_list.info")
    path = inst.get(sname,"structure_path")
#    diff_matrix_output = inst.get_with_default(sname,"diff_matrix_output","")
    processes = inst.get_processes_limit(sname)

    coll = misc.load_collection_with_inst(inst,sname)
    coll.sort(key=lambda struct: struct.struct_id)
    
    f = open(path,"r")
    struct = Structure()
    struct.loads(f.read())
    f.close()

    result = calculate_single_rcd_difference(struct, coll, key, ratio,
                                            select_pairs=pairs, 
                                            allow_enantiomer=not dis_en,
                                            processes=processes)
    
    f = open(diff_list_output,"a")
    for i in range(len(coll)):
        f.write("%s %s %f\n" % (struct.struct_id, 
                                coll[i].struct_id,
                                result[i]))
    f.close()

def rcd_difference_folder_inner(inst):
    sname = "rcd_difference_folder_inner"
    folder = inst.get(sname,"structure_dir")
    suffix = inst.get_with_default(sname,"structure_suffix","")

    key = inst.get_with_default(sname,"stored_property_key","RCD_vector")
    ratio = inst.get_with_default(sname,"contribution_ratio",1,eval=True)
    pairs = inst.get_with_default(sname,"select_pairs",4,eval=True)
    dis_en = inst.get_boolean(sname,"disable_enantiomer")

    calculate_folder_rcd_difference_worker(folder, suffix=suffix,
                                           key=key, ratio=ratio, select_pairs=pairs,
                                           allow_enantiomer=not dis_en)

        
def calculate_folder_rcd_difference_worker(folder, suffix="", key="RCD_vector", 
                                           ratio=1, select_pairs=4, allow_enantiomer=True):
    '''
    Given a folder, grab structure that has not been evaluated
    '''
    time.sleep(random.random()*5)
    f = open("./RCD_report.out","a")
    name = socket.gethostname()+"_"+misc.get_random_index()
    f.write("RCD_difference_worker %s reporting from host: %s\n" % 
            (name, socket.gethostname()))
    f.close()

    list_of_struct = _get_structure_list(folder,suffix)
    while True:
        struct, struct_path = _get_structure(folder, suffix,worker=socket.gethostname())
        if struct == False:
            break
        f = open("./RCD_report.out","a")
        f.write("RCD_difference_worker %s grabed structure %s\n" %
                (name,struct_path))
        f.close()

        diff_list = []
        for x in list_of_struct:
            compared_struct = Structure()
            compared_struct.loads_try_both(x)
            diff = calculate_rcd_difference(struct.properties[key],
                                            compared_struct.properties[key],
                                            ratio=ratio, select_pairs=select_pairs,
                                            allow_enantiomer=allow_enantiomer)
 
            diff_list.append((struct.struct_id, compared_struct.struct_id, diff))

        diff_list.sort(key=lambda x: x[1])
        f = open(struct_path+".rcd","w")
        for x in diff_list:
            f.write("%s %s %f\n" % (x[0], x[1], x[2]))
        f.close()
        try:
            os.remove(struct_path+".lock")
        except:
            pass

def calculate_collection_rcd_difference(coll, key="RCD_vector",
                                        ratio=1, select_pairs=4, allow_enantiomer=True,
                                        processes=1):
    '''
    Calculates the rcd difference between each pairs of structures
    Outputs: a tuple of difference list and difference matrix
    '''
    result = []
    for i in range(len(coll)):
        arglist = [(i, j, coll[i].properties[key], coll[j].properties[key],
                    ratio, select_pairs, allow_enantiomer)
                    for j in range(i,len(coll))]

        if processes > 1:
            p = multiprocessing.Pool(processes)
            result += p.map(_calculate_rcd_difference_multiprocess_wrapper,arglist)
            p.close()
            p.join()
        else:
            result += [_calculate_rcd_difference_multiprocess_wrapper(args) 
                                                    for args in arglist]

    diff_mat = [[0]*len(coll) for x in range(len(coll))]
    for k in result: 
        diff_mat[k[0]][k[1]] = k[2]
        diff_mat[k[1]][k[0]] = k[2]
    return result, diff_mat

def calculate_single_rcd_difference(struct, coll, key="RCD_vector",
                                        ratio=1, select_pairs=4, allow_enantiomer=True,
                                        processes=1):
    '''
    Calculates the rcd difference between each pairs of structures
    Outputs: a tuple of difference list and difference matrix
    '''
    result = []
    arglist = [(-1, j, struct.properties[key], coll[j].properties[key])
                for j in range(len(coll))]

    if processes > 1:
        p = multiprocessing.Pool(processes)
        result += p.map(_calculate_rcd_difference_multiprocess_wrapper,arglist)
        p.close()
        p.join()
    else:
        result += [_calculate_rcd_difference_multiprocess_wrapper(args) 
                                                for args in arglist]

    diff_list = [0]*len(coll)
    for k in result:
        diff_list[k[1]] = k[2]
    
    return diff_list


def calculate_collection_rcd_difference_2(coll, key="RCD_vector",
                                        ratio=1, select_pairs=4, allow_enantiomer=True,
                                        processes=1):
    '''
    Calculates the rcd difference between each pairs of structures
    Outputs: a tuple of difference list and difference matrix
    '''
    print "This is before arglist"
    arglist = [(i, j, coll[i].properties[key], coll[j].properties[key],
                ratio, select_pairs, allow_enantiomer)
                for i in range(len(coll)) for j in range(i,len(coll))]
    print "This is after arglist"

    if processes > 1:
        p = multiprocessing.Pool(processes)
        result = p.map(_calculate_rcd_difference_multiprocess_wrapper,arglist)
    else:
        result = [_calculate_rcd_difference_multiprocess_wrapper(args) 
                                                    for args in arglist]

    diff_mat = [[0]*len(coll) for x in range(len(coll))]
    for k in result: 
        diff_mat[k[0]][k[1]] = k[2]
        diff_mat[k[1]][k[0]] = k[2]
    return result, diff_mat

   

def calculate_rcd_difference(v1, v2, ratio=1, select_pairs=4, 
                             allow_enantiomer=True):
    '''
    Given two rcd vectors, calculate their RCD_vector difference
    ratio: The ratio between the contributions to the final difference
           by relative orientation difference and relative coordinate difference
    allow_enantiomer: if True, then s2 will be mirror reflected 
                      by flipping the sign of its third relative coordinate.
                      This is necessary b/c the third reference axes is
                      generated with cross product, and thus will acquire a 
                      sign flip under mirror reflection  
    '''

    v1 = deepcopy(v1)
    v2 = deepcopy(v2)
    result = _calculate_rcd_difference_2(v1, v2, ratio=ratio,
                                                select_pairs=select_pairs)
    if not allow_enantiomer:
        return result
    
    for x in v2:
        x[0][2] = - x[0][2]
    resulte = _calculate_rcd_difference_2(v1, v2, ratio=ratio,
                                                 select_pairs=select_pairs)

    return min(result, resulte)

def _calculate_rcd_difference_multiprocess_wrapper(arglist):
    return (arglist[0], arglist[1], calculate_rcd_difference(*arglist[2:]))

    

def _calculate_rcd_difference(v1, v2, ratio=1, select_pairs=4):
    '''
    Given 2 RCD vectors, return their difference
    v1: RCD vector 1
    v2: RCD vector 2
    ratio: The difference will be the normalized Euclidean distance of relative coordinates 
           + Euclidean distance of relative orientation
    '''
    
    dist = [(x,y,_calculate_diff(v1[x],v2[y],ratio)) for y in range(len(v2))
                                                     for x in range(len(v1))]

    dist.sort(key=lambda x: x[2])

    result = 0
    s1 = [False] * len(v1)
    s2 = [False] * len(v2)
    
    dist_iter = iter(dist)
    for l in range(select_pairs):
        try:
            while True:
                k = dist_iter.next()
                if not s1[k[0]] and not s2[k[1]]:
                #Avoiding selecting same molecule from cell
                    result += k[2]
                    s1[k[0]] = True
                    s2[k[1]] = True
#                    print "Here is the selection", k
                    break
        except:
            raise RuntimeError("Not enough close picks to select the amount of pairs")

    return result

def _calculate_rcd_difference_2(v1, v2, ratio=1, select_pairs=4):
    '''
    Given 2 RCD vectors, return their difference
    v1: RCD vector 1
    v2: RCD vector 2
    ratio: The difference will be the normalized Euclidean distance of relative coordinates 
           + Euclidean distance of relative orientation
    This version requires that the closest 4 molecules be accounted for
    The difference is taken as the average between the attempts
    '''
    
    dist = [(x,y,_calculate_diff(v1[x],v2[y],ratio)) for y in range(len(v2))
                                                     for x in range(select_pairs)]

    dist.sort(key=lambda x: x[2])

    result = 0
    s1 = [False] * len(v1)
    s2 = [False] * len(v2)
    
    dist_iter = iter(dist)
    for l in range(select_pairs):
        try:
            while True:
                k = dist_iter.next()
                if not s1[k[0]] and not s2[k[1]]:
                #Avoiding selecting same molecule from cell
                    result += k[2]
                    s1[k[0]] = True
                    s2[k[1]] = True
#                    print "Here is the selection", k
                    break
        except:
            raise RuntimeError("Not enough close picks to select the amount of pairs")

    dist = [(x,y,_calculate_diff(v1[x],v2[y],ratio)) for y in range(select_pairs)
                                                     for x in range(len(v1))]

    dist.sort(key=lambda x: x[2])

    s1 = [False] * len(v1)
    s2 = [False] * len(v2)
    
    dist_iter = iter(dist)
    for l in range(select_pairs):
        try:
            while True:
                k = dist_iter.next()
                if not s1[k[0]] and not s2[k[1]]:
                #Avoiding selecting same molecule from cell
                    result += k[2]
                    s1[k[0]] = True
                    s2[k[1]] = True
#                    print "Here is the selection", k
                    break
        except:
            raise RuntimeError("Not enough close picks to select the amount of pairs")


    return result/2




def _calculate_diff(v1, v2, ratio=1):
    relative_coordinate_diff = (np.linalg.norm(np.subtract(v1[0],v2[0]))**2
                                /np.linalg.norm(v1[0])/np.linalg.norm(v2[0]))
    relative_orientation_diff = np.linalg.norm(np.subtract(v1[1],v2[1]))**2/3
    #Divided by 6 to normalize (b/c the cosine value range from -1 to 1)

#    print "These are relative coordinates", v1[0], v2[0]
#    print "This is relative_coordinate_diff", relative_coordinate_diff
#
#    print "These are relative_orientations", v1[1], v2[1]
#    print "This is relative_orientation_diff", relative_orientation_diff
    return relative_coordinate_diff + ratio*relative_orientation_diff


 
def _get_structure(folder,structure_suffix="",worker=None):
    '''
    Finds a structure currently under the tmp directory 
    Places a lock on it for use
    '''
    f = open("./RCD_report.out","a")
    f.write("RCD_difference_worker "+str(worker)+" beginning to list files\n")
    f.close()

    lists = os.listdir(folder)
    f = open("./RCD_report.out","a")
    f.write("RCD_difference_worker "+str(worker)+" finished listing, beginning to filter\n")
    f.close()



#    list_of_files = [name for name in lists \
#                     if os.path.isfile(os.path.join(folder,name))]

#    lists = os.listdir(folder)
#    f = open("./RCD_report.out","a")
#    f.write("RCD_difference_worker finished making sure all files, beginning to filter\n")
#    f.close()

    list_of_files = [name for name in lists \
                     if name[len(name)-len(structure_suffix):]==structure_suffix and
                     name[len(name)-5:]!=".lock" and name[len(name)-4:]!=".rcd" and
                     not name+".lock" in lists and not name+".rcd" in lists]
#                     not os.path.isfile(os.path.join(folder,name+'.lock')) and 
#                     not os.path.isfile(os.path.join(folder,name+".rcd"))]

    f = open("./RCD_report.out","a")
    f.write("RCD_difference_worker "+str(worker)+" finished filtering files\n")
    f.close()


    found = False
    while len(list_of_files)>0:
        chosen = int(random.uniform(0,len(list_of_files)))
        name = list_of_files[chosen]
        f = open("./RCD_report.out","a")
        f.write("RCD_worker looking at item %i of a list of %i files\n" % 
                (chosen, len(list_of_files)))
        f.close()

        if not os.path.isfile(os.path.join(folder, name+".lock")):

            f = open(os.path.join(folder, name+".lock"),"w")
            f.write(str(time.time()))
            f.close()

            f = open(os.path.join(folder,name),"r")
            struct = Structure()
            try:
                struct.loads(f.read())
            except:
                try:
                    struct.build_geo_whole_atom_format(f.read())
                except:
                    print("Failed to load structure: " + 
                          os.path.join(folder, name))

            if struct.struct_id == None:
                l = name[:]
                while l[len(l)-5:]==".next":
                    l = l[:len(l)-5]
                struct.struct_id = l

            return struct, os.path.join(folder,name)

        else:
	    
            list_of_files=list_of_files[:chosen]+list_of_files[chosen+1:]

    return False, False

def _get_structure_list(folder,structure_suffix=""):
    f = open("./RCD_report.out","a")
    f.write("RCD_difference_worker beginning to get structure list\n")
    f.close()

    list_of_files = [os.path.join(folder,name) for name in os.listdir(folder) \
                     if os.path.isfile(os.path.join(folder,name)) and 
                     name[len(name)-len(structure_suffix):]==structure_suffix and
                     name[len(name)-5:]!=".lock" and name[len(name)-4:]!=".rcd" and
                     name[len(name)-4:]!=".fin"]
    return list_of_files

    f = open("./RCD_report.out","a")
    f.write("RCD_difference_worker finished getting structure list\n")
    f.close()



   

def test_create_ref_struct():
    import os
    direct = "./"
    path = [os.path.join(direct,x) for x in os.listdir(direct)]
#    path = ["./duplicate_pair_1/target_2_0983_b52d9d4a9b.json","./duplicate_pair_1/target_2_0986_33d0a828d6.json"]
    s = []
    for i in range(2):
        f = open(path[i],"r")
        struct=Structure()
        struct.build_geo_whole_atom_format(f.read())
#        calculate_molecule_axes(struct.geometry[0:11],[(1,10),(3,9)])
        generate_relative_coordinate_descriptor(struct,4,11,[(1,10),(3,9)],close_picks=4)
        s.append(struct)
    print calculate_rcd_difference(s[0].properties["RCD_vector"],s[1].properties["RCD_vector"],select_pairs=2)
#        create_ref_struct(struct,4,11,close_picks=8)
if __name__=="__main__":
    test_create_ref_struct()
            

