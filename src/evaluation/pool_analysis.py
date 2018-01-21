'''
Created by Patrick Kilecdi

This module outputs certain results of a pool of structure
'''

import numpy, multiprocessing, copy, math
import numpy as np
from utilities import write_log, misc
from core import structure_handling, pool_management
import random, os
from evaluation.evaluation_util import BatchSingleStructureOperation, \
        load_batch_single_structure_operation_keywords

def niggli_reduction_batch(inst):
    '''
    Conduct batch niggli reduction
    '''
    sname = "niggli_reduction_batch"
    napm = inst.get_eval(sname,"NAPM")
    kwargs = load_batch_single_structure_operation_keywords(
            inst, sname)

    op = BatchSingleStructureOperation(
            structure_handling.cell_niggli_reduction,
            name=sname, args=(napm,),
            kwargs={"create_duplicate" : False}, **kwargs)

    return op.run()

def specific_radius_batch(inst):
    '''
    Conduct batch calculation of specific radius
    '''
    sname = "specific_radius_batch"
    napm = inst.get_eval(sname, "NAPM")
    kwargs = load_batch_single_structure_operation_keywords(
            inst, sname)
    property_name = inst.get_with_default(sname, "property_name", "sr")
    op = BatchSingleStructureOperation(
            _specific_radius_calculation,
            name=sname, args=(napm,),
            kwargs={"property_name" : property_name},
            **kwargs)

    return op.run()

def _specific_radius_calculation(struct, napm, property_name="sr"):
    nmpc = struct.get_n_atoms() / napm
    struct.properties[property_name] = \
            structure_handling.specific_radius_check(struct, nmpc)
    return struct

def pool_single_structure_analysis(inst):
    '''
    A module that conducts batch operation on pools
    '''
    info_level = inst.get_info_level()

    sname = "pool_single_structure_analysis" #Module name
    if inst.has_option(sname,"structure_path"):
        struct_list = [inst.get(sname,"structure_path")]
    elif inst.has_option(sname,"structure_dir"):
        struct_list = misc.list_directory_file(inst.get(sname,"structure_dir"),\
        inst.get_with_default(sname,"structure_suffix",""),\
        inst.get_with_default(sname,"structure_dir_depth",0,eval=True))
    else:
        raise ValueError("Both structure_path and structure_dir missing")

    if info_level>=2 or (inst.has_procedure("Pool_Single_Structure_Analysis")\
    and info_level>=1):
        write_log.write_master_log(inst,"pool_single_structure_analysis has been called to evaluate %i structures" % len(struct_list))

    processes_limit = inst.get_processes_limit(sname)

    change_folder = False 
    try: #Determine whether overwriting original structure is called
        if os.path.abspath(inst.get(sname,"output_folder")) ==\
           os.path.abspath(inst.get(sname,"structure_dir")):
            change_folder = True
            if info_level>=3:
                write_log.write_master_log(inst,"Overwriting printing begins")
    except:
        pass

    if processes_limit > 1:
        mulp = multiprocessing.Pool(processes=processes_limit)
        inst_list = []
        new_inst = copy.deepcopy(inst)
        new_inst.set(sname,"processes_limit",str(info_level-1))

        for path in struct_list:
            new_inst = copy.deepcopy(new_inst)
            new_inst.set(sname,"structure_path",path)
            if change_folder:
                new_inst.set(sname,"output_folder",os.path.dirname(path))
            inst_list.append(new_inst)
        
        result = mulp.map(pool_single_structure_analysis, inst_list)
        return result
    ####################### Main FUnction ################################3
    functions = inst.get_list(sname,"functions")
    from evaluation import single_structure_analysis
    all_results = []
    for path in struct_list:
        struct = misc.load_structure_with_inst(inst,sname,"structure_format",path)
        result = []
        for function in functions:
            struct, nr = getattr(single_structure_analysis,function)(inst,struct)
            if isinstance(nr,tuple):
                result += list(nr)
            else:
                result.append(nr)
        if inst.has_option(sname,"output_info_file"):
            message = struct.struct_id + "  "+"  ".join(map(str,result))
            write_log.write_log(inst.get(sname,"output_info_file"),message,time_stamp=False)
        if inst.has_option(sname,"output_folder"):
            if change_folder:
                inst.set(sname,"output_folder",os.path.dirname(path))
            misc.dump_structure_with_inst(inst,sname,struct)

        all_results.append(result)


    if len(all_results)==1:
        return all_results[0]
    else:
        return all_results
            

def vector_distance_calculation(inst):
    '''
    Calculates and compare a given vector in the structure list
    '''
    sname = "vector_distance_calculation"
    coll = misc.load_collection_with_inst(inst,sname)
    coll.sort(key=lambda struct:struct.struct_id)

    vpn = inst.get(sname,"vector_property_name")
    dt = inst.get(sname,"distance_type")
    out_list = inst.get_with_default(sname,"dist_list_output","./vector_dist_list.info")
    out_matr = inst.get_with_default(sname,"dist_matrix_output",None)
    processes = inst.get_processes_limit(sname)
    info_level = inst.get_info_level(sname)

    if info_level > 3:
        verbose = True
    else:
        verbose = False

    if info_level > 1:
        print("Conducting vector distance calculation for %i structures" % len(coll))
        print("Using vector property: " + vpn)
    
    X = [struct.properties[vpn] for struct in coll]
    dist_list = get_distance_list(X, distance_type=dt, processes=processes, verbose=verbose)

    f = open(out_list,"a")
    for d in dist_list:
        print d
        f.write("%s %s %f\n" % (coll[d[0]].struct_id,
                                coll[d[1]].struct_id,
                                d[2]))
    f.close()

    if out_matr!=None:

        diff_mat = _convert_list_to_mat(dist_list, len(coll))

        f = open(out_matr,"a")       
        for l in diff_mat:
            f.write(" ".join(map(str,l))+"\n")
        f.close()

    if info_level > 1:
        print("Calculation completed!")
    

def get_distance_list(X, distance_type="euclidean",processes=1,verbose=True):
    
    x = len(X)
    dist_list = []
    for i in range(x-1):
        if processes==1:
            for j in range(i+1,x):
                dist_list.append(_calculate_vector_distance((X[i],X[j],i,j,distance_type)))
        else:
            m = multiprocessing.Pool(processes=processes)
            dist_list += m.map(_calculate_vector_distance, 
                               [(X[i],X[j],i,j,distance_type) for j in range(i+1,x)])
            m.close()
            m.join()
        if verbose:
            print("Distance calculation done for index: " + str(i))

    return dist_list


def _convert_list_to_mat(dist_list,x):

    dist_mat = [[None for i in range(x)] for j in range(x)]

    for i in range(x):
        dist_mat[i][i] = 0

    for l in dist_list:
        dist_mat[l[0]][l[1]] = l[2]
        dist_mat[l[1]][l[0]] = l[2]

    return dist_mat


def parallel_property_comparison(coll_1,coll_2,property_name="",distance_type="cosine",avoid_same_structure=True,processes_limit=20):
    '''
    Conducts comparison of vectors in parallel
    '''
    pairs = [[struct_1,struct_2,property_name,distance_type] for struct_1 in coll_1 for struct_2 in coll2 \
    if not avoid_same_structure or struct_1.struct_id==None or struct_2.struct_id==None
    or struct_1.struct_id!=struct_2.struct_id]
    m = multiprocessing.Pool(processes=processes_limit)
    return m.map(calculate_property_distance,pairs)
    

def calculate_property_distance(args):
    struct_1, struct_2, property_name, distance_type = tuple(args)
    result = [struct_1.struct_id, struct_2.struct_id]
    if distance_type == "cosine":
        result.append(1-np.dot(struct_1.properties[property_name],struct_2.properties[property_name])/np.linalg.norm(struct_1.properties[property_name])/np.linalg.norm(struct_2.properties[property_name]))
    elif distance_type == "euclidean":
        result.append(np.linalg.norm(np.subtract(struct_1.properties[property_name],struct_2.properties[property_name])))
    elif distance_type == "normalized euclidean":
        result.append(np.linalg.norm(np.subtract(struct_1.properties[property_name],struct_2.properties[property_name]))
                      /np.linalg.norm(struct_1.properties[property_name])/np.linalg.norm(struct_2.properties[property_name]))

    else:
        raise ValueError("Unsupported distance type; Only supporting 'cosine' and 'euclidean' now")

def _calculate_vector_distance(args):
    v1, v2, i1, i2, distance_type = tuple(args)
    if distance_type == "cosine":
        dist = 1-np.dot(v1, v2)/np.linalg.norm(v1)/np.linalg.norm(v2)
    elif distance_type == "euclidean":
        dist = np.linalg.norm(np.subtract(v1,v2))
    elif distance_type == "normalized euclidean":
        dist = np.linalg.norm(np.subtract(v1, v2))/np.linalg.norm(v1)/np.linalg.norm(v2)

    else:
        raise ValueError("Unsupported distance type; Only supporting 'cosine', 'euclidean' and 'normalized euclidean' now")

    return (i1, i2, dist)

def calculate_distance(v1, v2, distance_type="cosine"):     
    if distance_type == "cosine":
        return 1-np.dot(v1, v2)/np.linalg.norm(v1)/np.linalg.norm(v2)
    elif distance_type == "euclidean":
        return np.linalg.norm(np.subtract(v1, v2))
    elif distance_type == "normalized euclidean":
        return np.linalg.norm(np.subtract(v1,v2))/np.linalg.norm(v1)/np.linalg.norm(v2)
    else:
        raise ValueError("Unsupported distance type; Only supporting 'cosine', 'euclidean' and 'normalized euclidean' now")

