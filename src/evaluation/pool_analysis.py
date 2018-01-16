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

def interatomic_distance_evaluation(inst):
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
            r = min_distance(l)
            l.remove(r)
        result.append(r)

    return (struct, result)

def _interatomic_distance_evaluation_callback(result, output_file):
    message = ""
    for case in result:
        struct = case[0]
        case_result = case[1]
        message += struct.struct_id
        for r in case_result:
            message += " " + " ".join(map(str, r))
        message += "\n"
    write_log.write_log(output_file, message, time_stamp=False)
            
def radial_distribution_function(inst):
    '''
    Returns a list of interatomic distances between pairs of certain species

    '''
    info_level = inst.get_info_level()
        
    
    sname = "radial_distribution_function" #Module name
    if inst.has_option(sname,"structure_path"):
        struct_list = [inst.get(sname,"structure_path")]
    elif inst.has_option(sname,"structure_dir"):
        struct_list = misc.list_directory_file(inst.get(sname,"structure_dir"),\
        inst.get_with_default(sname,"structure_suffix",""),\
        inst.get_with_default(sname,"structure_dir_depth",0,eval=True))
    else:
        raise ValueError("radial_distribution_function: Both structure_path and structure_dir missing")

    if info_level>=2 or (inst.has_procedure("Radial_Distribution_Function")\
    and info_level>=1):
        write_log.write_master_log(inst,"radial_distribution_function has been called to evaluate %i structures" % len(struct_list))

    processes_limit = inst.get_processes_limit(sname)
    if processes_limit > 1:
        mulp = multiprocessing.Pool(processes=processes_limit)
        inst_list = []
        new_inst = copy.deepcopy(inst)
        new_inst.set(sname,"processes_limit",str(info_level-1))
        for path in struct_list:
            new_inst = copy.deepcopy(new_inst)
            new_inst.set(sname,"structure_path",path)
            inst_list.append(new_inst)
        
        result = mulp.map(radial_distribution_function, inst_list)
        return result

    #List of parameters
    reference_atom = inst.get_eval(sname,"reference_atom")
    target_atom_type = inst.get_eval(sname,"target_atom_type")
    dist = inst.get_eval(sname,"distance_range")
    bins = inst.get_with_default(sname,"bins",10,eval=True)
    normalized = inst.get_boolean(sname,"normalized")
    if inst.has_option(sname,"output_folder"): #Only outputs json
        inst.set(sname,"output_format","['json']")
        if inst.has_option(sname,"output_suffix"): #Need to put it into the list format
            inst.set(sname,"output_suffix","['"+inst.get(sname,"output_suffix")+"']")
        property_name = inst.get_with_default(sname,"property_name","RDF_vector")
        
    if inst.has_option(sname,"output_info_file"):
        output_file = inst.get(sname,"output_info_file")
    else:
        output_file = False

    result = []
    for path in struct_list:
        original_struct = misc.load_structure_with_inst(inst,sname,"structure_format",path)
        struct = copy.deepcopy(original_struct)

        #First figure out the cell extension needed
        structure_handling.cell_modification(struct,struct.get_n_atoms(),\
                            1, False)
        structure_handling.cell_lower_triangular(struct,False)
        a_ext = int(math.ceil(dist[1]/struct.properties["lattice_vector_a"][0]))
        b_ext = int(math.ceil(dist[1]/struct.properties["lattice_vector_b"][1]))
        c_ext = int(math.ceil(dist[1]/struct.properties["lattice_vector_c"][2]))
        structure_extension = [[x,y,z] for x in range (-a_ext,a_ext+1)\
                           for y in range (-b_ext,b_ext+1)\
                           for z in range (-c_ext,c_ext+1)]
        structure_handling.cell_extension(struct,extension=structure_extension,create_duplicate=False)

        if info_level >= 3:
            print(struct.get_geometry_atom_format())

        if isinstance(reference_atom,int):
            a1_range = [reference_atom]
        elif isinstance(reference_atom,list):
            a1_range = reference_atom
        elif reference_atom == "X":
            a1_range = range(struct.get_n_atoms())
        elif isinstance(reference_atom,str):
            a1_range = [i for i in range (original_struct.get_n_atoms())\
            if original_struct.geometry[i]["element"]==reference_atom]
        else:
            raise ValueError("Unsupported format of reference_atom input")

        if isinstance(target_atom_type,str):
            target_atoms = [target_atom_type]
        elif isinstance(target_atom_type,list):
            target_atoms = target_atom_type[:]
        else:
            raise ValueError("Unsupported format of target_atom_type input")
        
        vector_all = [struct.struct_id]
        for target_atom_type in target_atoms:
            result = all_interatomic_distances(struct, "X", target_atom_type,a1_range=a1_range)
            distances = [case[2] for case in result]
            vector, bins_1 = np.histogram(distances,bins=bins,range=dist)
            vector = list(vector)
            bins_1 = list(bins_1)
            vector = [x/len(a1_range) for x in vector]
            if normalized:
                for i in range (len(vector)):
                    vector[i] = vector[i]/(4/3.0*math.pi*(bins_1[i+1]**3-bins_1[i]**3))
            vector_all += vector[:]
    
        if not output_file == False:
            message = " ".join(map(str,vector_all))
            write_log.write_log(output_file,message,time_stamp=False)

        if inst.has_option(sname,"output_folder"):
            if property_name in original_struct.properties:
                original_struct.properties[property_name] += vector_all[1:]
            else:
                original_struct.properties[property_name] = vector_all[1:]
            if inst.get_boolean(sname,"add_description"):
                message = "reference atom indices: %s; target atoms: %s; bins: %s. " % (str(a1_range), str(target_atoms), str(bins_1))
                try:
                    original_struct.properties["RDF_vector_description"] += message
                except:
                    original_struct.properties["RDF_vector_description"] = message

            original_struct.properties[property_name] = vector_all[1:]
            misc.dump_structure_with_inst(inst,sname,original_struct)
        
        result.append(vector_all)
    return result       
        
def radial_distribution_function_1(inst):
    '''
    Returns a list of interatomic distances between pairs of certain species

    '''
    info_level = inst.get_info_level()
        
    
    sname = "radial_distribution_function_1" #Module name
    if inst.has_option(sname,"structure_path"):
        struct_list = [inst.get(sname,"structure_path")]
    elif inst.has_option(sname,"structure_dir"):
        struct_list = misc.list_directory_file(inst.get(sname,"structure_dir"),\
        inst.get_with_default(sname,"structure_suffix",""),\
        inst.get_with_default(sname,"structure_dir_depth",0,eval=True))
    else:
        raise ValueError("radial_distribution_function: Both structure_path and structure_dir missing")

    if info_level>=2 or (inst.has_procedure("Radial_Distribution_Function")\
    and info_level>=1):
        write_log.write_master_log(inst,"radial_distribution_function has been called to evaluate %i structures" % len(struct_list))

    processes_limit = inst.get_processes_limit(sname)
    if processes_limit > 1:
        mulp = multiprocessing.Pool(processes=processes_limit)
        inst_list = []
        new_inst = copy.deepcopy(inst)
        new_inst.set(sname,"processes_limit",str(info_level-1))
        change_folder = False
        try:
            if os.path.abspath(inst.get(sname,"output_folder")) ==\
               os.path.abspath(inst.get(sname,"structure_dir")):
                change_folder = True
        except:
            pass

        for path in struct_list:
            new_inst = copy.deepcopy(new_inst)
            new_inst.set(sname,"structure_path",path)
            if change_folder:
                new_inst.set(sname,"output_folder",os.path.dirname(path))
            inst_list.append(new_inst)
        
        result = mulp.map(radial_distribution_function_1, inst_list)
        return [x[0] for x in result]

    #List of parameters
    atomic_pairs = inst.get_eval(sname,"atomic_pairs")
    if isinstance(atomic_pairs[0],str):
        atomic_pairs = [atomic_pairs]
        
    dist = eval(inst.get(sname,"sampled_distances"))
    smoothing_parameter = inst.get_with_default(sname,"smoothing_parameter",1,eval=True)
    normalized = inst.get_boolean(sname,"normalized")
    if normalized:
        a = smoothing_parameter
        n_factor = [4*math.pi*(0.25*(math.pi/a**3)**0.5+\
        0.5*r**2*(math.pi/a)**0.5+3*r/a*math.exp(-a*r**2)+\
        (0.25*(math.pi/a**3)**0.5+0.5*r**2*(math.pi/a)**0.5)*math.erf(r*a**0.5)) \
        for r in dist]

    if inst.has_option(sname,"output_folder"): #Only outputs json
        inst.set(sname,"output_format","['json']")
        if inst.has_option(sname,"output_suffix"): #Need to put it into the list format
            inst.set(sname,"output_suffix","['"+inst.get(sname,"output_suffix")+"']")
        property_name = inst.get_with_default(sname,"property_name","RDF_vector")
        if os.path.abspath(inst.get(sname,"output_folder")) == os.path.abspath(inst.get(sname,"structure_dir")):
            over_write = True
            #If the output_folder is the same as input folder, then overwrites original structures
        else:
            over_write = False
        
    if inst.has_option(sname,"output_info_file"):
        output_file = inst.get(sname,"output_info_file")
    else:
        output_file = False


    result = []
    for path in struct_list:
        original_struct = misc.load_structure_with_inst(inst,sname,"structure_format",path)
        if info_level >= 3:
            print(original_struct)
        struct = copy.deepcopy(original_struct)

        #First figure out the cell extension needed
        structure_handling.cell_modification(struct,struct.get_n_atoms(),\
                            1, False)
        structure_handling.cell_lower_triangular(struct,False)
        a_ext = int(math.ceil(max(dist)/struct.properties["lattice_vector_a"][0]))
        b_ext = int(math.ceil(max(dist)/struct.properties["lattice_vector_b"][1]))
        c_ext = int(math.ceil(max(dist)/struct.properties["lattice_vector_c"][2]))
        structure_extension = [[x,y,z] for x in range (-a_ext-1,a_ext+2)\
                           for y in range (-b_ext-1,b_ext+2)\
                           for z in range (-c_ext-1,c_ext+2)]
        structure_handling.cell_extension(struct,extension=structure_extension,create_duplicate=False)

        if info_level >= 3:
            print(struct.get_geometry_atom_format())
            print(struct.get_n_atoms())

        vector_all = [struct.struct_id]

        for pairs in atomic_pairs:
            reference_atom = pairs[0]
            target_atom_type = pairs[1]
            if isinstance(reference_atom,int):
                a1_range = [reference_atom]
            elif isinstance(reference_atom,list):
                a1_range = reference_atom
            elif reference_atom == "X":
                a1_range = range(struct.get_n_atoms())
            elif isinstance(reference_atom,str):
                a1_range = [i for i in range (original_struct.get_n_atoms())\
                if original_struct.geometry[i]["element"]==reference_atom]
            else:
                raise ValueError("Unsupported format of reference_atom input")

            if not isinstance(target_atom_type,str):
                raise ValueError("Second atom in the atomic pair must by a specified species")
            rs = all_interatomic_distances(struct, "X", target_atom_type,a1_range=a1_range)
            distances = [case[2] for case in rs]
            g = [sum([math.exp(-smoothing_parameter*(r-r_ij)**2) for r_ij in distances])/len(a1_range) for r in dist]
            
            if normalized:
#               g = [g[i]/(4*math.pi*dist[i]**2) for i in range(len(dist))]
                g = [g[i]/n_factor[i] for i in range (len(dist))]
            vector_all += g[:]
    
        if not output_file == False:
            message = " ".join(map(str,vector_all))
            write_log.write_log(output_file,message,time_stamp=False)

        if inst.has_option(sname,"output_folder"):
            if property_name in original_struct.properties:
                original_struct.properties[property_name] += vector_all[1:]
            else:
                original_struct.properties[property_name] = vector_all[1:]
            if inst.get_boolean(sname,"add_description"):
                message = "atomic pairs = %s; sampled distances = %s. " % (str(atomic_pairs), str(dist))
                try:
                    original_struct.properties["RDF_vector_description"] += message
                except:
                    original_struct.properties["RDF_vector_description"] = message
            if over_write:
                inst.set(sname,"output_folder",os.path.dirname(path))
                #This honors the original structure of the input folder
            misc.dump_structure_with_inst(inst,sname,original_struct)

        result.append(vector_all)
    return result

# TODO: Refactor all interatomic distance calculation into an
# InteratomicDistanceCalculation class

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
    napc = len(struct.geometry)

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
            return all_interatomic_distances(struct, a1_element, a2,\
                a1_range=[a1], a2_range=range(0,napc), napm=napm)
        else:
            result = []
            # Iterate through all the molecules except the one a1 is in
            for i in range(0, a1//napm) + range(a1//napm+1, nmpc):
                result.append(min_distance(all_interatomic_distances(\
                struct,struct.geometry[a1]["element"],a2,\
                a1_range=[a1], a2_range=range(i*napm,(i+1)*napm))))
            return [x for x in result if x!=None]

    elif isinstance(a1,str) and isinstance(a2,int):
        a2_element = struct.geometry[a2]["element"]
        if not block_same_pair:
            if not allow_all_pairs:
                if a2 < napm:
                    # Immediately returns because atoms must be from different
                    # molecules. However, if a1 must be selected from the first
                    # molecule since allow_all_pairs=False, and a2 is from the
                    # first molecule as well, then no result can be returned
                    return []
                return all_interatomic_distances(struct, a1, a2_element, \
                        a1_range=range(0,napm), a2_range=[a2])
            else:
                # Avoids a1_range from including any atoms of a2 molecule
                return all_interatomic_distances(struct, a1, \
                        a2_element, a2_range=[a2],
                        a1_range=range(0,napm*(a2//napm)) + \
                                range(napm*(a2//napm+1),len(struct.geometry)))
        else:
            if not allow_all_pairs:
                if a2 < napm:
                    # Immediately returns because atoms must be from different
                    # molecules. However, if a1 must be selected from the first
                    # molecule since allow_all_pairs=False, and a2 is from the
                    # first molecule as well, then no result can be returned    
                    return []
                else:
                    return [min_distance(all_interatomic_distances(struct, a1, \
                            a2_element, a1_range=range(0,napm), a2_range=[a2]))]
            else:
                result = []
                for i in range (0,a2//napm)+range(a2//napm+1,nmpc):
                    result.append(min_distance(all_interatomic_distances(\
                    struct, a1, a2_element,
                    a1_range=range(i*napm,(i+1)*napm), a2_range=[a2])))
                return [x for x in result if x!=None]

    elif isinstance(a1,str) and isinstance(a2,str):
        if not block_same_pair and not allow_all_pairs:
            return all_interatomic_distances(struct, a1, a2, \
            a1_range=range(0,napm), a2_range=range(napm,len(struct.geometry)))
        elif not block_same_pair and allow_all_pairs:
            return all_interatomic_distances(struct, a1, a2,\
            a1_range=range(0,napc), a2_range=range(0,napc), napm=napm)
        elif block_same_pair and not allow_all_pairs:
            result = []
            for i in range (1,nmpc):
                result.append(min_distance(all_interatomic_distances(struct, a1, a2,\
                a1_range=range(0,napm), a2_range=range(i*napm,(i+1)*napm))))
            return [x for x in result if x!=None]
        elif a1 == a2:
            result = []
            for i in range (0, nmpc-1):
                for j in range (i,nmpc):
                    result.append(min_distance(all_interatomic_distances(struct,\
                    a1, a2, a1_range=range(i*napm,(i+1)*napm), a2_range=range(j*napm,(j+1)*napm))))
            return [x for x in result if x!=None]
        else:
            result = []
            for i in range (0,nmpc):
                for j in range(0,i)+range(i+1,nmpc):
                    result.append(min_distance(all_interatomic_distances(struct,\
                    a1, a2, a1_range=range(i*napm,(i+1)*napm), a2_range=range(j*napm,(j+1)*napm))))

            return [x for x in result if x!=None]
    
    else:
        raise ValueError("Invalid input of a1 and a2")

    

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

    
                
            

def min_distance(distance_list):
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

