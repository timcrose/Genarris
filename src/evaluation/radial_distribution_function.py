import copy, math
import numpy as np
from utilities import write_log, misc
from core import structure_handling
from evaluation.evaluation_util import BatchSingleStructureOperation, \
        load_batch_single_structure_operation_keywords
from evaluation.interatomic_distance_evaluation import all_interatomic_distances

def rdf_descriptor_by_bin(inst):
    '''
    Returns a list of interatomic distances between pairs of certain species
    '''
    sname = "rdf_descriptor_by_bin" #Module name
    kwargs = load_batch_single_structure_operation_keywords(inst, sname)

    atomic_pairs = inst.get_eval(sname, "atomic_pairs")
    distance_range = inst.get_eval(sname, "distance_range")
    bins = inst.get_with_default(sname, "bins", 10, eval=True)
    normalized = inst.get_boolean(sname, "normalized")
    property_name = inst.get_with_default(sname, "property_name",
                                          "rdf_vector")
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

def _rdf_descriptor_by_bin(struct, atomic_pairs,
        distance_range, bins=10, normalized=False, property_name="rdf_vector"):
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

    rdf_vector = []
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

        rdf_vector += vector[:]

    original_struct.set_property(property_name, rdf_vector)
    return (original_struct, rdf_vector)

def _radial_distribution_function_callback(result, output_info_file):
    result.sort(key=lambda x : x[0].struct_id)
    message = ""
    for case in result:
        struct = case[0]
        rdf_vector = case[1]
        message += struct.struct_id + " " + " ".join(map(str, rdf_vector))
        message += "\n"
    write_log.write_log(output_info_file, message, time_stamp=False)
        
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


