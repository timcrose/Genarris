'''

Created by Patrick Kilecdi (Xiayue Li) on 11/25/2015

This program contains different modules that calls upon shared.py to conduct structural generation

'''
from generation import shared
from utilities import misc, write_log,parallel_run
from core import pool_management
import multiprocessing, os, shutil
import copy

def structure_generation_single(inst,molecule=None):
    '''
    A module for testing purposes
    Generate and output one single structure
    '''
    info_level = inst.get_eval("Genarris_master","info_level")

    struct = shared.single_structure_generation(inst,molecule)

    if struct == False:
        if info_level >= 2 or (info_level >= 1 and inst.has_procedure("Structure_Generation_Single")):
            write_log.write_master_log(inst,"structure_generation_single failed")
        return False

    if inst.has_option("structure_generation_single","struct_id"):
        struct.struct_id = inst.get("structure_generation_single","struct_id")
    elif inst.has_option("structure_generation_single","struct_id_scheme"):
        misc.struct_id_assign(struct, inst.get_eval("structure_generation_single","struct_id_scheme"), inst.get_with_default("structure_generation_single","structure_index","0"))

    if inst.has_option("structure_generation_single","output_file") or \
    inst.has_option("structure_generation_single","output_folder"):
        misc.dump_structure_with_inst(inst,"structure_generation_single",struct)
    return struct

def single_structure_generation_call(info):
    '''
    Call the single structure generation with inst and molecul
    info[0] should be inst, and info[1] should be the molecule
    Used so that multiprocessing.map can send multiple arguments
    '''
    return structure_generation_single(info[0],info[1])


def structure_generation_batch(inst):
    '''
    Generates a batch of structure
    '''
    info_level = inst.get_eval("Genarris_master","info_level")
    number_of_structures = inst.get_eval("structure_generation_batch","number_of_structures")
    inst.set_default("structure_generation_batch","number_of_attempts",number_of_structures)
    number_of_attempts = inst.get_eval("structure_generation_batch","number_of_attempts")

    if info_level >= 1 and not inst.get_boolean("structure_generation_batch","im_not_master_process"):
        write_log.write_master_log(inst,"generation.generation_modules.structure_generation_batch is called to generate %i structure; with " % (number_of_structures), additional_item = [["structure_generation","molecule_path"],["structure_generation","unit_cell_volume"],["structure_generation","space_group"],["structure_generation","wyckoff_list"]])

    if info_level>=2 and inst.get_boolean("structure_generation_batch","im_not_master_process"):
        write_log.write_master_log(inst,"structure_generation_batch process reporting from %s" % inst.get_master_node())


    inst.set_default("structure_generation_batch","tmp_dir",os.path.join(inst.get_master_working_dir(),"gen_tmp_"+misc.get_random_index()))
    tmp_dir = inst.get("structure_generation_batch","tmp_dir")
    structure_index_dir = os.path.join(tmp_dir,"structure_index.info")
    attempt_counter_dir = os.path.join(tmp_dir,"attempt_counter.info")

    if not inst.get_boolean("structure_generation_batch","im_not_master_process"):      
    #Generate a file that indicates the next structure inde
    #This will be read by the other multiprocessing instances
        misc.safe_make_dir(tmp_dir)
        if not os.path.isfile(structure_index_dir) or misc.retrieve_integer_and_subtract(structure_index_dir,0)<0:
            f = open(structure_index_dir,"w")
            f.write(str(number_of_structures-1))
            f.close()

        if not os.path.isfile(attempt_counter_dir) or misc.retrieve_integer_and_subtract(attempt_counter_dir,0)<0:
            f = open(attempt_counter_dir,"w")
            f.write(str(number_of_attempts-1))
            f.close()

    
    if inst.has_option("structure_generation_batch","processes_limit"):
        processes_limit = inst.get_eval("structure_generation_batch","processes_limit")
    else:
        processes_limit = inst.get_eval("Genarris_master","processes_limit")

    spread_processes_enabled = inst.get_boolean("structure_generation_batch","spread_processes_across_nodes")
    if spread_processes_enabled:
        nodelist = inst.get_nodelist()
    if spread_processes_enabled and len(nodelist)>1:
        #Launch similar python modules through ssh
        nodes_and_insts = []
        sample_inst = copy.deepcopy(inst)
        sample_inst.set("Genarris_master","procedures",str(["Structure_Generation_Batch"]))
        sample_inst.set("structure_generation_batch","processes_limit",str(int(processes_limit/len(nodelist))))
        processes_limit = int(processes_limit/len(nodelist)) #This will influence the rest of the running on this node
        sample_inst.set("structure_generation_batch","im_not_master_process","TRUE")
        sample_inst.set("structure_generation_batch","tmp_dir",tmp_dir)
        sample_inst.remove_option("structure_generation_batch","spread_processes_across_nodes")
        inst.remove_option("structure_generation_batch","spread_processes_across_nodes")
        nodelist.remove(inst.get_master_node())
        for node in nodelist:
            new_inst = copy.deepcopy(sample_inst)
            new_inst.set("Genarris_master","master_node",node)
            nodes_and_insts.append([node,new_inst])
        subp = parallel_run.ssh_spawn_python_subprocess(nodes_and_insts,inst.get_keywords([["parallel_settings","modules_launch",[]]],eval=True)[0])    
        if info_level >= 1:
            write_log.write_master_log(inst,"structure_generation_batch ssh spawned additional python subprocesses at "+str(nodelist))
    
    if processes_limit > 1:
        new_inst = copy.deepcopy(inst)
        new_inst.set("structure_generation_batch","processes_limit","1")
        new_inst.set("structure_generation_batch","im_not_master_process","TRUE")
        inst_list = []
        for i in range (processes_limit-1):
            inst_list.append(copy.deepcopy(new_inst))

        mulp = multiprocessing.Pool(processes=processes_limit-1)
        mulp_result = mulp.map_async(structure_generation_batch,inst_list)
        if info_level >= 1:
            write_log.write_master_log(inst,"structure_generation_batch multiprocessed additional %i processes at node %s " % (processes_limit-1,inst.get_master_node()))

    ################## Parallalization completed #####################      
    molecule = misc.load_structure_with_inst(inst,"structure_generation","molecule_format",inst.get("structure_generation","molecule_path")) 
    #Load the molecule so structure_generation_single does not have to load it every time
    
    attempt_counter = misc.retrieve_integer_and_subtract(attempt_counter_dir)
    structure_index = misc.retrieve_integer_and_subtract(structure_index_dir,subtract=0)
    sample_inst = copy.deepcopy(inst)
    sample_inst.remove_section("structure_generation_single")
    sample_inst.add_section("structure_generation_single")
    sample_inst.transfer_keywords("structure_generation_batch","structure_generation_single",["output_format","output_suffix"])

    result = []
    while attempt_counter >= 0 and structure_index >= 0 :
        
        struct = structure_generation_single(sample_inst,molecule)
        if struct == False: #Generation failure
            attempt_counter = misc.retrieve_integer_and_subtract(attempt_counter_dir)
            structure_index = misc.retrieve_integer_and_subtract(structure_index_dir,subtract=0)
            continue
    
        #Post-process generated structure   
        structure_index = misc.retrieve_integer_and_subtract(structure_index_dir)
        if structure_index < 0:
            break
        struct_id_scheme = inst.get_with_default("structure_generation_batch","struct_id_scheme",["struct_index","_","random_index"],eval=True)
        misc.struct_id_assign(struct,struct_id_scheme,str(structure_index).zfill(len(str(number_of_structures))))
        

        if inst.has_option("structure_generation_batch","output_folder"):
            output_folder = inst.get("structure_generation_batch","output_folder")
            if inst.has_option("structure_generation_batch","output_folder_split"):
                output_folder_split = inst.get_eval("structure_generation_batch","output_folder_split")
                if number_of_structures%output_folder_split==0:
                    folder_index = structure_index/(number_of_structures/output_folder_split)
                else:
                    folder_index = structure_index/(number_of_structures/output_folder_split+1)
                output_folder = os.path.join(output_folder,str(folder_index))
            sample_inst.set("structure_generation_single","output_folder",output_folder)
            misc.dump_structure_with_inst(sample_inst,"structure_generation_single",struct)
            sample_inst.remove_option("structure_generation_single","output_folder")
            
        attempt_counter = misc.retrieve_integer_and_subtract(attempt_counter_dir)

    ########################## Collecting results #########3###############
    coll = [x for x in result if x!=False]
    if processes_limit > 1:
        if info_level >= 1:
            write_log.write_master_log(inst,"structure_generation_batch beginning to join multiprocesses")
        mulp_result = mulp_result.get()
        mulp.close()
        for l in mulp_result:
            coll += l
        if info_level >= 1:
            write_log.write_master_log(inst,"structure_generation_batch multiprocesses joined")

    if spread_processes_enabled and (len(nodelist)>1 or nodelist[0]!=inst.get_master_node()):
        if info_level >= 1:
            write_log.write_master_log(inst,"structure_generation_batch beginning to join ssh subprocesses")
        for sp in subp:
            sp.wait()
        if info_level >= 1:
            write_log.write_master_log(inst,"structure_generation_batch ssh subprocesses joined")

    structure_index = misc.retrieve_integer_and_subtract(structure_index_dir,0)

    if not inst.get_boolean("structure_generation_batch","im_not_master_process"):
        if misc.retrieve_integer_and_subtract(structure_index_dir,subtract=0)<0:
            shutil.rmtree(tmp_dir)
#       if info_level >= 1:
#           write_log.write_master_log(inst,"structure_generation_batch master processes collected %i structures" % len(coll))


    return coll, number_of_structures-max(-1,structure_index)-1

def estimate_unit_cell_volume(inst):
    '''
    Module for estimating the unit cell of a given molecule
    Using binary search and a fixed generating success ratio to estimate experimental structure's volume
    The ordinary parameter for each generation should be already specified in structure_generation section of the instruction file
    '''
    sname = "estimate_unit_cell_volume"
    volume_lower = inst.get_eval(sname,"volume_lower_limit")
    volume_upper = inst.get_eval(sname,"volume_upper_limit")
    iters = inst.get_eval(sname,"binary_search_iterations")
    trials = inst.get_with_default(sname,"attempts_per_iteration",500,eval=True)
    successes = inst.get_with_default(sname,"desired_success_per_iteration",50,eval=True)
    stol = inst.get_with_default(sname,"desired_success_tolerance",7,eval=True)

    info_level = inst.get_info_level()

    inst = copy.deepcopy(inst)
    inst.remove_section("structure_generation_batch")
    inst.add_section("structure_generation_batch")
    inst.transfer_keywords(sname,"structure_generation_batch",
                   ["output_folder_split",
                "output_id_scheme",
                "output_format",
                "processes_limit",
                "spread_processes_across_nodes"])
    inst.set("structure_generation_batch","number_of_attempts",str(trials))
    inst.set("structure_generation_batch","number_of_structures",str(trials))
    inst.set_info_level(0) #Mute output from other module
    
    wml = write_log.write_master_log
    if info_level >= 1:
        wml(inst,"Beginning unit cell volume estimate procedure")
        wml(inst,"Lower volume limit: %f, upper volume limit: %f" 
              % (volume_lower, volume_upper))
        wml(inst,"Number of mid point iterations: %i" % iters)
        wml(inst,"Number of trials per iteration: %i; required successes: %i"
                % (trials, successes))

    original_inst = copy.deepcopy(inst)
    for i in range (iters):
        inst = copy.deepcopy(original_inst)
        volume = (volume_lower + volume_upper) / 2
        inst.set("structure_generation","unit_cell_volume",str(volume))
        inst.set("structure_generation_batch",
             "tmp_dir",
             "step_%i_%s" % (i,misc.get_random_index()))

        if info_level>=1: 
            wml(inst,"Iteration %i: volume examined is %f" % (i,volume))
            wml(inst,"Temporary folder: " + inst.get("structure_generation_batch",
                                 "tmp_dir"))
    
        if inst.has_option(sname,"output_folder"):
            output_folder = os.path.join(inst.get(sname,"output_folder"),
                             str(i))
            inst.set("structure_generation_batch",
                 "output_folder",
                 output_folder)
        
        coll, counts = structure_generation_batch(inst)
        if info_level>=1:
            wml(inst,"Number of successful generations: " + str(counts))
        if abs(counts-successes) <= stol:
            if info_level>=1: 
                wml(inst,"Number of successes fall into the desired range")
                wml(inst,"Final volume result: " + str(volume))
            return volume

        if counts > successes:
            if info_level >= 1:
                wml(inst,"Too many successes; decreasing upper limit")
            volume_upper = volume
        else:
            if info_level >= 1:
                wml(inst,"Too few successes; increasing lower limit")
            volume_lower = volume

    volume = (volume_upper+volume_lower) / 2
    if info_level >= 1:
        wml(inst,"Final volume result: "+ str(volume))
        wml(inst,"WARNING: The final success rate has not converged to desired range")
        wml(inst,"Consider increasing iterations or success tolerance")


