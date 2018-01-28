'''

Created by Patrick Kilecdi (Xiayue Li) on 11/25/2015

This program contains different modules that calls upon shared.py to conduct structural generation

'''
from generation import shared
from utilities import misc, parallel_run
import multiprocessing, os, shutil, copy
from generation.generation_util import \
        get_structure_generator_from_inst
from utilities.parallel_master import get_parallel_run_args, \
        launch_parallel_run_single_inst
from utilities.write_log import print_time_log
from utilities.write_log import write_master_log as wml
from utilities.misc import safe_make_dir

def structure_generation_single(inst):
    '''
    A module for testing purposes
    Generate and output one single structure
    '''
    sname = "structure_generation_single"
    generator = get_structure_generator_from_inst(inst, sname)
    generator.generate_structure()

def structure_generation_batch(inst):
    '''
    Generates a batch of structure
    '''
    sname = "structure_generation_batch"
    print_time_log("Structure generation batch begins.")
    if inst.has_option(sname, "index_tracking_file"):
        print_time_log("Using existing index tracking file: " +
                inst.get(sname, "index_tracking_file"))
    else:
        _set_up_index_and_attempt_tracking_files(inst, sname)

    _make_path_absolute(inst, sname)
    _inst = copy.deepcopy(inst)
    _inst.set_procedures(["_Structure_Generation_Batch"])
    args = get_parallel_run_args(inst, sname)
    launch_parallel_run_single_inst(_inst, *args)
    _clean_up_index_and_attempt_tracking_files(inst, sname)

def _structure_generation_batch(inst):
    sname = "structure_generation_batch"
    generator = get_structure_generator_from_inst(inst, sname)
    generator.generate_structure_until_index_or_attempts_zero()

def _set_up_index_and_attempt_tracking_files(inst, sname):
    number_of_structures = inst.get_eval(sname, "number_of_structures")
    number_of_attempts = inst.get_or_none(sname, "number_of_attempts", None)

    tmp_dir = inst.get_tmp_dir(sname)
    safe_make_dir(tmp_dir)

    index_tracker = os.path.join(tmp_dir, "_structure_generation_index")
    f = open(index_tracker, "w")
    f.write(str(number_of_structures))
    f.close()
    inst.set(sname, "index_tracking_file", index_tracker)

    if number_of_attempts is None:
        return
    attempt_tracker = os.path.join(tmp_dir, "_structure_generation_attempt")
    f = open(os.path.join(tmp_dir, "_structure_generation_attempt"), "w")
    f.write(str(number_of_attempts))
    f.close()
    inst.set(sname, "attempt_tracking_file", attempt_tracker)

def _clean_up_index_and_attempt_tracking_files(inst, sname):
    tmp_dir = inst.get_tmp_dir(sname)
    index_tracker = os.path.join(tmp_dir, "_structure_generation_index")
    if os.path.isfile(index_tracker):
        os.remove(index_tracker)
    attempt_tracker = os.path.join(tmp_dir, "_structure_generation_attempt")
    if os.path.isfile(attempt_tracker):
        os.remove(attempt_tracker)

def _make_path_absolute(inst, sname):
    inst.make_path_absolute(sname, "molecule_path")
    inst.make_path_absolute(sname, "output_dir")
    inst.make_path_absolute(sname, "index_tracking_file")
    inst.make_path_absolute(sname, "attempt_tracking_file")

# TODO: Deprecate estimate_unit_cell_volume
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


