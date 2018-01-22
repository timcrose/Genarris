'''
This set of utilities functions contain the necessary functions that carry out parallelization of the code
'''
import os,socket,time
from utilities import misc
from utilities.write_log import print_time_log
import subprocess, math, inspect

def launch_parallel_run_single_inst(new_inst,
                                    parallel_run_command,
                                    parallel_run_instances,
                                    tmp_dir,
                                    clean_up):
    '''
    new_inst: A instruction file for the parallelized Genarris master run
    parallel_run_command: 
        A list of command arguments to be used for Subprocessing.
        Must contain $CONF to be replaced by 
    '''
    if not "$CONF" in parallel_run_command:
        raise ValueError(
            "parallel_run_command must contain '$CONF' to be replaced by"
            "temporary conf file path")
    misc.safe_make_dir(tmp_dir)
    print_time_log("Launching parallel instances of Genarris master; tmp_dir = " +
                   tmp_dir)
    genarris_master_path = \
            os.path.join(
                    os.path.dirname(
                        os.path.dirname(os.path.abspath(inspect.stack()[0][1]))),
                    "genarris_master.py")
    print_time_log("Auto-detected master path: " + genarris_master_path)

    processes = []
    conf_files = []
    log_files = []
    err_files = []
    for i in range(parallel_run_instances):
        instance_name = misc.get_random_index();
        instance_log = os.path.join(tmp_dir, instance_name + ".log")
        instance_err = os.path.join(tmp_dir, instance_name + ".err")
        new_inst.set("Genarris_master", "master_log_path", instance_log)
        new_inst.set("Genarris_master", "master_err_path", instance_err)
        log_files.append(instance_log)
        err_files.append(instance_err)

        inst_conf_file_path = os.path.join(tmp_dir, instance_name + ".conf")


        new_inst.write_to_file(inst_conf_file_path)
        conf_files.append(inst_conf_file_path)

        args = parallel_run_command[:]
        args[args.index("$CONF")] = inst_conf_file_path
        if "$MASTER" in args:
            args[args.index("$MASTER")] = genarris_master_path

        p = subprocess.Popen(args)
        processes.append(p)
        print_time_log("Launched instance %s using arguments %s"
                       % (instance_name, args))

    for p in processes:
        p.wait()

    print_time_log("Parallel instances of Genarris master have all completed")

    if clean_up:
        print_time_log("Cleaning up tmp files generated by the parallel launch")
        misc.safe_remove_files(conf_files)
        misc.safe_remove_files(log_files)
        misc.safe_remove_files(err_files)