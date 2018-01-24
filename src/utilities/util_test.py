'''
These are test modules for utils
'''
from core import instruct
from utilities.write_log import print_time_log
from utilities.parallel_master import launch_parallel_run_single_inst

def test_launch_parallel_run_single_inst_main(inst):
    sname = "test_launch_parallel_run_single_inst"
    # TODO: Update to use get_parallel_run_args util
    parallel_run_command = inst.get_with_default(
            sname, "parallel_run_command",
            ["python", "$MASTER", "$CONF"], eval=True)
    parallel_run_instances = inst.get_with_default(
            sname, "parallel_run_instances", 5, eval=True)
    tmp_dir = inst.get_tmp_dir()

    test_launch_parallel_run_single_inst(
        parallel_run_command,
        parallel_run_instances,
        tmp_dir)

def test_launch_parallel_run_single_inst(
    parallel_run_command, parallel_run_instances, tmp_dir):
    print_time_log("Launching test parallel instances of genarris master")
    
    new_inst = instruct.Instruct()
    new_inst.set("Genarris_master", "procedures",
                 str(["_Test_Launch_Parallel_Run_Single_Inst"]))

    launch_parallel_run_single_inst(
        new_inst, parallel_run_command, parallel_run_instances,
        tmp_dir, False)

    print_time_log("Completed test parallel instances of genarris master")


def _test_launch_parallel_run_single_inst(inst):
    '''
    Function called by the relaunched genarris_master instances
    '''
    print("Hello world!")
