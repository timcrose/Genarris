"""
If any part of this module is used for a publication please cite:

X. Li, F. Curtis, T. Rose, C. Schober, A. Vazquez-Mayagoitia, K. Reuter,
H. Oberhofer, and N. Marom "Genarris: Random Generation of Molecular Crystal 
Structures and Fast Screening with a Harris Approximation, ",
J. Chem. Phys., DOI: 10.1063/1.5014038; arXiv 1803.02145 (2018)
"""
'''
These are test modules for utils
'''
from Genarris.core import instruct
from Genarris.utilities.write_log import print_time_log
from Genarris.utilities.parallel_master import launch_parallel_run_single_inst


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
