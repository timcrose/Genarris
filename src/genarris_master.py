"""
If any part of this module is used for a publication please cite:

X. Li, F. Curtis, T. Rose, C. Schober, A. Vazquez-Mayagoitia, K. Reuter,
H. Oberhofer, and N. Marom "Genarris: Random Generation of Molecular Crystal 
Structures and Fast Screening with a Harris Approximation, ",
J. Chem. Phys., DOI: 10.1063/1.5014038; arXiv 1803.02145 (2018)
"""
"""
Created on Wed Jun 17 22:13:09 2015

@author: Patrick Kilecdi
"""

import os
import sys,socket
from core import instruct
from utilities import parallel_run, write_log
from mpi4py import MPI
import time, random

start_time = time.time()

__author__ = "Xiayue Li, Timothy Rose, Christoph Schober, and Farren Curtis"
__copyright__ = "Copyright 2018, Carnegie Mellon University and "+\
                "Fritz-Haber-Institut der Max-Planck-Gessellschaft"
__credits__ = ["Xiayue Li", "Luca Ghiringhelli", "Farren Curtis", "Tim Rose",
               "Christoph Schober", "Alvaro Vazquez-Mayagoita",
               "Karsten Reuter", "Harald Oberhofer", "Noa Marom"]
__license__ = "BSD-3"
__version__ = "180324"
__maintainer__ = "Timothy Rose"
__email__ = "trose@andrew.cmu.edu"
__url__ = "http://www.noamarom.com"

def main():
    '''
    When calling genarris_master.py, should specifiy the path to an instruction file
    '''
    
    external_libs_dir = os.path.join(sys.argv[0],"external_libs")
    sys.path.append(external_libs_dir)
    if len(sys.argv)==1:
        test()
    else:
        main_process = Genarris(sys.argv[-1])

def test():
    print("You are running genarris_master.py under testing mode"
          " (without conf file). Please specify testing code in the "
          "test() function")

class Genarris():
    '''
    This is the master class of Genarris
    Takes the path to a configuration file as a necessary input
    '''
    def __init__(self,inst_path):
        '''
        Interprets the instruction and calls the respective attributes of self.
        '''
        comm = MPI.COMM_WORLD
        world_rank = comm.Get_rank()
        world_size = comm.Get_size()
        #inst_path is sys.argv[-1] aka path/to/ui.conf
        self.inst_path = inst_path
        #Instruct object inherits from SafeConfigParser
        self.inst = instruct.Instruct()
        
        #enable getting values from the conf file. In other words, this function
        # will read the config file in a format that SafeConfigParser can
        # parse.
        self.inst.load_instruction_from_file(inst_path)

        #sname is the section name. [Genarris_master] is required in any conf
        # file.
        sname = "Genarris_master"

        #Add section sname and add option 'working_dir' with a value that is
        # the full path of the conf file if any of these parts DNE. Then,
        # get the value of the 'working_dir' option.
        self.working_dir = self.inst.get_with_default(
            sname, "working_dir",
            os.path.abspath(os.path.dirname(inst_path)))

        #Fix the working directory to this folder
        #From this point on, unless absolute necessary
        #The working directory will stay there  

        os.chdir(self.working_dir)

        #Add section sname and add option 'tmp_dir' with a value if DNE
        self.inst.set_default(
            sname, "tmp_dir",
            os.path.join(self.working_dir, "tmp"))
        self.inst.set_default(sname, "master_log_path",
                os.path.join(self.working_dir,"Genarris.log"))
        self.inst.set_default(sname, "master_err_path",
                os.path.join(self.working_dir,"Genarris.err"))
          
        sys.stdout = open(self.inst.get(sname,"master_log_path"),"a")
        sys.stderr = open(self.inst.get(sname,"master_err_path"),"a")
        write_log.set_global_output(sys.stdout, sys.stderr)

        self.inst.set_default(sname,"script_path",
                      os.path.realpath(__file__))
        self.inst.set_default(sname,"master_node",
                      socket.gethostname())
        
        procedures = self.inst.get_keywords([[sname,"procedures"]],True)[0]
        
        for section in self.inst.sections():
            #add default num_cores if section is a procedure
            if section in procedures:
                self.inst.get_with_default(section, 'num_cores', 1, eval=True)
            #print(section)
            #print(self.inst.options(section))
            #for option in self.inst.options(section):
            #    print(self.inst.get(section,option))
            #print(' ')
        
        #for procedure in list of procedures, run that procedure
        active_comm = None
        for procedure in procedures:
            comm.barrier()
            
            #free the active communicator: the communicator with ranks that 
            # execute a given procedure
            if active_comm is not None:
                active_comm.Free()
            #See if splitting the communicator is necessary. It will be necessary
            # if num_cores != world_size
            num_cores = int(self.inst.get(procedure.lower(), 'num_cores'))
            if num_cores != world_size:
                if world_rank < num_cores:
                    #This rank will be used in the procedure
                    color = 0
                else:
                    color = MPI.UNDEFINED
                #get the communicator corresponding to the non-communicating
                # processes only
                active_comm = comm.Split(color)
                #Only run with however many ranks desired
                try:
                    #rank belongs to the set of active ranks
                    world_rank = active_comm.Get_rank()
                except:
                    #rank doesn't belong to the set of active ranks
                    continue
                getattr(self, procedure)(active_comm)
            else:
                #If num_cores requested is all of them then don't need to split
                getattr(self, procedure)(comm)
            

        
        comm.barrier()
        end_time = time.time()
        if world_rank == 0:
            print('num_cores', int(self.inst.get(procedure.lower(), 'num_cores')), end_time - start_time)
        

    def Affinity_Propagation_Analyze_Preference(self, comm):
        from evaluation import affinity
        affinity.affinity_propagation_analyze_preference(self.inst)

    def Affinity_Propagation_Distance_Matrix(self, comm):
        from evaluation import affinity
        affinity.affinity_propagation_distance_matrix(self.inst)

    def Affinity_Propagation_Fixed_Clusters(self, comm):
        from evaluation import affinity
        affinity.affinity_propagation_fixed_clusters(self.inst)

    def Affinity_Propagation_Fixed_Silhouette(self, comm):
        from evaluation import affinity
        affinity.affinity_propagation_fixed_silhouette(self.inst)

    def Cluster_Based_Selection(self, comm):
        from evaluation import selection
        selection.cluster_based_selection(self.inst)

    def Estimate_Unit_Cell_Volume(self, comm):
        from generation import generation_modules
        generation_modules.estimate_unit_cell_volume(self.inst)

    def FHI_Aims_Single_Run(self, comm):
        from evaluation import run_aims
        run_aims.aims_single_run(self.inst, comm)

    def FHI_Aims_Batch_Run(self, comm):
        from evaluation import FHI_aims
        FHI_aims.fhi_aims_batch_run(self.inst)

    def FHI_Aims_Extract(self, comm):
        from evaluation import fhi_aims_modules
        fhi_aims_modules.fhi_aims_extract(self.inst)

    def FHI_Aims_Scavenge(self, comm):
        from evaluation import FHI_aims
        FHI_aims.fhi_aims_scavenge(self.inst)

    def Find_Duplicates(self, comm):
        from evaluation import duplicate
        duplicate.find_duplicates(self.inst)

    def Find_Duplicates_Distance_Matrix(self, comm):
        from evaluation import duplicate
        duplicate.find_duplicates_distance_matrix(self.inst)

    def Harris_Approximation_Batch(self, comm):
        from evaluation import harris_approximation
        harris_approximation.harris_approximation_batch(self.inst)

    def Harris_Approximation_Single(self, comm):
        from evaluation import harris_approximation
        harris_approximation.harris_approximation_single(self.inst)

    def Harris_Single_Molecule_Prep(self, comm):
        from evaluation import harris_approximation
        harris_approximation.harris_single_molecule_prep(self.inst)

    def Interatomic_Distance_Evaluation(self, comm):
        from evaluation import interatomic_distance_evaluation
        interatomic_distance_evaluation.main(self.inst)

    def Interatomic_Proximities(self, comm):
        from evaluation import interatomic_proximities
        interatomic_proximities.interatomic_proximities(self.inst)

    def K_Mean_Clustering(self, comm):
        from evaluation import k_mean_clustering
        k_mean_clustering.k_mean_clustering(self.inst)

    def Niggli_Reduction_Batch(self, comm):
        from evaluation import pool_analysis
        pool_analysis.niggli_reduction_batch(self.inst)

    def Pool_Single_Structure_Analysis(self, comm):
        from evaluation import pool_analysis
        pool_analysis.pool_single_structure_analysis(self.inst)

    def Random_Value_Assignment(self, comm):
        from evaluation.pool_analysis import random_value_assignment
        random_value_assignment(self.inst)

    def RDF_Descriptor_By_Bin(self, comm):
        from evaluation import radial_distribution_function
        radial_distribution_function.rdf_descriptor_by_bin(self.inst)

    def RDF_Descriptor_By_Point(self, comm):
        from evaluation import radial_distribution_function
        radial_distribution_function.rdf_descriptor_by_point(self.inst)

    def RCD_Calculation(self, comm):
        from evaluation import relative_coordinate_descriptor
        relative_coordinate_descriptor.rcd_calculation(self.inst)

    def RCD_Difference_Compare_Single(self, comm):
        from evaluation import relative_coordinate_descriptor as rcd
        rcd.rcd_difference_compare_single(self.inst)

    def RCD_Difference_Folder(self, comm):
        from evaluation import relative_coordinate_descriptor as rcd
        rcd.rcd_difference_folder(self.inst)

    def RCD_Difference_Folder_Inner(self, comm):
        from evaluation import relative_coordinate_descriptor as rcd
        rcd.rcd_difference_folder_inner(self.inst)

    def RCD_Difference_Calculation(self, comm):
        from evaluation import relative_coordinate_descriptor
        relative_coordinate_descriptor.rcd_difference_calculation(self.inst)

    def Reverse_Harris_Approximation(self, comm):
        from evaluation import harris_approximation
        harris_approximation.reverse_harris_approximation(self.inst)
        
    def Reverse_Harris_Approximation_Batch(self, comm):
        from evaluation import harris_approximation
        harris_approximation.reverse_harris_approximation_batch(self.inst)

    def Specific_Radius_Batch(self, comm):
        from evaluation import pool_analysis
        pool_analysis.specific_radius_batch(self.inst)

    def Structure_Generation_Single(self, comm):
        from generation import generation_modules
        generation_modules.structure_generation_single(self.inst)
    
    def Structure_Generation_Batch(self, comm):
        from generation import generation_modules
        generation_modules.structure_generation_batch(self.inst)

    def _Structure_Generation_Batch(self):
        from generation import generation_modules
        generation_modules._structure_generation_batch(self.inst)

    def Vector_Distance_Calculation(self, comm):
        from evaluation import pool_analysis
        pool_analysis.vector_distance_calculation(self.inst)

    # Below are modules for testing
    def Test_Launch_Parallel_Run_Single_Inst(self):
        from utilities import util_test
        util_test.test_launch_parallel_run_single_inst_main(self.inst)
    def _Test_Launch_Parallel_Run_Single_Inst(self):
        from utilities import util_test
        util_test._test_launch_parallel_run_single_inst(self.inst)

if __name__ == "__main__":
    main()

