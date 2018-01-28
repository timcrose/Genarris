# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 22:13:09 2015

@author: Patrick Kilecdi
"""

import os
import sys,socket
from core import instruct
from utilities import parallel_run, write_log

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
    Takes the path to a configuration file as an necessary input
    '''
    def __init__(self,inst_path):
        '''
        Interprets the instruction and calls the respective attributes of self
        '''
        self.inst_path = inst_path
        self.inst = instruct.Instruct()
        self.inst.load_instruction_from_file(inst_path)

        sname = "Genarris_master"
        if self.inst.has_section("ISGEP_master"):
        #Some examples may be still using the old name
            self.inst.transfer_keywords("ISGEP_master",
            sname,self.inst.options("ISGEP_master"))

        self.working_dir = self.inst.get_with_default(
            sname, "working_dir",
            os.path.abspath(os.path.dirname(inst_path)))

        #Fix the working directory to this folder
        #From this point on, unless absolute necessary
        #The working directory will stay there  
        if not self.inst.has_option(sname,"working_dir"):
            self.inst.set(sname,"working_dir",self.working_dir)
        os.chdir(self.working_dir)

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

        job_scheduler = self.inst.get_with_default("parallel_settings",
                               "job_scheduler",
                               "SLURM")

        self.inst.set_default("parallel_settings","nodelist",
            str(parallel_run.get_nodelist(job_scheduler)))

        try:
            nodelist = parallel_run.get_nodelist(self.inst.get("parallel_settings","job_scheduler"))
            self.inst.set_default("parallel_settings","nodes",str(nodelist))
            write_log.write_master_log("Genarris launched on nodes "+str(nodelist))

        except:
            pass
        self.inst.set_default("parallel_settings","processes_per_node","20")
        

        procedures = self.inst.get_keywords([[sname,"procedures"]],True)[0]
        for procedure in procedures:
            getattr(self,procedure)()

    def Affinity_Propagation_Analyze_Preference(self):
        from evaluation import affinity
        affinity.affinity_propagation_analyze_preference(self.inst)

    def Affinity_Propagation_Distance_Matrix(self):
        from evaluation import affinity
        affinity.affinity_propagation_distance_matrix(self.inst)

    def Affinity_Propagation_Fixed_Clusters(self):
        from evaluation import affinity
        affinity.affinity_propagation_fixed_clusters(self.inst)

    def Cluster_Based_Selection(self):
        from evaluation import selection
        selection.cluster_based_selection(self.inst)

    def Estimate_Unit_Cell_Volume(self):
        from generation import generation_modules
        generation_modules.estimate_unit_cell_volume(self.inst)

    def FHI_Aims_Single_Run(self):
        from evaluation import run_aims
        run_aims.aims_single_run(self.inst)

    def FHI_Aims_Batch_Run(self):
        from evaluation import FHI_aims
        FHI_aims.fhi_aims_batch_run(self.inst)

    def FHI_Aims_Extract(self):
        from evaluation import fhi_aims_modules
        fhi_aims_modules.fhi_aims_extract(self.inst)

    def FHI_Aims_Scavenge(self):
        from evaluation import FHI_aims
        FHI_aims.fhi_aims_scavenge(self.inst)

    def Find_Duplicates(self):
        from evaluation import duplicate
        duplicate.find_duplicates(self.inst)

    def Find_Duplicates_Distance_Matrix(self):
        from evaluation import duplicate
        duplicate.find_duplicates_distance_matrix(self.inst)

    def Harris_Approximation_Batch(self):
        from evaluation import harris_approximation
        harris_approximation.harris_approximation_batch(self.inst)

    def Harris_Approximation_Single(self):
        from evaluation import harris_approximation
        harris_approximation.harris_approximation_single(self.inst)

    def Harris_Single_Molecule_Prep(self):
        from evaluation import harris_approximation
        harris_approximation.harris_single_molecule_prep(self.inst)

    def Interatomic_Distance_Evaluation(self):
        from evaluation import interatomic_distance_evaluation
        interatomic_distance_evaluation.main(self.inst)

    def Interatomic_Proximities(self):
        from evaluation import interatomic_proximities
        interatomic_proximities.interatomic_proximities(self.inst)

    def K_Mean_Clustering(self):
        from evaluation import k_mean_clustering
        k_mean_clustering.k_mean_clustering(self.inst)

    def Niggli_Reduction_Batch(self):
        from evaluation import pool_analysis
        pool_analysis.niggli_reduction_batch(self.inst)

    def Pool_Single_Structure_Analysis(self):
        from evaluation import pool_analysis
        pool_analysis.pool_single_structure_analysis(self.inst)

    def RDF_Descriptor_By_Bin(self):
        from evaluation import radial_distribution_function
        radial_distribution_function.rdf_descriptor_by_bin(self.inst)

    def RDF_Descriptor_By_Point(self):
        from evaluation import radial_distribution_function
        radial_distribution_function.rdf_descriptor_by_point(self.inst)

    def RCD_Calculation(self):
        from evaluation import relative_coordinate_descriptor
        relative_coordinate_descriptor.rcd_calculation(self.inst)

    def RCD_Difference_Compare_Single(self):
        from evaluation import relative_coordinate_descriptor as rcd
        rcd.rcd_difference_compare_single(self.inst)

    def RCD_Difference_Folder(self):
        from evaluation import relative_coordinate_descriptor as rcd
        rcd.rcd_difference_folder(self.inst)

    def RCD_Difference_Folder_Inner(self):
        from evaluation import relative_coordinate_descriptor as rcd
        rcd.rcd_difference_folder_inner(self.inst)

    def RCD_Difference_Calculation(self):
        from evaluation import relative_coordinate_descriptor
        relative_coordinate_descriptor.rcd_difference_calculation(self.inst)

    def Reverse_Harris_Approximation(self):
        from evaluation import harris_approximation
        harris_approximation.reverse_harris_approximation(self.inst)
    def Reverse_Harris_Approximation_Batch(self):
        from evaluation import harris_approximation
        harris_approximation.reverse_harris_approximation_batch(self.inst)

    def Specific_Radius_Batch(self):
        from evaluation import pool_analysis
        pool_analysis.specific_radius_batch(self.inst)

    def Structure_Generation_Single(self):
        from generation import generation_modules
        generation_modules.structure_generation_single(self.inst)
    
    def Structure_Generation_Batch(self):
        from generation import generation_modules
        generation_modules.structure_generation_batch(self.inst)

    def _Structure_Generation_Batch(self):
        from generation import generation_modules
        generation_modules._structure_generation_batch(self.inst)

    def Vector_Distance_Calculation(self):
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
