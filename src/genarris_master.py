# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 22:13:09 2015

@author: Patrick Kilecdi
"""

import os
import sys,socket
#import core.file_handler as fh
from core import instruct
from utilities import parallel_run, write_log

#src_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),os.pardir))
#sys.path.append(src_dir)
#(options,argv)=fh.argument_opt()
#if options.generate: #generation of new structure is called

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
    from evaluation import relative_coordinate_descriptor as rcd
    rcd.test_create_ref_struct()
    

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

		self.working_dir = self.inst.get_with_default(sname,
		"working_dir",os.path.abspath(os.path.dirname(inst_path)))

		if not self.inst.has_option(sname,"working_dir"):
			self.inst.set(sname,"working_dir",self.working_dir)
		os.chdir(self.working_dir) 
		#Fix the working directory to this folder
		#From this point on, unless absolute necessary
		#The working directory will stay there
		
		if not self.inst.has_option(sname,"master_log_path"):
			self.inst.set(sname,"master_log_path",
			os.path.join(self.working_dir,"Genarris.log"))

		if not self.inst.has_option(sname,"master_err_path"):
			self.inst.set(sname,"master_err_path",
			os.path.join(self.working_dir,"Genarris.err"))
		
		sys.stdout = open(self.inst.get(sname,"master_log_path"),"a")
		sys.stderr = open(self.inst.get(sname,"master_err_path"),"a")

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
			self.inst.set_default("parllel_settings","nodes",str(nodelist))
			write_log.write_master_log("Genarris launched on nodes "+str(nodelist))

		except:
			pass
		self.inst.set_default("parallel_settings","processes_per_node","20")
		

		procedures = self.inst.get_keywords([[sname,"procedures"]],True)[0]
		for procedure in procedures:
			getattr(self,procedure)()

	def Estimate_Unit_Cell_Volume(self):
		from generation import generation_modules
		generation_modules.estimate_unit_cell_volume(self.inst)

	def Harris_Approximation_Batch(self):
#		os.system("module load gcc/4.8.2")
#		os.system("module load anaconda")
		from evaluation import harris_approximation
		harris_approximation.harris_approximation_batch(self.inst)

	def Harris_Approximation_Single(self):
		from evaluation import harris_approximation
		harris_approximation.harris_approximation_single(self.inst)

	def Harris_Single_Molecule_Prep(self):
		from evaluation import harris_approximation
		harris_approximation.harris_single_molecule_prep(self.inst)

	def Interatomic_Distance_Evaluation(self):
		from evaluation import pool_analysis
		pool_analysis.interatomic_distance_evaluation(self.inst)

        def Interatomic_Proximities(self):
                from evaluation import interatomic_proximities
                interatomic_proximities.interatomic_proximities(self.inst)

	def Aims_Single_Run(self):
		from evaluation import run_aims
		run_aims.aims_single_run(self.inst)

        def FHI_Aims_Batch_Run(self):
                from evaluation import FHI_aims
                FHI_aims.fhi_aims_batch_run(self.inst)

        def FHI_Aims_Scavenge(self):
                from evaluation import FHI_aims
                FHI_aims.fhi_aims_scavenge(self.inst)

	def K_Mean_Clustering(self):
		from evaluation import k_mean_clustering
		k_mean_clustering.k_mean_clustering(self.inst)

        def Affinity_Propagation_Distance_Matrix(self):
                from evaluation import affinity
                affinity.affinity_propagation_distance_matrix(self.inst)

        def Affinity_Propagation_Analyze_Preference(self):
                from evaluation import affinity
                affinity.affinity_propagation_analyze_preference(self.inst)

        def Affinity_Propagation_Fixed_Clusters(self):
                from evaluation import affinity
                affinity.affinity_propagation_fixed_clusters(self.inst)

        def Cluster_Based_Selection(self):
                from evaluation import selection
                selection.cluster_based_selection(self.inst)

        def Find_Duplicates(self):
                from evaluation import duplicate
                duplicate.find_duplicates(self.inst)

        def Find_Duplicates_Distance_Matrix(self):
                from evaluation import duplicate
                duplicate.find_duplicates_distance_matrix(self.inst)

        def Niggli_Reduction_Batch(self):
                from evaluation import pool_analysis
                pool_analysis.niggli_reduction_batch(self.inst)

        def Specific_Radius_Batch(self):
                from evaluation import pool_analysis
                pool_analysis.specific_radius_batch(self.inst)

	def Pool_Single_Structure_Analysis(self):
		from evaluation import pool_analysis
		pool_analysis.pool_single_structure_analysis(self.inst)
		
	def Radial_Distribution_Function(self):
		from evaluation import pool_analysis
		pool_analysis.radial_distribution_function(self.inst)
	def Radial_Distribution_Function_1(self):
		from evaluation import pool_analysis
		pool_analysis.radial_distribution_function_1(self.inst)

        def RCD_Vector_Calculation(self):
                from evaluation import relative_coordinate_descriptor
                relative_coordinate_descriptor.rcd_vector_calculation(self.inst)

        def RCD_Difference_Calculation(self):
                from evaluation import relative_coordinate_descriptor
                relative_coordinate_descriptor.rcd_difference_calculation(self.inst)
        
        def RCD_Difference_Compare_Single(self):
                from evaluation import relative_coordinate_descriptor as rcd
                rcd.rcd_difference_compare_single(self.inst)

        def RCD_Difference_Folder(self):
                from evaluation import relative_coordinate_descriptor as rcd
                rcd.rcd_difference_folder(self.inst)

        def RCD_Difference_Folder_Inner(self):
                from evaluation import relative_coordinate_descriptor as rcd
                rcd.rcd_difference_folder_inner(self.inst)

               

	def Reverse_Harris_Approximation(self):
		from evaluation import harris_approximation
		harris_approximation.reverse_harris_approximation(self.inst)
	def Reverse_Harris_Approximation_Batch(self):
		from evaluation import harris_approximation
		harris_approximation.reverse_harris_approximation_batch(self.inst)

	def Structure_Generation_Single(self):
		from generation import generation_modules
		generation_modules.structure_generation_single(self.inst)
	
	def Structure_Generation_Batch(self):
		from generation import generation_modules
		generation_modules.structure_generation_batch(self.inst)

        def Vector_Distance_Calculation(self):
                from evaluation import pool_analysis
                pool_analysis.vector_distance_calculation(self.inst)
	

if __name__ == "__main__":
	main()	
