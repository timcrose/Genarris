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
from Genarris.core import instruct
#from utilities import parallel_run, write_log
from Genarris.utilities import write_log, mpi_utils
from mpi4py import MPI
import time, random

import json
from glob import glob

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
    Master class of Genarris. It controls all aspects of the 
    Genarris workflow which can be executed individually or sequantially. 
    Begins by reading and intepreting the configuration file. 
    Calls the defined procedures with the options specified in the 
    configuration file. Some options may be inferred from previous sections 
    if they are not present in every section. 
    
    '''
    def __init__(self,inst_path):
        '''
        Defines abstract workflow of Genarris. 
        
        '''
        start_time = time.time()
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
        procedures_with_master_slave = ['Run_FHI_Aims_Batch', 
            'Harris_Single_Molecule_Prep', 'Harris_Approximation_Batch', 
            'Relax_Single_Molecule', 'FHI_Aims_Energy_Evaluation']
        
        for section in self.inst.sections():
            print(section, flush=True)
            print(self.inst.options(section), flush=True)
            for option in self.inst.options(section):
                print(self.inst.get(section,option), flush=True)
            print(' ', flush=True)
        
        #for procedure in list of procedures, run that procedure
        active_comm = None
        split_comm = None
        for i, procedure in enumerate(procedures):
            if i > 0 and procedure[i - 1] == 'Pygenarris_Structure_Generation' and \
                    self.inst.get_with_default('pygenarris_structure_generation', 'omp_num_threads', 1, eval=True) > 1:
                mpi_utils.barrier(comm, tag=2277437, sleep=60)
            else:
                comm.barrier()
            if comm.rank == 0:
                data = self.inst
            else:
                data = None
            self.inst = comm.bcast(data, root=0)
            
            #free the active communicator: the communicator with ranks that 
            # execute a given procedure
            if active_comm is not None:
                try:
                    active_comm.Free()
                except:
                    pass
            if split_comm is not None:
                try:
                    split_comm.Free()
                except:
                    pass
            #See if splitting the communicator is necessary. It will be necessary
            # if num_cores != world_size
            #default num_cores is all available cores
            num_cores = self.inst.get_with_default(procedure.lower(), 'num_cores', comm.size, eval=True)
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
                    active_comm_rank = active_comm.Get_rank()
                except:
                    #rank doesn't belong to the set of active ranks
                    continue
                #Split into the replicas desired, if any
                try:
                    num_replicas = int(self.inst.get(procedure.lower(), 'num_replicas'))
                except:
                    num_replicas = 1
                if num_replicas == 1:
                    if procedure in procedures_with_master_slave:
                        if world_rank == 0:
                            color = MPI.UNDEFINED
                        else:
                            color = (active_comm.rank - 1) % num_replicas
                        split_comm = active_comm.Split(color)
                        getattr(self, procedure)(split_comm, comm, MPI.ANY_SOURCE, num_replicas)
                    else:
                        getattr(self, procedure)(active_comm)
                elif num_replicas > num_cores:
                    raise Exception('Cannot run more replicas than number of cores requested. ',
                                    'num_replicas = ', num_replicas, ' num_cores = ', num_cores)
                elif num_replicas == num_cores:
                    if procedure in procedures_with_master_slave:
                        raise Exception('Running Aims requires num_replicas to be less than num_cores '+
                                        'because one core is used for controlling the queue of jobs.')
                    else:
                        color = active_comm.rank % num_replicas
                        split_comm = active_comm.Split(color)
                        getattr(self, procedure)(split_comm)
                else: # 1 < num_replicas < num_cores
                    if procedure in procedures_with_master_slave:
                        if world_rank == 0:
                            color = MPI.UNDEFINED
                        else:
                            color = (active_comm.rank - 1) % num_replicas
                        split_comm = active_comm.Split(color)
                        getattr(self, procedure)(split_comm, comm, MPI.ANY_SOURCE, num_replicas)
                    else:
                        color = active_comm.rank % num_replicas
                        split_comm = active_comm.Split(color)
                        getattr(self, procedure)(split_comm)
            else: # num_cores == world_size
                try:
                    num_replicas = int(self.inst.get(procedure.lower(), 'num_replicas'))
                except:
                    num_replicas = 1
                if num_replicas == 1:
                    #If num_cores requested is all of then only need to split if replicas > 1
                    if procedure in procedures_with_master_slave:
                        if world_rank == 0:
                            color = MPI.UNDEFINED
                        else:
                            color = (comm.rank - 1) % num_replicas
                        split_comm = comm.Split(color)
                        getattr(self, procedure)(split_comm, comm, MPI.ANY_SOURCE, num_replicas)
                    else:
                        getattr(self, procedure)(comm)
                elif num_replicas > num_cores:
                    raise Exception('Cannot run more replicas than number of cores requested. ',
                                    'num_replicas = ', num_replicas, ' num_cores = ', num_cores)
                elif num_replicas == num_cores:
                    if procedure in procedures_with_master_slave:
                        raise Exception('Running Aims requires num_replicas to be less than num_cores '+
                                        'because one core is used for controlling the queue of jobs.')
                    else:
                        color = comm.rank % num_replicas
                        split_comm = comm.Split(color)
                        getattr(self, procedure)(split_comm)
                else: # 1 < num_replicas < num_cores
                    if procedure in procedures_with_master_slave:
                        if world_rank == 0:
                            color = MPI.UNDEFINED
                        else:
                            color = (comm.rank - 1) % num_replicas
                        split_comm = comm.Split(color)
                        getattr(self, procedure)(split_comm, comm, MPI.ANY_SOURCE, num_replicas)
                    else:
                        color = comm.rank % num_replicas
                        split_comm = comm.Split(color)
                        getattr(self, procedure)(split_comm)
            

        comm.barrier()
        end_time = time.time()
        if world_rank == 0:
            print('num_cores', int(self.inst.get(procedure.lower(), 'num_cores')), end_time - start_time, flush=True)
        

    def Affinity_Propagation_Analyze_Preference(self, comm):
        from evaluation import affinity
        affinity.affinity_propagation_analyze_preference(self.inst)

    def Affinity_Propagation_Distance_Matrix(self, comm):
        from evaluation import affinity
        affinity.affinity_propagation_distance_matrix(self.inst)

    def Affinity_Propagation_Fixed_Clusters(self, comm):
        """
        AP that explores the setting of preference in order to generate
        desired number of clusters.
        
        Configuration File Options
        --------------------------
        output_dir : str
            Path to the directory where the chosen structures will be stored. 
        preference_range : list
            List of two values as the [min, max] of the range of allowable 
            preference values.
        structure_dir : str, inferred
            Path to the directory of files to be used for the calculation. 
            Default is to infer this value from the previous section.
        dist_mat_input_file : str, inferred
            Path to the distance matrix output from the descriptor calculation.
            Default is to infer this value from the previous sections.
        output_format : str, optional
            Format the structure files should be saved as. Default is both.
        cluster_on_energy : bool, optional 
            Uses energy values to determine examplars. Structures with the 
            lowest energy values from each cluster are selected. 
            Default is False.
        plot_histograms : bool, optional
            If histogram plots should be created of the volume and space
            groups. Default is False.
        num_of_clusters : int or float, optional
            Float, must be less than 0. Selects a fraction of the structures. 
            Int, selects specific number of structures equal to int. 
            Default is 0.1.
        num_of_clusters_tolerance : int, optional
            Algorithm will stop if it has generated the number of clusters 
            within the number of desired clusters and this tolerance. 
            Default is 0.
        max_sampled_preferences : int, optional
            Maximum number of preference values to try. 
        output_without_success : bool, optional
            Whether to perform output procedures if the algorithm has reached
            the maximum number of sampled preferences without finding the 
            correct number of clusters. Default is False.
        affinity_type : list, optional
            List of [type of afinity, value] argument Scikit-Learn AP alogrithm. 
        affinity_matrix_path : str, optional
            Path to the affinity matrix to use for the AP algorithm.
            Default is ``affinity_matrix.dat``.
        damping : float, optional
            damping argument for Scikit-Learn AP algorithm. Default is 0.5.
        convergence_iter : int, optional
            convergence_iter argument for Scikit-Learn AP algorithm.
            Default is 15. 
        max_iter : int, optional
            max_iter argument for Scikit-Learn AP algorithm. Default is 1000.
        preference : int, optional
            preference argument for Scikit-Learn AP algorithm. Default is None.
        verbose_output : bool, optional
            verbose argument for Scikit-Learn AP algorithm. Default is False.
        property_key : str, optional
            Key which the AP cluster will be stored in the properties of 
            each structure object. Default is ``AP_cluster``. 
        output_file : str, optional
            Path where info about the AP alogrithm execution will be stored.
            Default is ``./AP_cluster.info``.
        exemplars_output_dir : str, optional
            If provided, will output the examplars of each cluster to this 
            folder. Default is None.
        exemplars_output_format : str, optional
            File format of structures to be output. Default is both.
        structure_suffix : str, optional
            Suffix to apply to structure files which are written. 
            Default is ``.json``.
            
        output_dir_2: str, inferred
            Code automatically looks for the option output_dir_2 if the 
            output directory already exists. This is how the code currently 
            identifies that AP is running for a second time. Default behavior
            is to not use this option if output_dir does not already exist.
        num_of_clusters_2: int or float, optional
            num_of_clusters for second clustering step. Default value is 0.1.
        output_file_2 : str, inferred
            Use if running AP algorithm twice, such as in the Robust workflow.
            Default is to use output_file.
        exemplars_output_dir_2 : str, inferred
            Exemplars output directory if second clustering step is used. 
            Default is to use exemplars_output_dir.
        cluster_on_energy_2 : str, inferred
            How to choose examplars for the second clustering step. Default 
            is to use cluster_on_energy value. 
        energy_name_2 : str, inferred
            Energy name to use for second clustering step. Default is to use 
            energy_name.
    
        """
        
        from evaluation import affinity
        start_time = time.time()
        affinity.affinity_propagation_fixed_clusters(self.inst, comm)
        if comm.rank == 0:
            print('time taken for Affinity_Propagation_Fixed_Clusters', time.time() - start_time, flush=True)

    def Affinity_Propagation_Fixed_Silhouette(self, comm):
        from evaluation import affinity
        affinity.affinity_propagation_fixed_silhouette(self.inst)

    def Cluster_Based_Selection(self, comm):
        from evaluation import selection
        selection.cluster_based_selection(self.inst)

    def Estimate_Unit_Cell_Volume(self, comm):
        """
        Performs volume estimation using a machine learned model train on the 
        CSD and based on Monte Carlo volume integration and topological 
        molecular fragments. See Genarris 2.0 paper for full description.
            
        Configuration File Options
        --------------------------
        volume_mean : float, optional
            If provided, uses this value as the volume generation mean without
            using the ML model to etimate the volume.
        volume_std : float, optional
            If provided, uses this value for structure generation, otherwise 
            a default value of 0.075 multiplied by the prediction volume per 
            unit cell is provided.
        
        """
        from utilities import volume_estimator
        start_time = time.time()
        self.inst = volume_estimator.estimate_unit_cell_volume(self.inst, comm)
        if comm.rank == 0:
            print('time taken for Estimate_Unit_Cell_Volume', 
                  time.time() - start_time, flush=True)

    def Estimate_Unit_Cell_Volume_v1(self, comm):
        from generation import generation_modules
        generation_modules.estimate_unit_cell_volume(self.inst)

    def FHI_Aims_Single_Run(self, comm):
        from evaluation import run_aims
        run_aims.fhi_aims_single_run(self.inst, comm)

    def FHI_Aims_Batch_Run(self, comm):
        from evaluation import FHI_aims
        FHI_aims.fhi_aims_batch_run(self.inst)

    def FHI_Aims_Extract(self, comm):
        from evaluation import fhi_aims_modules
        fhi_aims_modules.fhi_aims_extract(self.inst)

    def FHI_Aims_Energy_Evaluation(self, comm, world_comm, MPI_ANY_SOURCE, 
                                   num_replicas):
        """
        Runs Self-Consistent Field calculation on a pool of structures. 
        
        Configuration File Options
        --------------------------
        See :meth:`Run_FHI_Aims_Batch`
        
        """
        from evaluation import run_fhi_aims
        start_time = time.time()
        run_fhi_aims.run_fhi_aims_batch(comm, world_comm, MPI_ANY_SOURCE, num_replicas, inst=self.inst, sname='fhi_aims_energy_evaluation')
        if world_comm.rank == 0:
            print('time taken for FHI_Aims_Energy_Evaluation', time.time() - start_time, flush=True)

    def FHI_Aims_Scavenge(self, comm):
        from evaluation import FHI_aims
        FHI_aims.fhi_aims_scavenge(self.inst)

    def Find_Duplicates(self, comm):
        from evaluation import duplicate
        duplicate.find_duplicates(self.inst)

    def Find_Duplicates_Distance_Matrix(self, comm):
        from evaluation import duplicate
        duplicate.find_duplicates_distance_matrix(self.inst)

    def Harris_Approximation_Batch(self, comm, world_comm, MPI_ANY_SOURCE, num_replicas):
        from evaluation import harris_approximation
        start_time = time.time()
        harris_approximation.harris_approximation_batch(comm, world_comm, MPI_ANY_SOURCE, num_replicas, self.inst)
        if world_comm.rank == 0:
            print('time taken for Harris_Approximation_Batch', time.time() - start_time, flush=True)

    def Harris_Approximation_Single(self, comm):
        from evaluation import harris_approximation
        harris_approximation.harris_approximation_single(self.inst)

    def Harris_Single_Molecule_Prep(self, comm, world_comm, MPI_ANY_SOURCE, num_replicas):
        from evaluation import harris_approximation
        start_time = time.time()
        harris_approximation.harris_single_molecule_prep(comm, world_comm, MPI_ANY_SOURCE, self.inst)
        if world_comm.rank == 0:
            print('time taken for Harris_Single_Molecule_Prep', time.time() - start_time, flush=True)

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

    def Pygenarris_Structure_Generation(self, comm):
        """
        Uses the Genarris module written in C to perform structure generation. 
        This module enables generation on special positions. 
        
        Configuration File Options
        --------------------------
        molecule_path : str
            Path to the relaxed molecule geometry.
        output_format : str, optional, default="json"
            Determines the type of file which will be output for each 
            structure. Can be one of: json, geo, both. 
        output_dir : str
            Path to the directory which will contain all generated structures
            which pass the intermolecular distance checks.
        
        
        num_structures : int
            Target number of structures to generate.
        Z : int
            Number of molecules per cell to generate. 
        volume_mean : float, optional
            See :meth:`Estimate_Unit_Cell_Volume`
        volume_std : float, optional
            See :meth:`Estimate_Unit_Cell_Volume`
        sr : float, optional
            Defines the minimum intermolecular distance that is considered
            physical by multiplying the sum of the van der Waals radii of the
            interacting atoms by sr. Default value is 0.85. 
        tol : float, optional
            Tolerance to be used to identify space groups compatible with the 
            input molecule.
        
        max_attempts_per_spg_per_rank : int
            Defines the maximum number of attempts the structure generator 
            makes before moving on to the next space group.
        num_structures_per_allowed_SG_per_rank : int
            Number of structures per space group per rank which will be 
            generated by Pygenarris.
        geometry_out_filename : str
            Filename where all structures generated by Pygenarris will be found.
        
        omp_num_threads : int
            Number of OpenMP threads to pass into Pygenarris
        truncate_to_num_structures : bool
            If true, will reduce pool to exactly the number defined by 
            num_structures. 
        
        
        """
        from generation import run_pygenarris
        start_time = time.time()
        run_pygenarris.pygenarris_structure_generation(inst=self.inst, comm=comm)
        if comm.rank == 0:
            print('time taken for Pygenarris_Structure_Generation', time.time() - start_time, flush=True)

    def Random_Value_Assignment(self, comm):
        from evaluation.pool_analysis import random_value_assignment
        random_value_assignment(self.inst)

    def RCD_Calculation(self, comm):
        from evaluation import relative_coordinate_descriptor
        start_time = time.time()
        relative_coordinate_descriptor.rcd_calculation(self.inst, comm)
        if comm.rank == 0:
            print('time taken for RCD_Calculation', time.time() - start_time, flush=True)

    def RCD_Difference_Compare_Single(self, comm):
        from evaluation import relative_coordinate_descriptor as rcd
        rcd.rcd_difference_compare_single(self.inst)

    def RCD_Difference_Folder(self, comm):
        from evaluation import relative_coordinate_descriptor as rcd
        rcd.rcd_difference_folder(self.inst)

    def RCD_Difference_Folder_Inner(self, comm):
        from evaluation import relative_coordinate_descriptor as rcd
        start_time = time.time()
        rcd.rcd_difference_folder_inner(self.inst, comm)
        if comm.rank == 0:
            print('time taken for RCD_Difference_Folder_Inner', time.time() - start_time, flush=True)

    def RCD_Difference_Calculation(self, comm):
        from evaluation import relative_coordinate_descriptor
        relative_coordinate_descriptor.rcd_difference_calculation(self.inst)

    def Run_Rdf_Calc(self, comm):
        """
        Runs RDF calculation for the pool of generated structures. RDF 
        descriptor is similar to that described in Behler and Parrinello 2007.
        Then calculates the structure difference matrix.        
        
        Configuration File Options
        --------------------------
        structure_dir : str, inferred
            Path to the directory of structures to evaluate. 
        dist_mat_fpath : str
            Path to file to write distance matrix to.
        output_dir : str
            Path of directory to write structures to (will create if it DNE). 
            If 'no_new_output_dir' then input structures will be overwritten.
        normalize_rdf_vectors: bool,optional
            Whether to normalize the rdf vectors over the columns of the 
            feature matrix before using them to compute the distance matrix. 
            Default is Falase. 
        standardize_distance_matrix: bool
            If True, standardizes the distance matrix. The method is to divide 
            all elements by the max value in the distance matrix. 
            Because it is a distance matrix and thus all elements are positive,
            the standardized elements will be in the range [0, 1].
            Default is False. 
        save_envs: bool, optional
            Whether to save the environment vectors calculated by the RDF 
            method in the output structure files. Default is False.
        cutoff : float, optional
            Cutoff radius to apply to the atom centered symmetry function.
            Default is 12.
        n_D_inter : int, optional
            Number of dimensions to use for each type of pair-wise 
            interatomic interaction found in the structure. Default is 12.
        init_scheme : str, optional
            Can be centered or shifted, as described in Gastegger et al. 2018.
            Default is shifted.
        eta_range : list, optional
            List of two floats which define the range for eta parameter in 
            Gastegger et al. 2018. Default is [0.05,0.5].
        Rs_range : list, optional
            List of two floats which define the range for Rs parameter in 
            Gastegger et al. 2018. Default is [[0.1,12].
        pdist_distance_type : str,optional
            Input parameter for the pdist function. Default is Euclidean. 
        
        """
        from evaluation import rdf_calc
        start_time = time.time()
        rdf_calc.run_rdf_calc(self.inst, comm)
        if comm.rank == 0:
            print('time taken for Run_Rdf_Calc', time.time() - start_time, flush=True)

    def RDF_Descriptor_By_Bin(self, comm):
        from evaluation import radial_distribution_function
        radial_distribution_function.rdf_descriptor_by_bin(self.inst)

    def RDF_Descriptor_By_Point(self, comm):
        from evaluation import radial_distribution_function
        radial_distribution_function.rdf_descriptor_by_point(self.inst)

    def Relax_Single_Molecule(self, comm, world_comm, MPI_ANY_SOURCE, num_replicas):
        """
        Calls run_fhi_aims_batch using the provided single molecule path. 
        
        Configuration File Options
        --------------------------
        See :meth:`Run_FHI_Aims_Batch`
        
        """
        from evaluation import run_fhi_aims
        start_time = time.time()
        run_fhi_aims.run_fhi_aims_batch(comm, world_comm, MPI_ANY_SOURCE, num_replicas, inst=self.inst, sname='relax_single_molecule')
        if world_comm.rank == 0:
            print('time taken for Relax_Single_Molecule', time.time() - start_time, flush=True)

    def Reverse_Harris_Approximation(self, comm):
        from evaluation import harris_approximation
        harris_approximation.reverse_harris_approximation(self.inst)
        
    def Reverse_Harris_Approximation_Batch(self, comm):
        from evaluation import harris_approximation
        harris_approximation.reverse_harris_approximation_batch(self.inst)

    def Run_FHI_Aims_Batch(self, comm, world_comm, MPI_ANY_SOURCE, num_replicas):
        """
        Runs FHI-aims calculations on a pool of structures using num_replicas.
        
        Configuration File Options
        --------------------------
        verbose : bool
            Controls verbosity of output.
        energy_name : str
            Property name which the calculated energy will be stored with in the 
            Structure file.
        output_dir : str
            Path to the directory where the output structure file will be saved.
        aims_output_dir : str
            Path where the aims calculation will take place.
        aims_lib_dir : str, inferred
            Path to the location of the directory containing the FHI-aims library 
            file. 
        molecule_path : str
            Path to the geometry.in file of the molecule to be calculated if 
            called using harris_single_molecule_prep or relax_single_molecule.
        structure_dir : str, inferred
            Path to the directory of structures to be calculated if calculation
            was called not using harris_single_molecule_prep or 
            relax_single_molecule.
        Z : int, inferred
            Number of molecules per cell. 
            
        """
        from evaluation import run_fhi_aims
        start_time = time.time()
        run_fhi_aims.run_fhi_aims_batch(comm, world_comm, MPI_ANY_SOURCE, num_replicas, inst=self.inst, sname='run_fhi_aims_batch')
        if world_comm.rank == 0:
            print('time taken for Run_FHI_Aims_Batch', time.time() - start_time, flush=True)

    def Specific_Radius_Batch(self, comm):
        from evaluation import pool_analysis
        pool_analysis.specific_radius_batch(self.inst)

    def Structure_Generation_Single(self, comm):
        from generation import generation_modules
        generation_modules.structure_generation_single(self.inst)
    
    def Structure_Generation_Batch(self, comm):
        from generation import generation_modules
        start_time = time.time()
        generation_modules.structure_generation_batch(self.inst, comm)
        if comm.rank == 0:
            print('time taken for Structure_Generation_Batch', time.time() - start_time, flush=True)

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

