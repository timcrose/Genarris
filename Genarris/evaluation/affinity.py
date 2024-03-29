"""
If any part of this module is used for a publication please cite:

X. Li, F. Curtis, T. Rose, C. Schober, A. Vazquez-Mayagoitia, K. Reuter,
H. Oberhofer, and N. Marom "Genarris: Random Generation of Molecular Crystal 
Structures and Fast Screening with a Harris Approximation, ",
J. Chem. Phys., DOI: 10.1063/1.5014038; arXiv 1803.02145 (2018)
"""
'''
Created by Patrick Kilecdi on 12/31/2016

Conducts Affinity Propagation for a given collection
'''

import os, random, psutil
from sklearn.cluster import AffinityPropagation
from sklearn.metrics import silhouette_score
import numpy as np
from Genarris.utilities.misc import output_pool
from Genarris.utilities import file_utils, time_utils, math_utils
from Genarris.core.structure import StoicDict, get_struct_coll
from bisect import bisect
from Genarris.core.instruct import get_last_active_procedure_name


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

iter_n = 0

class APHandler():
    def __init__(self, inst, sname, comm):
        self.comm = comm
        self.rank = comm.Get_rank()
        self.size = comm.Get_size()
        self.distance_matrix = None
        self.num_of_clusters = inst.get_with_default(sname,"num_of_clusters", 0.1, eval=True)
        self.preference_range = inst.get_eval(sname, "preference_range")
        self.num_of_clusters_tolerance = inst.get_with_default(
                sname, "num_of_clusters_tolerance", 0, eval=True)
        self.max_sampled_preferences = inst.get_with_default(
                sname, "max_sampled_preferences", 100, eval=True)
        self.output_without_success = inst.get_boolean(
                sname, "output_without_success")
        self.run_num = 1
        # Parameters for AP calculation
        self.affinity_type = inst.get_with_default(
                sname, "affinity_type", ["multiplicative",1000], eval=True)
        self.damping = inst.get_with_default(
                sname, "damping", 0.5, eval=True)
        self.convergence_iter = inst.get_with_default(
                sname, "convergence_iter", 15, eval=True)
        self.max_iter = inst.get_with_default(
                sname, "max_iter", 1000, eval=True)
        self.preference = inst.get_with_default(
                sname, 'preference', None, eval=True)
        self.property_key = inst.get_with_default(
                sname, "property_key", "AP_cluster")
        self.output_file = inst.get_with_default(
                sname, "output_file", "./AP_cluster.info")
        if os.path.exists(self.output_file):
            self.run_num = 2
            self.output_file = inst.get_with_default(
                sname, "num_of_clusters_2", self.output_file)
        self.exemplars_output_dir = inst.get_or_none(
                sname, "exemplars_output_dir")
        if os.path.exists(self.exemplars_output_dir):
            self.run_num = 2
            self.exemplars_output_dir = inst.get_with_default(
                sname, "exemplars_output_dir_2", self.exemplars_output_dir)
        self.exemplars_output_format = inst.get_with_default(
                sname, "exemplars_output_format", "both")
        self.verbose_output = inst.get_boolean(sname, "verbose")
        
        # Parameters common to pool operations
        
        self.structure_suffix = inst.get_with_default(sname, 'structure_suffix', '.json')
        self.output_dir = inst.get(sname, 'output_dir')
        if os.path.exists(self.output_dir):
            self.run_num = 2
            self.output_dir = inst.get(sname, 'output_dir_2')
            self.num_of_clusters = inst.get_with_default(sname,"num_of_clusters_2", 0.1, eval=True)

        self.output_format = inst.get_with_default(sname, 'output_format', 'both')

        if self.run_num == 1:
            self.cluster_on_energy = inst.get_boolean(sname, 'cluster_on_energy')
            if self.cluster_on_energy:
                sname_list = [sname, 'fhi_aims_energy_evaluation', 'harris_approximation_batch']
                self.energy_name = inst.get_inferred(sname, sname_list, ['energy_name'] * len(sname_list))
        elif self.run_num == 2:
            self.cluster_on_energy = inst.get_inferred(sname, [sname] * 2, ['cluster_on_energy_2', 'cluster_on_energy'], type_=bool)
            if self.cluster_on_energy:
                sname_list = [sname, sname, 'fhi_aims_energy_evaluation', 'harris_approximation_batch']
                self.energy_name = inst.get_inferred(sname, sname_list, ['energy_name_2'] + ['energy_name'] * (len(sname_list) - 1))
        #print('self.cluster_on_energy', self.cluster_on_energy, flush=True)
        #print('self.energy_name', self.energy_name, flush=True)
        self.dist_mat_input_file = inst.get_inferred(sname, [sname, 'run_rsf_calc', 'rcd_difference_folder_inner', 'rcd_calculation'], 
                                                    ['dist_mat_input_file', 'dist_mat_fpath', 'diff_matrix_output', 'diff_matrix_output'], type_='file')

        self.affinity_matrix_path = inst.get_with_default(sname, 'affinity_matrix_path', 'affinity_matrix.dat')
        if self.run_num == 1:
            last_section = get_last_active_procedure_name(inst, sname, iteration=0)
            sname_list = [sname, last_section, 'run_rsf_calc', 'rcd_difference_folder_inner', 'rcd_calculation']
            self.structure_dir = inst.get_inferred(sname, sname_list, ['structure_dir'] + (4 * ['output_dir']), type_='dir')
        elif self.run_num == 2:
            last_section = get_last_active_procedure_name(inst, sname, iteration=1)
            sname_list = [sname, last_section, 'run_rsf_calc', 'rcd_difference_folder_inner', 'rcd_calculation']
            self.structure_dir_for_ap1 = inst.get_inferred(sname, sname_list, ['structure_dir'] + (4 * ['output_dir']), type_='dir')
            sname_list = [sname, last_section, 'fhi_aims_energy_evaluation', last_section, last_section, 'run_rsf_calc', 'rcd_difference_folder_inner', 'rcd_calculation']
            self.structure_dir = inst.get_inferred(sname, sname_list, ['structure_dir', 'output_dir', 'output_dir', 'exemplars_output_dir'] + (4 * ['output_dir']), type_='dir')
            ext_pos = self.dist_mat_input_file.find('.')
            self.dist_mat_input_file_1 = self.dist_mat_input_file
            self.dist_mat_input_file = self.dist_mat_input_file[:ext_pos] + '1' + self.dist_mat_input_file[ext_pos:]
            ext_pos = self.affinity_matrix_path.find('.')
            self.affinity_matrix_path = self.affinity_matrix_path[:ext_pos] + '1' + self.affinity_matrix_path[ext_pos:]

        if self.num_of_clusters < 1:
                self.num_of_clusters = int(self.num_of_clusters * len(file_utils.glob(os.path.join(self.structure_dir, '*.json'))))


        self.plot_histograms = inst.get_boolean(sname, 'plot_histograms')
        if self.plot_histograms:
            self.prop = inst.get_with_default(sname, 'prop', 'unit_cell_volume')
            if self.run_num == 1:
                self.prop_figname = inst.get_with_default(sname, 'prop_figname', 'raw_pool_volume_histogram.pdf')
            else:
                self.prop_figname = inst.get_with_default(sname, 'prop_figname_2', 'raw_pool_volume_histogram.pdf')
            self.prop_xlabel = inst.get_with_default(sname, 'prop_xlabel', 'Structure Volume, $\AA^3$')
            self.prop_ylabel = inst.get_with_default(sname, 'prop_ylabel', 'Counts')
            self.prop_figure_size = inst.get_with_default(sname, 'prop_figure_size', (12,8), eval=True)
            self.prop_label_size = inst.get_with_default(sname, 'prop_label_size', 24, eval=True)
            self.prop_tick_size = inst.get_with_default(sname, 'prop_tick_size', 18, eval=True)
            self.prop_tick_width = inst.get_with_default(sname, 'prop_tick_width', 3, eval=True)
            self.prop_GAtor_IP = inst.get_boolean(sname, 'prop_GAtor_IP')

            self.pygenarris_outfile = inst.get_with_default(sname, 'pygenarris_outfile', 'outfile')
            if self.run_num == 1:
                self.spg_bar_chart_fname = inst.get_with_default(sname, 'spg_bar_chart_fname', 'raw_pool_spg_bar_chart.pdf')
            else:
                self.spg_bar_chart_fname = inst.get_with_default(sname, 'spg_bar_chart_fname_2', 'raw_pool_spg_bar_chart.pdf')
            self.spg_bar_width = inst.get_with_default(sname, 'spg_bar_width', 0.5, eval=True)
            self.spg_bar_xlabel = inst.get_with_default(sname, 'spg_bar_xlabel', 'Allowed space groups')
            self.spg_bar_ylabel = inst.get_with_default(sname, 'spg_bar_ylabel', 'Count')
            if self.run_num == 1:
                self.spg_bar_title = inst.get_or_none(sname, 'spg_bar_title')
            else:
                self.spg_bar_title = inst.get_or_none(sname, 'spg_bar_title_2')
            self.spg_bar_tick_rotation = inst.get_with_default(sname, 'spg_bar_tick_rotation', 'vertical')

        if comm.rank == 0:
            print(sname, 'using structure_dir:', self.structure_dir, flush=True)

    def get_affinity_matrix(self):
        self.affinity_matrix = np.memmap(self.affinity_matrix_path, dtype='float32', mode='r', shape=self.distance_matrix_shape)


    def get_affinity_and_distance_matrix(self):
        #print("dist_mat_input_file",self.dist_mat_input_file,flush=True)
        
        #distance_matrix_1= np.reshape(np.memmap(self.dist_mat_input_file_1, dtype='float32', mode='r'),)
        #distance_matrix_shape_1=distance_matrix_1.shape
        if self.distance_matrix is None:
            self.distance_matrix = np.memmap(self.dist_mat_input_file, dtype='float32', mode='r')
            self.distance_matrix.resize((int(np.sqrt(len(self.distance_matrix))), int(np.sqrt(len(self.distance_matrix)))))
        self.distance_matrix_shape = self.distance_matrix.shape
            #print("the shape of distance matrix:", self.distance_matrix.shape, flush=True)

        if len(self.distance_matrix) != len(self.distance_matrix[0]):
            raise ValueError("Distance matrix is not a square matrix: "
                    + self.dist_mat_input_file)

        self.get_affinity_matrix()
        

    def make_affinity_matrix(self):
        if self.verbose_output:
            print('self.dist_mat_input_file used to create affinity matrix', self.dist_mat_input_file, flush=True)
        distance_matrix = np.memmap(self.dist_mat_input_file, dtype='float32', mode='r')

        m = np.sqrt(len(distance_matrix))
        if m != int(m):
            raise ValueError("Distance matrix is not a square matrix: "
                + self.dist_mat_input_file)
        m = int(m)

        distance_matrix.resize((m,m))
        
        if len(distance_matrix) != len(distance_matrix[0]):
            raise ValueError("Distance matrix is not a square matrix: "
                    + self.dist_mat_input_file)

        if self.affinity_type[0]=="exponential":
            affinity_matrix = -np.exp(distance_matrix * self.affinity_type[1])
        elif self.affinity_type[0] == "power":
            affinity_matrix = -distance_matrix ** self.affinity_type[1]
        elif self.affinity_type[0] == 'multiplicative':
            affinity_matrix = -distance_matrix * self.affinity_type[1]
        else:
            raise Exception('Unsuppored affinity_type. Got:', self.affinity_type[0])

        self.distance_matrix_shape = distance_matrix.shape
        if self.verbose_output:
            print('self.affinity_matrix_path', self.affinity_matrix_path, flush=True)
        # Write affinity matrix
        fp = np.memmap(self.affinity_matrix_path, dtype='float32', mode='write', shape=affinity_matrix.shape)
        fp[:] = affinity_matrix[:]
        if self.verbose_output:
            print('writing affinity matrix...', flush=True)
        affinity_matrix_filesize = os.path.getsize(self.affinity_matrix_path)
        time_utils.sleep(5)
        while os.path.getsize(self.affinity_matrix_path) != affinity_matrix_filesize:
            affinity_matrix_filesize = os.path.getsize(self.affinity_matrix_path)
            time_utils.sleep(3)
        if self.verbose_output:
            print('affinity matrix written', flush=True)

    def get_new_pref_range(self, num_of_clusters_and_pref_list, prev_l, prev_u, iter_n):
        '''
        Determine if AP generated an acceptable number of clusters and if not,
        return the new preference range to try as well as the preference range 
        of the previous iteration.
        
        num_of_clusters_and_pref_list : array-like, shape (num_cores, 2)
        
            Matrix where each row has form [num_clusters, pref] where num_clusters
            is the number of clusters produced by AP when preference pref is used.
            The reason these are grouped is to reduce MPI communication latency.
            
        prev_l : float
        
            The lower-bound on the preference range used on the previous iteration.
            
        prev_u : float
        
            The upper-bound on the preference range used on the previous iteration.
            
        iter_n : int
        
            The iteration number.
        '''
        
        #Separate the matrix into two lists, one list of the number of clusters
        # produced by each core, and the other are the preference values used 
        # by each core. These lists are ordered by MPI rank.
        num_clusters_list, pref_list = zip(*num_of_clusters_and_pref_list)
        num_clusters_list = list(num_clusters_list)
        
        #Reverse these lists because that's how the following code was written.
        #You could rewrite this function to work with the forward direction if
        # you want instead.
        num_clusters_list.reverse()
        pref_list = list(pref_list)
        pref_list.reverse()
        
        for n, num_clusters in enumerate(num_clusters_list):
            #If you've found the desired number of clusters or want to output 
            # without success on the last iteration:
            if abs(num_clusters - self.num_of_clusters) <= self.num_of_clusters_tolerance or\
                    (self.output_without_success and iter_n == self.max_sampled_preferences - 1):
                success = True
                
                if self.output_without_success and iter_n == self.max_sampled_preferences - 1:
                    #n is the index of the original list. Thus the len... - ... - 1
                    #Output the result with the closest number of clusters to the 
                    # desired amount.
                    n = len(num_clusters_list) -  num_clusters_list.index(
                            min(num_clusters_list, key=lambda x:abs(x - self.num_of_clusters))) - 1
                else:
                    n = len(num_clusters_list) - n - 1
                
                #In the case of success, the last four return values are just
                # for formatting
                return success, pref_list[n], n, 0, 0, 0, 0
        
        success = False
        
        #Get the index of the value that the desired number of clusters would 
        # be in if the list num_clusters_list were to maintain being ordered.
        #Warning: for this to work properly, num_clusters_list should be 
        # ordered at this point. As preference decreases, the number of clusters
        # generated by AP should also decrease. However, this is not always
        # strictly true. In general, don't demand AP to give you exactly a certain
        # number of clusters.
        i = bisect(num_clusters_list, self.num_of_clusters)
        
        if i == 0:
            #self.num_of_clusters is lower than any value obtained. Set the 
            # preference bounds to be the previous lower bound and the lowest
            # bound on this iteration to supposedly find a better region.
            pref_l, pref_u = prev_l, pref_list[0]
        elif i == len(pref_list):
            #self.num_of_clusters is higher than any value obtained. Set the 
            # preference bounds to be the previous upper bound and the highest
            # bound on this iteration to supposedly find a better region.
            pref_l, pref_u = pref_list[-1], prev_u
        else:
            #self.num_of_clusters is in between values obtained. Set the 
            # preference bounds to be the preference values on either side of 
            # where self.num_of_clusters would be inserted.
            pref_l, pref_u = pref_list[i - 1], pref_list[i]

        if self.verbose_output:
                print('i', i, 'pref_l',pref_l, 'pref_u', pref_u, 'prev_l', prev_l, 'prev_u', prev_u, 'pref_list', pref_list, flush=True)
        if self.size != 1:
            prev_l, prev_u = min(pref_list), max(pref_list)
            resample_prob = 0.01
        else:
            resample_prob = 0.001
        if sorted(num_clusters_list) != num_clusters_list:
            resample_prob *= 2
        if (num_clusters_list == [num_clusters_list[0]] * len(num_clusters_list) and\
                    self.size != 1) or\
                 random.random() < resample_prob or\
                     pref_u - pref_l < 0.000001:
            pref_l -= 50.0 + random.random()
            pref_u += 50.0 + random.random()
            prev_l -= 50.0 + random.random()
            prev_u += 50.0 + random.random()
            if self.verbose_output:
                print('changing pref values to not get stuck', flush=True)

        if self.verbose_output:
            print('about to exit get_new_pref_range', 'pref_l',pref_l, 'pref_u', pref_u, 'prev_l', prev_l, 'prev_u', prev_u, flush=True)
        
        #In this case, the 0 values are just there for formatting
        return success, 0, 0, pref_l, pref_u, prev_l, prev_u 

    def get_allowed_ranks_per_node(self, allowed_number_of_ranks_per_node, rank_hostnames):
        '''
        Determine which ranks may do AP (True) and which should not (False).
        '''
        unique_hostnames = set(rank_hostnames)
        assignments = []
        assignments_dct = {hostname:0 for hostname in unique_hostnames}
        for hostname in rank_hostnames:
            if assignments_dct[hostname] > allowed_number_of_ranks_per_node:
                assignments.append(False)
            else:
                assignments.append(True)
                assignments_dct[hostname] = assignments_dct[hostname] + 1
        return assignments
        
    def get_allowed_number_of_ranks_per_node(self, available_mem):
        '''
        Each node might only be able to have a limited number of ranks
        run AP simultaneously because of memory accrued in the AffinityPropagation.fit()
        function.
        '''
        n = self.distance_matrix_shape[0]
        # Memory needed by the AffinityPropagation.fit() function determined empirically in terms
        # of the number of rows in the matrix (units kB)
        mem_needed = 0.351 * (n ** 2) + 9 * n + 63100
        return int(available_mem / mem_needed)
                            
    def run_fixed_num_of_clusters(self):
        '''
        Run the AP algorithm in parallel until either an acceptable number of 
        of clusters has been generated or the number of iterations asked for 
        has expired.
        '''
        #print('inside run_fixed_num_of_clusters', flush=True)
        
        if self.num_of_clusters > len(self.coll):
            raise ValueError("Cannot cluster pool into more clusters " +
                    "than number of structures in it. self.num_of_clusters = " + 
                    str(self.num_of_clusters) + ", len(self.coll) = " +
                    str(len(self.coll)))
        #print('self.num_of_clusters', self.num_of_clusters, flush=True)
        #print('len(self.coll)', len(self.coll), flush=True)
        if self.max_sampled_preferences < 1:
            raise ValueError("max_sampled_preferences must be >= 1")
        #print('self.max_sampled_preferences', self.max_sampled_preferences, flush=True)
        if self.rank == 0:
            print("Running Affinity Propagation with fixed number " +
                "of clusters " + str(self.num_of_clusters), flush=True)
        self.enable_subdir_output = False
        #print('self.enable_subdir_output', self.enable_subdir_output, flush=True)
        #print('self.rank', self.rank, flush=True)
        iter_n = 0
        pref_l, pref_u = self.preference_range
        prev_l, prev_u = self.preference_range
        closest_result = None
        
        while iter_n < self.max_sampled_preferences:
            #print('iter_n', iter_n, flush=True)
            '''
            if self.rank == 0:
                print('Beginning new iteration with iter_n being ' + str(iter_n), flush=True)
            '''
            self.preference = \
                    float(pref_l - pref_u) * float(self.rank + 1) / float(self.size + 1) + pref_u
            #print('self.preference', self.preference, flush=True)
            #print('self.size', self.size, flush=True)
            result = self._run()
            #print('type(result)', type(result), flush=True)
            num_of_clusters_gotten = result['num_of_clusters']
            #print('num_of_clusters_gotten', num_of_clusters_gotten, 'self.rank', self.rank, flush=True)
            #Package the number of clusters and preference value used by this rank
            # so that MPI communication latency is minimized.
            num_of_clusters_and_pref = [num_of_clusters_gotten, self.preference]
            
            num_of_clusters_and_pref_list = self.comm.gather(
                    num_of_clusters_and_pref, root=0)
            
            if self.rank == 0:
                if self.size == 1:
                    prev_l = pref_l
                    prev_u = pref_u
                new_pref_range_result = self.get_new_pref_range(
                        num_of_clusters_and_pref_list, prev_l, prev_u, iter_n)
            else:
                new_pref_range_result = None
            
            #Get the result of the clustering (success or failure and some 
            # other useful values)
            new_pref_range_result = self.comm.bcast(new_pref_range_result, root=0)
            '''
            if self.rank == 0:
                print('num_of_clusters_and_pref_list', num_of_clusters_and_pref_list, flush=True)
            '''
            #print('new_pref_range_result', new_pref_range_result, flush=True)
            success, place_holder_preference, n, pref_l, pref_u, prev_l, prev_u = new_pref_range_result
            self.preference_range = [pref_l, pref_u]
            
            result_num = result["num_of_clusters"]

            if success or (self.output_without_success and iter_n == self.max_sampled_preferences - 1):
                self.preference = place_holder_preference
                #print('success', flush=True)
                #If not the rank that was successful:
                if n != self.rank:
                    self.successful_rank = False
                    #print('rank ' + str(self.rank) + ' is returning', flush=True)
                    return
                else:
                    self.successful_rank = True
                
                #Print results! You're done!
                break
            else:
                self.successful_rank = False
            
            print("Iteration %i. Preference used: %f. "
                    "Clusters generated: %i."
                    % (iter_n, self.preference, result_num))
            self._print_result_summary(result)

            iter_n += 1
            
        if success:
            print("Affinity Propagation with fixed number of clusters "+
                    "succeeded!", flush=True)
        else:
            print("Failed to cluster to %i clusters with tolerance %i"
                    % (self.num_of_clusters, self.num_of_clusters_tolerance), flush=True)

        if success:
            self._print_results(result, self.verbose_output)
        elif self.output_without_success:
            closest_result = result
            self._print_results(closest_result, self.verbose_output)
        output_pool(result['exemplars'], self.exemplars_output_dir, self.exemplars_output_format)
            
        self.result = result
        
        #Unzip the 2-dimensional list which is the structure collection. The
        # second element in is the Structure object.
        pool = list(zip(*self.coll))[1]
        #In the future, parallelize the writing of structures to output_dir
        output_pool(pool, self.output_dir, self.output_format)
            
    def _run(self):
        #print('type(self.affinity_matrix)', type(self.affinity_matrix), flush=True)
        return self._affinity_propagation()

    def _print_results(self, result, verbose=False):
        if not self.output_file is None:
            #print(self.output_file, "Outputted iteration info:", flush=True)
            
            if verbose:
                #self._print_result_verbose(result)
                self._print_result_summary(result)
            else:
                self._print_result_summary(result)

    def _print_result_verbose(self, result):
        if self.output_file is None:
            raise IOError('output_file is None, cannot write output')

        self._print_result_summary(result)
        d = result
        assigned_cluster = d["assigned_cluster"]
        assigned_exemplar_id = d["assigned_exemplar_id"]
        distances_to_exemplar = d["distances_to_exemplar"]

        st = "List of assigned cluster, exemplar id, " \
                "and distance to exemplar:\n"
        for i in range(len(self.coll)):
            st += "%s %i %s %f\n" % (
                    self.coll[i][1].struct_id, assigned_cluster[i],
                    assigned_exemplar_id[i], distances_to_exemplar[i])

        exemplar_ids = d["exemplar_ids"]
        st += "Set of selected exemplars:\n"
        st += "\n".join(exemplar_ids)
        print(st, flush=True)

    def _print_result_summary(self, result):
        if self.output_file is None:
            raise IOError('output_file is None, cannot write output')

        d = result
        st = "Preference: %f. Number of clusters: %i.\n" \
                % (d["preference"], d["num_of_clusters"])
        #st += "Silhouette score: %f\n" % d["silhouette_score"]
        #st += "Mean distance to exemplar: %f\n" % d["avg_distance"]
        #st += "STD of distance to exemplar: %f\n" % d["std_distance"]
        #st += "Max distance to exemplar: %f\n" % d["max_distance"]
        #print(self.output_file, st, flush=True)


    def create_distance_matrix_from_exemplars(self, struct_ids, struct_ids_for_ap1):
        '''
        Purpose: After doing AP, you may want to do AP again on
            the exemplars. Therefore, you need to extract out the
            rows and columns associated with the exemplars in the 
            distance matrix to have a new distance matrix for 
            that second AP run.
        '''
        struct_ids_indices = [i for i in range(len(struct_ids_for_ap1)) if struct_ids_for_ap1[i] in struct_ids]
        self.distance_matrix = np.memmap(self.dist_mat_input_file_1, dtype='float32', mode='r')
        self.distance_matrix.resize((int(np.sqrt(len(self.distance_matrix))), int(np.sqrt(len(self.distance_matrix)))))
        self.distance_matrix = self.distance_matrix[struct_ids_indices,:][:,struct_ids_indices]

        fp = np.memmap(self.dist_mat_input_file, dtype='float32', mode='w+', shape=self.distance_matrix.shape)
        fp[:] = self.distance_matrix[:]

    def clean_up_tmp_affinity_matrix(self):
        file_utils.rm(self.my_affinity_matrix_fpath)

    def copy_affinity_matrix_on_disk(self):
        my_affinity_matrix_dir = os.path.basename(self.affinity_matrix_path)
        my_affinity_matrix_dir_fname = file_utils.fname_from_fpath(self.affinity_matrix_path) + '_' + str(self.rank) + '.dat'
        self.my_affinity_matrix_fpath = os.path.join(my_affinity_matrix_dir, my_affinity_matrix_dir_fname)
        file_utils.cp(self.affinity_matrix_path, my_affinity_matrix_dir, dest_fname=my_affinity_matrix_dir_fname, fail_if_cant_rm=True, overwrite=True)
        
    def _affinity_propagation(self):
        # my_affinity_matrix is altered during the course of the fit procedure so
        # retain a copy of the original on disk. We are making a copy for every rank
        # because then each rank can modify its own matrix. We are storing on disk
        # and using memmap to save on RAM usage.
        self.copy_affinity_matrix_on_disk()
        my_affinity_matrix = np.memmap(self.my_affinity_matrix_path, dtype='float32', mode='r+', shape=self.distance_matrix_shape)
        
        ap_iter_i = 0
        while ap_iter_i < self.max_iter:
            ap = AffinityPropagation(damping=self.damping, max_iter=self.max_iter, 
                                convergence_iter=self.convergence_iter,
                                copy=False, preference=self.preference,
                                affinity="precomputed",verbose=self.verbose_output)

            #print('inside _affinity_propagation', flush=True)
            #print('len(affinity_matrix)', len(affinity_matrix), flush=True)
            ##print('affinity_matrix', affinity_matrix, flush=True)
            if self.rank == 0:
                print('Fitting affinity propagation model...', flush=True)
            result = ap.fit(my_affinity_matrix)
        
            #print('result of fit', type(result), flush=True)
            assigned_cluster = result.labels_
            if assigned_cluster[0] != -1:
                # AP successful in converging on a clustering
                break
            if self.verbose_output:
                tmp_preference = math_utils.randrange_float(self.preference_range[0], self.preference_range[-1], 0.00001, num_decimal_places=5)
                print('preference of ' + str(self.preference) + ' failed to converge. Trying a random preference: '+ str(tmp_preference), flush=True)
                self.preference = tmp_preference
            ap_iter_i += 1
        if ap_iter_i == self.max_iter:
            raise Exception('AP was not able to find a preference that yielded a converged set of clsuters in',
                            self.max_iter, 'iterations. Check your affinity_type setting and your distance matrix values')

        num_of_clusters = len(result.cluster_centers_indices_)
        #print('num_of_clusters gotten on this fit', num_of_clusters, 'preference', self.preference, flush=True)
        
        if self.cluster_on_energy and not (iter_n != self.max_sampled_preferences - 1 and \
                (num_of_clusters > self.num_of_clusters + self.num_of_clusters_tolerance or \
                num_of_clusters < self.num_of_clusters - self.num_of_clusters_tolerance)):

            print('Using energy values to determine exemplars', flush=True)
            exemplar_indices = np.zeros(num_of_clusters)
            for i in range(num_of_clusters):
                lowest_energy = None
                for j,k in enumerate(assigned_cluster):
                    if i == k:
                        energy = float(self.coll[j][1].properties[self.energy_name])
                        #print('i', i, 'j', j, 'k', k, 'energy', energy, 'lowest_energy', lowest_energy, flush=True)
                        if lowest_energy is None or energy < lowest_energy:
                            lowest_energy = energy
                            #print('updated lowest_energy to', lowest_energy, flush=True)
                            exemplar_indices[i] = j

        else:
            exemplar_indices = result.cluster_centers_indices_
        #print('len(exemplar_indices)', len(exemplar_indices), flush=True)
        #print('len(self.coll)', len(self.coll), flush=True)
        #print('len(assigned_cluster)', len(assigned_cluster), flush=True)
        #print('assigned_cluster', assigned_cluster, flush=True)
        self.exemplar_indices = np.array(exemplar_indices, dtype='int')
        #print('self.exemplar_indices', self.exemplar_indices, flush=True)
        #print('self.coll', self.coll)
        exemplar_ids = [self.coll[x][1].struct_id for x in self.exemplar_indices]
        #print(coll, self.exemplar_indices, len(coll), len(self.exemplar_indices), flush=True)
        #print('len(coll[0])', len(coll[0]), flush=True)
        ##print('exemplar_ids', exemplar_ids, flush=True)
        
        #print('num_of_clusters',num_of_clusters, flush=True)
        #print('exemplar_ids', exemplar_ids, flush=True)

        assigned_exemplar_index = [self.exemplar_indices[x] for x in assigned_cluster]
        assigned_exemplar_id = [exemplar_ids[x] for x in assigned_cluster]
        #print('len(assigned_exemplar_id)', len(assigned_exemplar_id), flush=True)
        for x in range(len(self.coll)):
            self.coll[x][1].properties[self.property_key] = assigned_cluster[x]

        exemplars = [self.coll[x][1] for x in self.exemplar_indices]
        
        result_dict = {
                "num_of_clusters": num_of_clusters,
                "assigned_cluster": assigned_cluster,
                "assigned_exemplar_index": assigned_exemplar_index,
                "assigned_exemplar_id": assigned_exemplar_id,
                "exemplars": exemplars,
                "exemplar_indices": self.exemplar_indices,
                "exemplar_ids": exemplar_ids,
                "preference": self.preference}
        #output_pool(exemplars, exemplars_output_dir, exemplars_output_format)
        '''
        distances_to_exemplar = \
                [self.distance_matrix[x]
                        [self.exemplar_indices[result.labels_[x]]]
                        for x in range(len(self.coll))]
        avg_distance = np.mean(distances_to_exemplar)
        std_distance = np.std(distances_to_exemplar)
        max_distance = max(distances_to_exemplar)

        try:
            sil_score = silhouette_score(self.distance_matrix,
                    result.labels_, metric="precomputed")
        except:
            sil_score = 1.0
        
        result_dict = {
                "num_of_clusters": num_of_clusters,
                "assigned_cluster": assigned_cluster,
                "assigned_exemplar_index": assigned_exemplar_index,
                "assigned_exemplar_id": assigned_exemplar_id,
                "distances_to_exemplar": distances_to_exemplar,
                "exemplars": exemplars,
                "exemplar_indices": self.exemplar_indices,
                "exemplar_ids": exemplar_ids,
                "avg_distance": avg_distance,
                "std_distance": std_distance,
                "max_distance": max_distance,
                "silhouette_score": sil_score,
                "preference": self.preference}
        '''
        
        return result_dict


def affinity_propagation_fixed_clusters(inst, comm):
    '''
    AP that explores the setting of preference in order to generate
    desired number of clusters
    '''
    sname = "affinity_propagation_fixed_clusters"

    aph = APHandler(inst, sname, comm)
    print('got aph', flush=True)
    stoic = StoicDict(int)
    jsons_dir = aph.structure_dir
    #print('jsons_dir', jsons_dir, flush=True)
    aph.coll, struct_ids = get_struct_coll(jsons_dir, stoic)
    ##print('struct_ids', struct_ids, flush=True)
    #print('len(struct_ids)', len(struct_ids), flush=True)
    #print('type(aph.coll)', type(aph.coll), flush=True)
    ##print('aph.coll', aph.coll, flush=True)
    if aph.verbose_output:
        print('len(aph.coll)', len(aph.coll), flush=True)
    if aph.run_num == 2 and aph.rank == 0:
        _, struct_ids_for_ap1 = get_struct_coll(aph.structure_dir_for_ap1, stoic)
        aph.create_distance_matrix_from_exemplars(struct_ids,struct_ids_for_ap1)
        #if aph.verbose_output:
        #    print('len(aph.distance_matrix)0', len(aph.distance_matrix), flush=True)
    if aph.rank == 0:
        aph.make_affinity_matrix()
        time_utils.sleep(5)
    else:
        aph.distance_matrix_shape = None
    if aph.verbose_output:
        print('about to enter a barrier', flush=True)
    comm.barrier()
    if aph.verbose_output:
        print('got through the barrier', flush=True)
    aph.distance_matrix_shape = comm.bcast(aph.distance_matrix_shape, root=0)
    if aph.verbose_output:
        print('len(aph.distance_matrix_shape)', len(aph.distance_matrix_shape), flush=True)
    if not (aph.dist_mat_input_file).endswith('.dat'):
        raise Exception('Only supporting distance matrices saved as an np memmap with .dat extension')
    if aph.rank == 0:
        aph.get_affinity_and_distance_matrix()
    comm.barrier()
    #if aph.verbose_output:
    #    print('len(aph.distance_matrix)2', len(aph.distance_matrix), flush=True)
    rank_hostnames = comm.gather(gethostname(), root=0)
    if aph.rank == 0:
        # units kB
        available_mem = dict(psutil.virtual_memory()._asdict())['available'] / 1024.0
        allowed_number_of_ranks_per_node = aph.get_allowed_number_of_ranks_per_node(available_mem)
        assignments = aph.get_allowed_ranks_per_node(allowed_number_of_ranks_per_node, rank_hostnames)
    else:
        assignments = None

    assignment = comm.scatter(assignments, root=0)
    if assignment:    
        aph.run_fixed_num_of_clusters()
        #print('ran aph.run_fixed_num_of_clusters', flush=True)
        aph.clean_up_tmp_affinity_matrix()
    
    comm.barrier()
