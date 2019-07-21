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

import os, random
from sklearn.cluster import AffinityPropagation
from sklearn.metrics import silhouette_score
import numpy as np
from Genarris.utilities.misc import output_pool
from Genarris.utilities import file_utils, time_utils
from Genarris.evaluation.evaluation_util import PoolOperation, \
        load_pool_operation_keywords
from Genarris.core.structure import StoicDict, get_struct_coll
from bisect import bisect



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
        self.num_of_clusters = inst.get_eval(sname,"num_of_clusters")
        self.preference_range = inst.get_eval(sname, "preference_range")
        self.num_of_clusters_tolerance = inst.get_with_default(
                sname, "num_of_clusters_tolerance", 0, eval=True)
        self.max_sampled_preferences = inst.get_with_default(
                sname, "max_sampled_preferences", 10, eval=True)
        self.output_without_success = inst.get_boolean(
                sname, "output_without_success")
        self.run_num = 1
        # Parameters for AP calculation
        self.affinity_type = inst.get_with_default(
                sname, "affinity_type", ["exponential",1], eval=True)
        self.damping = inst.get_with_default(
                sname, "damping", 0.5, eval=True)
        self.convergence_iter = inst.get_with_default(
                sname, "convergence_iter", 15, eval=True)
        self.max_iter = inst.get_with_default(
                sname, "max_iter", 200, eval=True)
        self.preference = inst.get_with_default(
                sname, 'preference', None, eval=True)
        self.property_key = inst.get_with_default(
                sname, "property_key", "AP_cluster")
        self.output_file = inst.get_with_default(
                sname, "output_file", "./AP_cluster.info")
        if os.path.exists(self.output_file):
            self.run_num = 2
            self.output_file = inst.get_with_default(
                sname, "output_file_2", self.output_file)
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
        
        self.structure_dir_depth = inst.get_with_default(sname, 'structure_dir_depth' , 0, eval=True)
        self.structure_suffix = inst.get_with_default(sname, 'structure_suffix', '.json')
        self.output_dir = inst.get(sname, 'output_dir')
        if os.path.exists(self.output_dir):
            self.run_num = 2
            self.output_dir = inst.get(sname, 'output_dir_2')
            self.num_of_clusters = inst.get_eval(sname,"num_of_clusters_2")

        self.dist_mat_input_file = inst.get_inferred(sname, [sname, 'rcd_difference_folder_inner', 'rcd_calculation'], 
                                                    ['dist_mat_input_file', 'diff_matrix_output', 'diff_matrix_output'])

        if self.run_num == 1:
            self.structure_dir = inst.get_inferred(sname, [sname, 'rcd_difference_folder_inner', 'rcd_calculation'], ['structure_dir'] + (2 * ['output_dir']), type_='dir')
        elif self.run_num == 2:
            self.structure_dir = inst.get(sname, 'exemplars_output_dir')
            if not os.path.isdir(self.structure_dir):
                raise Exception('self.structure_dir', self.structure_dir, 'path DNE')
            ext_pos = self.dist_mat_input_file.find('.')
            self.dist_mat_input_file = self.dist_mat_input_file[:ext_pos] + '1' + self.dist_mat_input_file[ext_pos:]
        self.output_format = inst.get_with_default(sname, 'output_format', 'both')

        if self.run_num == 1:
            self.cluster_on_energy = inst.get_boolean(sname, 'cluster_on_energy')
            sname_list = [sname, 'fhi_aims_energy_evaluation', 'harris_approximation_batch']
            self.energy_name = inst.get_inferred(sname, sname_list, ['energy_name'] * len(sname_list))
        elif self.run_num == 2:
            self.cluster_on_energy = inst.get_inferred(sname, [sname] * 2, ['cluster_on_energy_2', 'cluster_on_energy'], type_=bool)
            sname_list = [sname, sname, 'fhi_aims_energy_evaluation', 'harris_approximation_batch']
            self.energy_name = inst.get_inferred(sname, sname_list, ['energy_name_2'] + ['energy_name'] * (len(sname_list) - 1))
        #print('self.cluster_on_energy', self.cluster_on_energy, flush=True)
        #print('self.energy_name', self.energy_name, flush=True)
        
        #Implement the affinity type desired
        self._initialize_affinity_matrix()

    def _initialize_affinity_matrix(self):
        f = open(self.dist_mat_input_file, "r")
        lines = f.read().split("\n")
        while lines[-1] == "":
            lines.pop()
        dist_mat = [[float(x) for x in y.split()] for y in lines]
        
        if len(dist_mat) != len(dist_mat[0]):
            raise ValueError("Distance matrix is not a square matrix: "
                    + self.dist_mat_input_file)
        
        m = len(dist_mat)
        self.distance_matrix = dist_mat

        self.affinity_matrix = [[0] * m for x in range(m)]
        
        for i in range(m):
            for j in range(m):
                if self.affinity_type[0]=="exponential":
                    self.affinity_matrix[i][j] = \
                            -np.exp(dist_mat[i][j] *
                                    self.affinity_type[1])
                elif self.affinity_type[0] == "power":
                    self.affinity_matrix[i][j] = \
                            -dist_mat[i][j] ** self.affinity_type[1]
        
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
        if self.size != 1:
            prev_l, prev_u = min(pref_list), max(pref_list)
            resample_prob = 0.2
        else:
            resample_prob = 0.02
        if ((sorted(num_clusters_list) != num_clusters_list or \
                num_clusters_list == [num_clusters_list[0]] * len(num_clusters_list)) and\
                    self.size != 1) or\
                 random.random() < resample_prob or\
                     pref_u - pref_l < 0.00001:
            pref_l -= 50.0
            pref_u += 50.0
            prev_l -= 50.0
            prev_u += 50.0
        
        #In this case, the 0 values are just there for formatting
        return success, 0, 0, pref_l, pref_u, prev_l, prev_u 
                            
    def run_fixed_num_of_clusters(self):
        '''
        Run the AP algorithm in parallel until either an acceptable number of 
        of clusters has been generated or the number of iterations asked for 
        has expired.
        '''
        #print('inside run_fixed_num_of_clusters', flush=True)
        
        if self.num_of_clusters > len(self.struct_coll):
            raise ValueError("Cannot cluster pool into more clusters " +
                    "than number of structures in it. self.num_of_clusters = " + 
                    str(self.num_of_clusters) + ", len(self.struct_coll) = " +
                    str(len(self.struct_coll)))
        #print('self.num_of_clusters', self.num_of_clusters, flush=True)
        #print('len(self.struct_coll)', len(self.struct_coll), flush=True)
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
            if self.rank == 0:
                print('Beginning new iteration with iter_n being ' + str(iter_n), flush=True)
            print('self.preference', self.preference, flush=True)
            self.preference = \
                    float(pref_l - pref_u) * float(self.rank + 1) / float(self.size + 1) + pref_u
            #print('self.size', self.size, flush=True)
            self.struct_coll_cld, result = self._run()
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
            
            if self.rank == 0:
                print('num_of_clusters_and_pref_list', num_of_clusters_and_pref_list, flush=True)
            print('new_pref_range_result', new_pref_range_result, flush=True)
            success, self.preference, n, pref_l, pref_u, prev_l, prev_u = new_pref_range_result
            
            result_num = result["num_of_clusters"]

            if success:
                #print('success', flush=True)
                #If not the rank that was successful:
                if n != self.rank:
                    #print('rank ' + str(self.rank) + ' is returning', flush=True)
                    return
                
                #Print results! You're done!
                break
            
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
        pool = list(zip(*self.struct_coll_cld))[1]
        #In the future, parallelize the writing of structures to output_dir
        output_pool(pool, self.output_dir, self.output_format)
            
    def _run(self):
        #print('type(self.affinity_matrix)', type(self.affinity_matrix), flush=True)
        return self._affinity_propagation(self.struct_coll, 
            self.distance_matrix, self.affinity_matrix,
            self.damping, self.convergence_iter, self.max_iter, self.preference,
            self.property_key,
            self.exemplars_output_dir, self.exemplars_output_format)

    def _print_results(self, result, verbose=False):
        if not self.output_file is None:
            print(self.output_file, "Outputted iteration info:", flush=True)
            
            if verbose:
                self._print_result_verbose(result)
            else:
                self._print_result_summary(result)

    def _print_result_verbose(self, result):
        if self.output_file is None:
            raise IOError('output_file is None, cannot write output')

        self._print_result_summary(result)
        coll = self.struct_coll_cld; d = result
        assigned_cluster = d["assigned_cluster"]
        assigned_exemplar_id = d["assigned_exemplar_id"]
        distances_to_exemplar = d["distances_to_exemplar"]

        st = "List of assigned cluster, exemplar id, " \
                "and distance to exemplar:\n"
        for i in range(len(coll)):
            st += "%s %i %s %f\n" % (
                    coll[i][1].struct_id, assigned_cluster[i],
                    assigned_exemplar_id[i], distances_to_exemplar[i])

        exemplar_ids = d["exemplar_ids"]
        st += "Set of selected exemplars:\n"
        st += "\n".join(exemplar_ids)
        print(self.output_file, st, flush=True)

    def _print_result_summary(self, result):
        if self.output_file is None:
            raise IOError('output_file is None, cannot write output')

        d = result
        st = "Preference: %f. Number of clusters: %i.\n" \
                % (d["preference"], d["num_of_clusters"])
        st += "Silhouette score: %f\n" % d["silhouette_score"]
        st += "Mean distance to exemplar: %f\n" % d["avg_distance"]
        st += "STD of distance to exemplar: %f\n" % d["std_distance"]
        st += "Max distance to exemplar: %f\n" % d["max_distance"]
        print(self.output_file, st, flush=True)


    def create_distance_matrix_from_exemplars(self):
        '''
        Purpose: After doing AP, you may want to do AP again on
            the exemplars. Therefore, you need to extract out the
            rows and columns associated with the exemplars in the 
            distance matrix to have a new distance matrix for 
            that second AP run.
        '''
        f = open(self.dist_mat_input_file, "r")
        lines = f.read().split("\n")
        while lines[-1] == "":
            lines.pop()
        dist_mat = [np.array([float(x) for x in y.split()]) for y in lines]
        f.close()
        
        if len(dist_mat) != len(dist_mat[0]):
            raise ValueError("Distance matrix is not a square matrix: "
                    + self.dist_mat_input_file)
        
        dist_mat = np.array(dist_mat)
        dist_mat = dist_mat[self.exemplar_indices,:][:,self.exemplar_indices]
        dist_mat = list(map(list, dist_mat))

        ext_pos = self.dist_mat_input_file.find('.')
        self.dist_mat_input_file_2 = self.dist_mat_input_file[:ext_pos] + '1' + self.dist_mat_input_file[ext_pos:]
        file_utils.write_rows_to_csv(self.dist_mat_input_file_2, dist_mat, mode='w', delimiter=' ')


    def _affinity_propagation(self, coll, distance_matrix, affinity_matrix,
            damping=0.5, convergence_iter=15, max_iter=20, preference=None,
            property_key="AP_cluster",
            exemplars_output_dir=None, exemplars_output_format="json"):
        ap = AffinityPropagation(damping=damping, max_iter=max_iter, 
                                convergence_iter=convergence_iter,
                                copy=True, preference=preference,
                                affinity="precomputed",verbose=False)
        #print('inside _affinity_propagation', flush=True)
        #print('len(affinity_matrix)', len(affinity_matrix), flush=True)
        ##print('affinity_matrix', affinity_matrix, flush=True)
        print('Fitting affinity propagation model...', flush=True)
        
        result = ap.fit(affinity_matrix)
        print('AP fit complete.', flush=True)
        #print('result of fit', type(result), flush=True)
        num_of_clusters = len(result.cluster_centers_indices_)
        
        assigned_cluster = result.labels_
        
        if self.cluster_on_energy and not (iter_n != self.max_sampled_preferences - 1 and \
                (num_of_clusters > self.num_of_clusters + self.num_of_clusters_tolerance or \
                num_of_clusters < self.num_of_clusters - self.num_of_clusters_tolerance)):

            print('Using energy values to determine exemplars', flush=True)
            exemplar_indices = np.zeros(num_of_clusters)
            for i in range(num_of_clusters):
                lowest_energy = None
                for j,k in enumerate(assigned_cluster):
                    if i == k:
                        energy = float(coll[j][1].properties[self.energy_name])
                        #print('i', i, 'j', j, 'k', k, 'energy', energy, 'lowest_energy', lowest_energy, flush=True)
                        if lowest_energy is None or energy < lowest_energy:
                            lowest_energy = energy
                            #print('updated lowest_energy to', lowest_energy, flush=True)
                            exemplar_indices[i] = j

        else:
            exemplar_indices = result.cluster_centers_indices_
        #print('len(exemplar_indices)', len(exemplar_indices), flush=True)
        #print('len(coll)', len(coll), flush=True)
        #print('len(assigned_cluster)', len(assigned_cluster), flush=True)
        self.exemplar_indices = np.array(exemplar_indices, dtype='int')
        #print('self.exemplar_indices', self.exemplar_indices, flush=True)
        ##print('self.exemplar_indices', self.exemplar_indices, flush=True)
        exemplar_ids = [coll[x][1].struct_id for x in self.exemplar_indices]
        #print(coll, self.exemplar_indices, len(coll), len(self.exemplar_indices), flush=True)
        #print('len(coll[0])', len(coll[0]), flush=True)
        ##print('exemplar_ids', exemplar_ids, flush=True)
        
        #print('num_of_clusters',num_of_clusters, flush=True)
        #print('exemplar_ids', exemplar_ids, flush=True)

        assigned_exemplar_index = [self.exemplar_indices[x] for x in assigned_cluster]
        assigned_exemplar_id = [exemplar_ids[x] for x in assigned_cluster]

        for x in range(len(coll)):
            coll[x][1].properties[property_key] = assigned_cluster[x]

        exemplars = [coll[x][1] for x in self.exemplar_indices]
        
        #output_pool(exemplars, exemplars_output_dir, exemplars_output_format)

        distances_to_exemplar = \
                [distance_matrix[x]
                        [self.exemplar_indices[result.labels_[x]]]
                        for x in range(len(coll))]
        avg_distance = np.mean(distances_to_exemplar)
        std_distance = np.std(distances_to_exemplar)
        max_distance = max(distances_to_exemplar)

        try:
            sil_score = silhouette_score(np.array(distance_matrix),
                    result.labels_, metric="precomputed")
        except:
            sil_score = 1.0

        if preference is None:
            preference = np.median([x for y in affinity_matrix for x in y])

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
                "preference": preference}
        
        return coll, result_dict


def affinity_propagation_distance_matrix(inst):
    '''
    Main AP distance matrix module
    '''
    sname = "affinity_propagation_distance_matrix"
    preference = inst.get_or_none(sname, "preference", eval=True)
    executor = get_affinity_propagation_executor(inst, sname)
    executor.run_single(preference)

def affinity_propagation_analyze_preference(inst):
    '''
    Run AP with a list of preferences
    '''
    sname = "affinity_propagation_analyze_preference"
    preference_list = inst.get_eval(sname,"preference_list")
    verbose_output = inst.get_boolean(sname, "verbose_output")
    executor = get_affinity_propagation_executor(inst, sname)
    executor.run_batch(preference_list, verbose_output=verbose_output)

def affinity_propagation_fixed_clusters(inst, comm):
    '''
    AP that explores the setting of preference in order to generate
    desired number of clusters
    '''
    sname = "affinity_propagation_fixed_clusters"
    #print('inside ap', flush=True)
    aph = APHandler(inst, sname, comm)
    #print('got aph', flush=True)
    stoic = StoicDict(int)
    jsons_dir = aph.structure_dir
    #print('jsons_dir', jsons_dir, flush=True)
    aph.struct_coll, struct_ids = get_struct_coll(jsons_dir, stoic)
    ##print('struct_ids', struct_ids, flush=True)
    #print('len(struct_ids)', len(struct_ids), flush=True)
    #print('type(aph.struct_coll)', type(aph.struct_coll), flush=True)
    ##print('aph.struct_coll', aph.struct_coll, flush=True)
    #print('len(aph.struct_coll)', len(aph.struct_coll), flush=True)
    
    aph.run_fixed_num_of_clusters()
    #print('ran aph.run_fixed_num_of_clusters', flush=True)
    if aph.run_num == 1:
        #print('creating distance matrix', flush=True)
        aph.create_distance_matrix_from_exemplars()
        #print('created dist mat', flush=True)
    
    #aph.struct_coll_cld, aph.result
    
    
    
    '''
    executor = get_affinity_propagation_executor(inst, sname)
    executor.run_fixed_num_of_clusters(
            num_of_clusters, preference_range,
            num_of_clusters_tolerance=num_of_clusters_tolerance,
            max_sampled_preferences=max_sampled_preferences,
            output_without_success=output_without_success)
    '''
    

def affinity_propagation_fixed_silhouette(inst):
    '''
    Run AP on the given preference list and return the one within
    the tolerance of the target_silhouette
    '''
    sname = "affinity_propagation_fixed_silhouette"
    target_silhouette = inst.get_eval(sname, "target_silhouette")
    preference_list = inst.get_eval(sname, "preference_list")
    silhouette_tolerance = inst.get_with_default(
            sname, "silhouette_tolerance", 0.01, eval=True)
    output_without_success = inst.get_boolean(
            sname, "output_without_success")

    executor = get_affinity_propagation_executor(inst, sname)
    executor.run_with_fixed_silhouette(
            target_silhouette, preference_list,
            silhouette_tolerance=silhouette_tolerance,
            output_without_success=output_without_success)

def get_affinity_propagation_executor(inst, section):
    dist_mat_input_file = inst.get(section, "dist_mat_input_file")
    affinity_type = inst.get_with_default(
            section, "affinity_type", ["exponential",1], eval=True)
    damping = inst.get_with_default(
            section, "damping", 0.5, eval=True)
    convergence_iter = inst.get_with_default(
            section, "convergence_iter", 15, eval=True)
    max_iter = inst.get_with_default(
            section, "max_iter", 200, eval=True)
    property_key = inst.get_with_default(
            section, "property_key", "AP_cluster")
    output_file = inst.get_with_default(
            section, "output_file", "./AP_cluster.info")
    exemplars_output_dir = inst.get_or_none(
            section, "exemplars_output_dir")
    exemplars_output_format = inst.get_with_default(
            section, "exemplars_output_format", "both")
    pool_operation_keywords = load_pool_operation_keywords(
            inst, section)

    return AffinityPropagationExecutor(
            dist_mat_input_file, affinity_type=affinity_type,
            damping=damping, convergence_iter=convergence_iter,
            max_iter=max_iter, property_key=property_key,
            output_file=output_file,
            exemplars_output_dir=exemplars_output_dir,
            exemplars_output_format=exemplars_output_format,
            **pool_operation_keywords)

class AffinityPropagationExecutor(PoolOperation):
    def __init__(self, dist_mat_input_file,
            affinity_type=["exponential",1], damping=0.5,
            convergence_iter=15, max_iter=200,
            property_key="AP_cluster",
            output_file="./AP_cluster.info",
            exemplars_output_dir=None, exemplars_output_format="json",
            # Below are variables inherited from the PoolOperation class
            structure_list=None, structure_dir=None,
            structure_dir_depth=0, structure_suffix="",
            output_dir=None, output_format="json",
            processes_limit=1):
        self._dist_mat_input_file = dist_mat_input_file
        self._affinity_type = affinity_type
        self._damping = damping
        self._convergence_iter = convergence_iter
        self._max_iter = max_iter
        self._property_key = property_key
        self._output_file = output_file
        self._exemplars_output_dir = exemplars_output_dir
        self._exemplars_output_format = exemplars_output_format
        PoolOperation.__init__(self,
                structure_list=structure_list,
                structure_dir=structure_dir,
                structure_dir_depth=structure_dir_depth,
                structure_suffix=structure_suffix,
                output_dir=output_dir,
                output_format=output_format,
                enable_subdir_output=False,
                processes_limit=processes_limit)

        self._initialize_affinity_matrix()
        if len(self._affinity_matrix) != \
                len(self._structure_list):
            raise ValueError("Distance matrix row number does not equal "
                    "to number of input structures")

        self._structure_list.sort(key=lambda x: x.struct_id)

    def _initialize_affinity_matrix(self):
        f = open(self._dist_mat_input_file, "r")
        lines = f.read().split("\n")
        while lines[-1] == "":
            lines.pop()
        dist_mat = [[float(x) for x in y.split()] for y in lines]

        if len(dist_mat) != len(dist_mat[0]):
            raise ValueError("Distance matrix is not a square matrix: "
                    + self._dist_mat_input_file)

        m = len(dist_mat)
        self._distance_matrix = dist_mat
        self._affinity_matrix = [[0]*m for x in range(m)]

        for i in range(m):
            for j in range(m):
                if self._affinity_type[0]=="exponential":
                    self._affinity_matrix[i][j] = \
                            -np.exp(dist_mat[i][j] *
                                    self._affinity_type[1])
                elif self._affinity_type[0] == "power":
                    self._affinity_matrix[i][j] = \
                            -dist_mat[i][j] ** self._affinity_type[1]
    
    def run_single(self, preference):
        print("Running Affinity Propagation with preference: "
                + str(preference), flush=True)
        # No need of subdir if only one run
        self._enable_subdir_output = False
        result = self._run(preference, enable_output=True)
        self._print_result_verbose(result)
        print("Affinity Propagation calculation complete!", flush=True)

    def run_batch(self, preference_list, verbose_output=False):
        print("Running Affinity Propagation with preference list: "
                + str(preference_list), flush=True)
        enable_output = False
        if not self._output_dir is None:
            self._enable_subdir_output = True
            enable_output = True

        for preference in preference_list:
            self._run_async(preference, enable_output)

        results = self.wait_for_all()
        self._print_results(results, verbose_output)
        print("Affinity Propagation batch calculation complete!", flush=True)

    def run_fixed_num_of_clusters(self, num_of_clusters,
            preference_range, num_of_clusters_tolerance=0,
            max_sampled_preferences=10, output_without_success=False):
        if num_of_clusters > len(self._structure_list):
            raise ValueError("Cannot cluster pool into more clusters "
                    "than number of structures in it")
        if max_sampled_preferences < 1:
            raise ValueError("max_sampled_preferences must be >= 1")

        print("Running Affinity Propagation with fixed number "
                "of clusters " + str(num_of_clusters), flush=True)
        self._enable_subdir_output = False

        iter_n = 0
        pref_l, pref_u = preference_range
        closest_result = None
        while iter_n < max_sampled_preferences:
            pref = (pref_l + pref_u) / 2.0
            struct_coll_cld, result = self._run(pref, enable_output=False)
            if result == False:
                continue
            closest_result = self._closer(
                    closest_result, result,
                    "num_of_clusters", num_of_clusters)

            result_num = result[1]["num_of_clusters"]
            print("Iteration %i. Preference used: %f. "
                    "Clusters generated: %i."
                    % (iter_n, pref, result_num), flush=True)
            self._print_result_summary(result)

            if result_num == num_of_clusters:
                break
            if result_num > num_of_clusters:
                pref_u = pref
            else:
                pref_l = pref

            iter_n += 1

        success = False
        if abs(closest_result[1]["num_of_clusters"] - num_of_clusters) \
                <= num_of_clusters_tolerance:
            print("Affinity Propagation with fixed number of clusters "
                    "succeeded!", flush=True)
            success = True
        else:
            print("Failed to cluster to %i clusters with tolerance %i"
                    % (num_of_clusters, num_of_clusters_tolerance), flush=True)
            if output_without_success:
                print("Outputing clustered collection although "
                        "fixed clustering failed.", flush=True)

        if success or output_without_success:
            self._output_closest_result(closest_result)

    def run_with_fixed_silhouette(self, target_silhouette,
            preference_list, silhouette_tolerance=0.01,
            output_without_success=False):
        if len(preference_list) == 0:
            raise ValueError("Cannot run fixed silhouette with empty "
                    "preference_list")

        print("Running Affinity Propagation to reach target "
                "silhouette: " + str(target_silhouette), flush=True)
        self._enable_subdir_output = False
        
        for preference in preference_list:
            self._run_async(preference, enable_output=False)

        results = self.wait_for_all()
        results.sort(key=lambda x: x[1]["preference"])
        closest_result = None
        for result in results:
            closest_result = self._closer(
                    closest_result, result, "silhouette_score",
                    target_silhouette)
            self._print_result_summary(result)

        success = False
        if abs(closest_result[1]["silhouette_score"] - target_silhouette) \
                <= silhouette_tolerance:
            print("Affinity Propagation with fixed silhouette score "
                    "succeeded!", flush=True)
            success = True
        else:
            print("Failed to cluster to silhouette score %f "
                    "with tolerance %f"
                    % (target_silhouette, silhouette_tolerance), flush=True)
            if output_without_success:
                print("Outputing clustered collection although "
                        "fixed clustering failed.", flush=True)

        if success or output_without_success:
            self._output_closest_result(closest_result)

    def _run(self, preference, enable_output=True):
        name = "AffinityPropagationWithPreference" + str(preference)

        exemplars_output_dir = self._get_exemplars_output_dir(
                name, enable_output)

        return self.run_operation(
                _affinity_propagation,
                name=name,
                args=self._get_args(),
                kwargs=self._get_kwargs(preference, exemplars_output_dir),
                enable_output=enable_output)

    def _run_async(self, preference, enable_output=True):
        name = "AffinityPropagationWithPreference" + str(preference)

        exemplars_output_dir = self._get_exemplars_output_dir(
                name, enable_output)

        return self.run_operation_async(
                _affinity_propagation,
                name=name,
                args=self._get_args(),
                kwargs=self._get_kwargs(preference, exemplars_output_dir),
                enable_output=enable_output)
  
    def _get_args(self):
        return (self._distance_matrix, self._affinity_matrix)

    def _get_kwargs(self, preference, exemplars_output_dir):
        return {"damping": self._damping,
                "convergence_iter": self._convergence_iter,
                "max_iter": self._max_iter,
                "preference": preference,
                "property_key": self._property_key,
                "exemplars_output_dir": exemplars_output_dir,
                "exemplars_output_format": self._exemplars_output_format}

    def _get_exemplars_output_dir(self, name, enable_output=True):
        if enable_output and not self._exemplars_output_dir is None:
            if self._enable_subdir_output:
                output_dir = os.path.join(
                        self._exemplars_output_dir, name)
            else:
                output_dir = self._exemplars_output_dir
            print("Exemplar output directory set as: " + output_dir, flush=True)
        else:
            output_dir = None
        return output_dir

    def _print_results(self, results, verbose=False):
        results.sort(key=lambda result: result[1]["preference"])
        for result in results:
            if verbose:
                self._print_result_verbose(result)
            else:
                self._print_result_summary(result)

    def _print_result_verbose(self, result):
        if self._output_file is None: return

        self._print_result_summary(result)
        coll = result[0]; d = result[1]
        assigned_cluster = d["assigned_cluster"]
        assigned_exemplar_id = d["assigned_exemplar_id"]
        distances_to_exemplar = d["distances_to_exemplar"]

        st = "List of assigned cluster, exemplar id, " \
                "and distance to exemplar:\n"
        for i in range(len(coll)):
            st += "%s %i %s %f\n" % (
                    coll[i].struct_id, assigned_cluster[i],
                    assigned_exemplar_id[i], distances_to_exemplar[i])

        exemplar_ids = d["exemplar_ids"]
        st += "Set of selected exemplars:\n"
        st += "\n".join(exemplar_ids)
        print(self._output_file, st, flush=True)

    def _print_result_summary(self, result):
        if self._output_file is None: return

        d = result[1]
        st = "Preference: %f. Number of clusters: %i.\n" \
                % (d["preference"], d["num_of_clusters"])
        st += "Silhouette score: %f\n" % d["silhouette_score"]
        st += "Mean distance to exemplar: %f\n" % d["avg_distance"]
        st += "STD of distance to exemplar: %f\n" % d["std_distance"]
        st += "Max distance to exemplar: %f\n" % d["max_distance"]
        print(self._output_file, st, flush=True)

    def _closer(self, result_1, result_2, key, target):
        if result_1 is None:
            return result_2

        diff_1 = abs(result_1[1][key] - target)
        diff_2 = abs(result_2[1][key] - target)
        if diff_1 <= diff_2:
            return result_1
        else:
            return result_2

    def _output_closest_result(self, closest_result):
        if not self._output_file is None:
            print(self._output_file, "Outputted iteration info:", flush=True)
        self._print_result_verbose(closest_result)
        output_pool(closest_result[0],
                self._output_dir, self._output_format)

        output_pool(closest_result[1]["exemplars"],
                self._exemplars_output_dir,
                self._exemplars_output_format)




