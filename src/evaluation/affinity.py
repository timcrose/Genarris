'''
Created by Patrick Kilecdi on 12/31/2016

Conducts Affinity Propagation for a given collection
'''

import os
from sklearn.cluster import AffinityPropagation
from sklearn.metrics import silhouette_score
import numpy as np
from utilities.misc import output_pool
from utilities.write_log import print_time_log, write_log, \
    print_time_warning
from evaluation.evaluation_util import PoolOperation, \
        load_pool_operation_keywords

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

def affinity_propagation_fixed_clusters(inst):
    '''
    AP that explores the setting of preference in order to generate
    desired number of clusters
    '''
    sname = "affinity_propagation_fixed_clusters"
    num_of_clusters = inst.get_eval(sname,"num_of_clusters")
    preference_range = inst.get_eval(sname, "preference_range")
    num_of_clusters_tolerance = inst.get_with_default(
            sname, "num_of_clusters_tolerance", 0, eval=True)
    max_sampled_preferences = inst.get_with_default(
            sname, "max_sampled_preferences", 10)
    output_without_success = inst.get_boolean(
            sname, "output_without_success")

    executor = get_affinity_propagation_executor(inst, sname)
    executor.run_fixed_num_of_clusters(
            num_of_clusters, preference_range,
            num_of_clusters_tolerance=num_of_clusters_tolerance,
            max_sampled_preferences=max_sampled_preferences,
            output_without_success=output_without_success)
    

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
            section, "exemplars_output_format", "json")
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
        print_time_log("Running Affinity Propagation with preference: "
                + str(preference))
        # No need of subdir if only one run
        self._enable_subdir_output = False
        result = self._run(preference, enable_output=True)
        self._print_result_verbose(result)
        print_time_log("Affinity Propagation calculation complete!")

    def run_batch(self, preference_list, verbose_output=False):
        print_time_log("Running Affinity Propagation with preference list: "
                + str(preference_list))
        enable_output = False
        if not self._output_dir is None:
            self._enable_subdir_output = True
            enable_output = True

        for preference in preference_list:
            self._run_async(preference, enable_output)

        results = self.wait_for_all()
        self._print_results(results, verbose_output)
        print_time_log("Affinity Propagation calculation complete!")

    def run_fixed_num_of_clusters(self, num_of_clusters,
            preference_range, num_of_clusters_tolerance=0,
            max_sampled_preferences=10, output_without_success=False):
        if num_of_clusters > len(self._structure_list):
            raise ValueError("Cannot cluster pool into more clusters "
                    "than number of structures in it")
        if max_sampled_preferences < 1:
            raise ValueError("max_sampled_preferences must be >= 1")

        print_time_log("Running Affinity Propagation with fixed number "
                "of clusters " + str(num_of_clusters))
        self._enable_subdir_output = False

        iter_n = 0
        pref_l, pref_u = preference_range
        closest_result = None
        while iter_n < max_sampled_preferences:
            pref = (pref_l + pref_u) / 2.0
            result = self._run(pref, enable_output=False)
            closest_result = self._closer(
                    closest_result, result,
                    "num_of_clusters", num_of_clusters)

            result_num = result[1]["num_of_clusters"]
            print_time_log("Iteration %i. Preference used: %f. "
                    "Clusters generated: %i."
                    % (iter_n, pref, result_num))
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
            print_time_log("Affinity Propagation with fixed number of clusters "
                    "succeeded!")
            success = True
        else:
            print_time_log("Failed to cluster to %i clusters with tolerance %i"
                    % (num_of_clusters, num_of_clusters_tolerance))
            if output_without_success:
                print_time_warning("Outputing clustered collection although "
                        "fixed clustering failed.")

        if success or output_with_success:
            self._output_closest_result(closest_result)

    def run_with_fixed_silhouette(self, target_silhouette,
            preference_list, silhouette_tolerance=0.01,
            output_without_success=False):
        if len(preference_list) == 0:
            raise ValueError("Cannot run fixed silhouette with empty "
                    "preference_list")

        print_time_log("Running Affinity Propagation to reach target "
                "silhouette: " + str(target_silhouette))
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

        if abs(closest_result[1]["silhouette_score"] - target_silhouette) \
                <= silhouette_tolerance:
            print_time_log("Affinity Propagation with fixed silhouette score "
                    "succeeded!")
            success = True
        else:
            print_time_log("Failed to cluster to silhouette score %f "
                    "with tolerance %f"
                    % (target_silhouette, silhouette_tolerance))
            if output_without_success:
                print_time_warning("Outputing clustered collection although "
                        "fixed clustering failed.")

        if success or output_with_success:
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
            print_time_log("Exemplar output directory set as: " + output_dir)
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
        write_log(self._output_file, st, time_stamp=False)

    def _print_result_summary(self, result):
        if self._output_file is None: return

        d = result[1]
        st = "Preference: %f. Number of clusters: %i.\n" \
                % (d["preference"], d["num_of_clusters"])
        st += "Silhouette score: %f\n" % d["silhouette_score"]
        st += "Mean distance to exemplar: %f\n" % d["avg_distance"]
        st += "STD of distance to exemplar: %f\n" % d["std_distance"]
        st += "Max distance to exemplar: %f\n" % d["max_distance"]
        write_log(self._output_file, st, time_stamp=False)

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
            write_log(self._output_file, "Outputted iteration info:\n")
        self._print_result_verbose(closest_result)
        output_pool(closest_result[0],
                self._output_dir, self._output_format)
        output_pool(closest_result[1]["exemplars"],
                self._exemplars_output_dir,
                self._exemplars_output_format)


def _affinity_propagation(coll, distance_matrix, affinity_matrix,
        damping=0.5, convergence_iter=15, max_iter=20, preference=None,
        property_key="AP_cluster",
        exemplars_output_dir=None, exemplars_output_format="json"):
    ap = AffinityPropagation(damping=damping, max_iter=max_iter, 
                             convergence_iter=convergence_iter,
                             copy=True, preference=preference,
                             affinity="precomputed",verbose=False)
    
    result = ap.fit(affinity_matrix)

    exemplar_indices = result.cluster_centers_indices_
    exemplar_ids = [coll[x].struct_id for x in exemplar_indices]
    num_of_clusters = len(exemplar_indices)

    assigned_cluster = result.labels_
    assigned_exemplar_index = [exemplar_indices[x] for x in assigned_cluster]
    assigned_exemplar_id = [exemplar_ids[x] for x in assigned_cluster]

    for x in range(len(coll)):
        coll[x].properties[property_key] = assigned_cluster[x]

    exemplars = [coll[x] for x in result.cluster_centers_indices_]
    output_pool(exemplars, exemplars_output_dir, exemplars_output_format)

    distances_to_exemplar = \
            [distance_matrix[x]
                    [result.cluster_centers_indices_[result.labels_[x]]]
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
            "exemplar_indices": exemplar_indices,
            "exemplar_ids": exemplar_ids,
            "avg_distance": avg_distance,
            "std_distance": std_distance,
            "max_distance": max_distance,
            "silhouette_score": sil_score,
            "preference": preference}
    
    return coll, result_dict

