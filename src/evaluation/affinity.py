'''
Created by Patrick Kilecdi on 12/31/2016

Conducts Affinity Propagation for a given collection
'''

import os
from sklearn.cluster import AffinityPropagation
from sklearn.metrics import silhouette_score
import numpy as np
from utilities import misc
from utilities.misc import output_pool
from utilities.write_log import print_time_log, write_log
from evaluation.evaluation_util import PoolOperation, \
        load_pool_operation_keywords

def affinity_propagation_distance_matrix(inst):
    '''
    Main AP distance matrix module
    '''
    sname = "affinity_propagation_distance_matrix"
    executor = get_affinity_propagation_executor(inst, sname)
    executor.run_single()

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
    AP that explores the setting of preference in order to generate desired number of clusters
    '''
    sname = "affinity_propagation_fixed_clusters"
    dist_mat_path = inst.get(sname,"dist_mat_input_file")
    n_clusters = inst.get_eval(sname,"number_of_clusters")
    pref_range = inst.get_eval(sname,"preference_range")
    pref_iter = inst.get_with_default(sname,"pref_test_max_iters",10,eval=True)
    affinity_type = inst.get_with_default(sname,"affinity_type",["exponential",1],eval=True)
    damping = inst.get_with_default(sname,"damping",0.5,eval=True)
    convergence_iter = inst.get_with_default(sname,"convergence_iter",15,eval=True)
    max_iter = inst.get_with_default(sname,"max_iter",200,eval=True)
#    preference = inst.get_with_default(sname,"preference",None,eval=True)
    stored_property_key = inst.get_with_default(sname,"stored_property_key","AP_cluster")
    cluster_output_file = inst.get_with_default(sname,"cluster_output_file","./AP_cluster.info")
    
    coll = misc.load_collection_with_inst(inst,sname)
    coll.sort(key=lambda struct: struct.struct_id)
    
    f = open(dist_mat_path,"r")
    lines = f.read().split("\n")
    while lines[-1]=="":
        lines.pop()
    dist_mat = [[float(x) for x in y.split()] for y in lines]

    pref_l, pref_u = pref_range

    iter_n = 0
    while pref_iter > iter_n: #Iterate to search for the good preference that gives the amount of clusters
        pref = (pref_l + pref_u) / 2.0
        preference = _generate_preference_array(coll,["constant",pref])
    
        coll, centers = AP_distance_matrix(coll, dist_mat, affinity_type=affinity_type,
                                           damping=damping, convergence_iter=convergence_iter,
                                           max_iter=max_iter, preference=preference,
                                           stored_property_key=stored_property_key)

        labels = np.array([struct.properties[stored_property_key] for struct in coll])
     
        f = open(cluster_output_file,"a")

        if len(centers) < len(coll):
             f.write("Iteration %i; preference used: %f; clusters generated %i; silhouette score: %f\n" 
                    % (iter_n, pref, len(centers), 
                       silhouette_score(np.array(dist_mat),labels,metric="precomputed")))
        else:
             f.write("Iteration %i; preference used: %f; clusters generated %i\n" 
                    % (iter_n, pref, len(centers)))

        f.close()
           
        iter_n += 1


        if len(centers) == n_clusters:
            break

        if len(centers) > n_clusters:
            pref_u = pref
        elif len(centers) < n_clusters:
            pref_l = pref   


    dists = [dist_mat[i][coll.index(centers[coll[i].properties[stored_property_key]])]
                for i in range(len(coll))]

    f = open(cluster_output_file,"a")
  
    f.write("\nFinal iteration clustering result:\n")
    f.write(_print_ap_cluster_information(centers, dists, coll, stored_property_key, dist_mat))
    f.close()

    if inst.has_option(sname,"output_folder"):
        misc.dump_collection_with_inst(inst,sname,coll)

    if inst.has_option(sname,"exemplar_output_folder"):
        op = inst.get(sname,"output_folder")
        inst.set(sname,"output_folder",inst.get(sname,"exemplar_output_folder"))
        misc.dump_collection_with_inst(inst,sname,centers)
        inst.set(sname,"output_folder",op) 

def affinity_propagation_fixed_silhouette(inst):
    '''
    AP that explores the setting of preference in order to generate desired number of clusters
    '''
    sname = "affinity_propagation_fixed_clusters"
    dist_mat_path = inst.get(sname,"dist_mat_input_file")
    desire_sil = inst.get_eval(sname,"desired_silhouette")
    pref_range = inst.get_eval(sname,"preference_range")
    pref_iter = inst.get_with_default(sname,"pref_test_max_iters",10,eval=True)
    affinity_type = inst.get_with_default(sname,"affinity_type",["exponential",1],eval=True)
    damping = inst.get_with_default(sname,"damping",0.5,eval=True)
    convergence_iter = inst.get_with_default(sname,"convergence_iter",15,eval=True)
    max_iter = inst.get_with_default(sname,"max_iter",200,eval=True)
#    preference = inst.get_with_default(sname,"preference",None,eval=True)
    stored_property_key = inst.get_with_default(sname,"stored_property_key","AP_cluster")
    cluster_output_file = inst.get_with_default(sname,"cluster_output_file","./AP_cluster.info")
    
    coll = misc.load_collection_with_inst(inst,sname)
    coll.sort(key=lambda struct: struct.struct_id)
    
    f = open(dist_mat_path,"r")
    lines = f.read().split("\n")
    while lines[-1]=="":
        lines.pop()
    dist_mat = [[float(x) for x in y.split()] for y in lines]

    pref_l, pref_u = pref_range

    iter_n = 0
    while pref_iter > iter_n: #Iterate to search for the good preference that gives the amount of clusters
        pref = (pref_l + pref_u) / 2.0
        preference = _generate_preference_array(coll,["constant",pref])
    
        coll, centers = AP_distance_matrix(coll, dist_mat, affinity_type=affinity_type,
                                           damping=damping, convergence_iter=convergence_iter,
                                           max_iter=max_iter, preference=preference,
                                           stored_property_key=stored_property_key)

        labels = np.array([struct.properties[stored_property_key] for struct in coll])
     
        sil_score = silhouette_score(np.array(dist_mat),labels,metric="precomputed")
        f = open(cluster_output_file,"a")

        if sil_score < len(coll):
             f.write("Iteration %i; preference used: %f; clusters generated %i; silhouette score: %f\n" 
                    % (iter_n, pref, len(centers), 
                       silhouette_score(np.array(dist_mat),labels,metric="precomputed")))
        else:
             f.write("Iteration %i; preference used: %f; clusters generated %i\n" 
                    % (iter_n, pref, len(centers)))

        f.close()
           
        iter_n += 1


        if len(centers) == n_clusters:
            break

        if len(centers) > n_clusters:
            pref_u = pref
        elif len(centers) < n_clusters:
            pref_l = pref   


    dists = [dist_mat[i][coll.index(centers[coll[i].properties[stored_property_key]])]
                for i in range(len(coll))]

    f = open(cluster_output_file,"a")
  
    f.write("\nFinal iteration clustering result:\n")
    f.write("Number of clusters generated: %i\n" % len(centers))
    f.write("Average distance to centers: %f\n" % np.mean(dists))
    f.write("STD of distance to centers: %f\n" % np.std(dists))
    f.write("List of cluster centers:\n" + 
            "\n".join(map(str,[struct.struct_id for struct in centers])) + "\n")
    f.write("Assigned cluster labels:\n")
    f.write("\n".join(map(str,[struct.struct_id + " " + 
                               str(struct.properties[stored_property_key])+" "+
                               str(dist_mat[coll.index(struct)]
                               [coll.index(centers[struct.properties[stored_property_key]])])
                               for struct in coll])))
    f.write("\n")
    f.close()

    if inst.has_option(sname,"output_folder"):
        misc.dump_collection_with_inst(inst,sname,coll)

    if inst.has_option(sname,"exemplar_output_folder"):
        op = inst.get(sname,"output_folder")
        inst.set(sname,"output_folder",inst.get(sname,"exemplar_output_folder"))
        misc.dump_collection_with_inst(inst,sname,centers)
        inst.set(sname,"output_folder",op) 


   
def affinity_propagation_fixed_clusters(inst):
    '''
    AP that explores the setting of preference in order to generate desired number of clusters
    '''
    sname = "affinity_propagation_fixed_clusters"
    dist_mat_path = inst.get(sname,"dist_mat_input_file")
    n_clusters = inst.get_eval(sname,"number_of_clusters")
    pref_range = inst.get_eval(sname,"preference_range")
    pref_iter = inst.get_with_default(sname,"pref_test_max_iters",10,eval=True)
    affinity_type = inst.get_with_default(sname,"affinity_type",["exponential",1],eval=True)
    damping = inst.get_with_default(sname,"damping",0.5,eval=True)
    convergence_iter = inst.get_with_default(sname,"convergence_iter",15,eval=True)
    max_iter = inst.get_with_default(sname,"max_iter",200,eval=True)
#    preference = inst.get_with_default(sname,"preference",None,eval=True)
    stored_property_key = inst.get_with_default(sname,"stored_property_key","AP_cluster")
    cluster_output_file = inst.get_with_default(sname,"cluster_output_file","./AP_cluster.info")
    
    coll = misc.load_collection_with_inst(inst,sname)
    coll.sort(key=lambda struct: struct.struct_id)
    
    f = open(dist_mat_path,"r")
    lines = f.read().split("\n")
    while lines[-1]=="":
        lines.pop()
    dist_mat = [[float(x) for x in y.split()] for y in lines]

    pref_l, pref_u = pref_range

    iter_n = 0
    while pref_iter > iter_n: #Iterate to search for the good preference that gives the amount of clusters
        pref = (pref_l + pref_u) / 2.0
        preference = _generate_preference_array(coll,["constant",pref])
    
        coll, centers = AP_distance_matrix(coll, dist_mat, affinity_type=affinity_type,
                                           damping=damping, convergence_iter=convergence_iter,
                                           max_iter=max_iter, preference=preference,
                                           stored_property_key=stored_property_key)

        labels = np.array([struct.properties[stored_property_key] for struct in coll])
     
        f = open(cluster_output_file,"a")

        if len(centers) < len(coll):
             f.write("Iteration %i; preference used: %f; clusters generated %i; silhouette score: %f\n" 
                    % (iter_n, pref, len(centers), 
                       silhouette_score(np.array(dist_mat),labels,metric="precomputed")))
        else:
             f.write("Iteration %i; preference used: %f; clusters generated %i\n" 
                    % (iter_n, pref, len(centers)))

        f.close()
           
        iter_n += 1


        if len(centers) == n_clusters:
            break

        if len(centers) > n_clusters:
            pref_u = pref
        elif len(centers) < n_clusters:
            pref_l = pref   


    dists = [dist_mat[i][coll.index(centers[coll[i].properties[stored_property_key]])]
                for i in range(len(coll))]

    f = open(cluster_output_file,"a")
  
    f.write("\nFinal iteration clustering result:\n")
    f.write(_print_ap_cluster_information(centers, dists, coll, stored_property_key, dist_mat))
    f.close()

    if inst.has_option(sname,"output_folder"):
        misc.dump_collection_with_inst(inst,sname,coll)

    if inst.has_option(sname,"exemplar_output_folder"):
        op = inst.get(sname,"output_folder")
        inst.set(sname,"output_folder",inst.get(sname,"exemplar_output_folder"))
        misc.dump_collection_with_inst(inst,sname,centers)
        inst.set(sname,"output_folder",op) 


def AP_distance_matrix(coll, dist_mat, affinity_type=["exponential",1],
                       damping=0.5, convergence_iter=15, max_iter=200,
                       preference=None, stored_property_key="AP_cluster"):
    '''
    Given a collection of structures and their distance matrix
    Conduct the Affinity Propagation on the collection
    '''
    affinity_mat = _convert_distance_to_affinity_mat(dist_mat, affinity_type)
    return AP_affinity_matrix(coll, affinity_mat, damping=damping,
                              convergence_iter=convergence_iter, 
                              max_iter=max_iter, preference=preference,
                              stored_property_key=stored_property_key)


def AP_affinity_matrix(coll, affinity_mat, damping=0.5, convergence_iter=15,
                       max_iter=200, preference=None, stored_property_key="AP_cluster"):
    '''
    Given a collection of structure and their affinity matrix
    Conduct the Affinity Propagation on the collection
    Returns a full list of collection with cluster assigned,
    as well as a list of exemplars
    '''
    ap = AffinityPropagation(damping=damping, max_iter=max_iter, 
                             convergence_iter=convergence_iter,
                             copy=True, preference=preference,
                             affinity="precomputed",verbose=False)

    
    result = ap.fit(affinity_mat)
    for x in range(len(coll)):
        coll[x].properties[stored_property_key] = result.labels_[x]

    centers = [coll[x] for x in result.cluster_centers_indices_]
    return coll, centers

def AP_vector_list(coll, vector_list, damping=0.5, convergence_iter=15,
                   max_iter=200, preference=None, stored_property_key="AP_cluster"):
    '''
    Given a collection of structure and their feature vectors
    Conduct the Affinity Propagation on the collection
    Returns a full list of collection with cluster assigned,
    as well as a list of exemplars
    '''
    ap = AffinityPropagation(damping=damping, max_iter=max_iter, 
                             convergence_iter=convergence_iter,
                             copy=True, preference=preference,
                             affinity="euclidean",verbose=False)

   


def _convert_distance_to_affinity_mat(dist_mat, affinity_type=["exponential",1]):
    '''
    Converts a given distance matrix to affinity matrix
    '''
    m = len(dist_mat); n = len(dist_mat[0])
    affinity_mat = [[0]*n for x in range(m)]

    for i in range(m):
        for j in range(n):
            if affinity_type[0]=="exponential":
                affinity_mat[i][j] = - np.exp(dist_mat[i][j]*affinity_type[1])
            elif affinity_type[0]=="power":
                affinity_mat[i][j] = - dist_mat[i][j]**affinity_type[1]

    return affinity_mat

def _generate_preference_array(coll, pref):
    '''
    Generate preference array according to preference
    '''
    if pref==None:
        return None
    if pref[0]=="constant":
        return [pref[1]]*len(coll)


def _print_ap_cluster_information(centers, dists, coll, stored_property_key, dist_mat):
    '''
    Generates standard cluster information
    '''
    st = ""
    st += "Number of clusters generated: %i\n" % len(centers)
    st += "Average distance to centers: %f\n" % np.mean(dists)
    st += "STD of distance to centers: %f\n" % np.std(dists)

    labels = np.array([struct.properties[stored_property_key] for struct in coll])
    llabels = list(labels)
    n_struct = [llabels.count(struct.properties[stored_property_key]) for struct in centers]
    
    st += "STD of number of structures in clusters: %f\n" % np.std(n_struct)

    try:
        sil_score = silhouette_score(np.array(dist_mat),labels,metric="precomputed") 
    except:
        sil_score = 1.0
    st += "Silhouette score: %f\n\n" % sil_score 


    st += "Center and number of structures of each cluster:\n" 

    labels = list(labels)
    for struct in centers:
         st += "Cluster %i: %s %i\n" % (struct.properties[stored_property_key],
                               struct.struct_id,
                               labels.count(struct.properties[stored_property_key]))


    st += "\n"
    st += "Assigned cluster label and distance to center:\n"
    st += "\n".join(map(str,[struct.struct_id + " " + 
                             str(struct.properties[stored_property_key])+" "+
                             str(dist_mat[coll.index(struct)]
                             [coll.index(centers[struct.properties[stored_property_key]])])
                             for struct in coll]))
    st += "\n"

    return st


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
    preference = inst.get_or_none(section, "preference", eval=True)
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
            max_iter=max_iter, preference=preference,
            property_key=property_key, output_file=output_file,
            exemplars_output_dir=exemplars_output_dir,
            exemplars_output_format=exemplars_output_format,
            **pool_operation_keywords)

class AffinityPropagationExecutor(PoolOperation):
    def __init__(self, dist_mat_input_file,
            affinity_type=["exponential",1], damping=0.5,
            convergence_iter=15, max_iter=200,
            preference=None,
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
        self._preference = preference
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
    
    def run_single(self):
        # No need of subdir if only one run
        self._enable_subdir_output = False
        result = self._run(self._preference, enable_output=True)
        self._print_result_verbose(result)

    def run_batch(self, preference_list, verbose_output=False):
        enable_output = False
        if not self._output_dir is None:
            self._enable_subdir_output = True
            enable_output = True

        for preference in preference_list:
            self._run_async(preference, enable_output)

        results = self.wait_for_all()
        self._print_results(results, verbose_output)

    def run_fixed_num_of_clusters(self, num_of_clusters):
        pass

    def run_with_target_silhouette_score(self, target_silhouette_score):
        pass

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
            "exemplar_indices": exemplar_indices,
            "exemplar_ids": exemplar_ids,
            "avg_distance": avg_distance,
            "std_distance": std_distance,
            "max_distance": max_distance,
            "silhouette_score": sil_score,
            "preference": preference}
    
    return coll, result_dict

