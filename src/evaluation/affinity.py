'''
Created by Patrick Kilecdi on 12/31/2016

Conducts Affinity Propagation for a given collection
'''

from sklearn.cluster import AffinityPropagation
from sklearn.metrics import silhouette_score
import numpy as np
from utilities import misc


def affinity_propagation_distance_matrix(inst):
    '''
    Main AP distance matrix module
    '''
    sname = "affinity_propagation_distance_matrix"
    dist_mat_path = inst.get(sname,"dist_mat_input_file")
    affinity_type = inst.get_with_default(sname,"affinity_type",["exponential",1],eval=True)
    damping = inst.get_with_default(sname,"damping",0.5,eval=True)
    convergence_iter = inst.get_with_default(sname,"convergence_iter",15,eval=True)
    max_iter = inst.get_with_default(sname,"max_iter",200,eval=True)
    preference = inst.get_with_default(sname,"preference",None,eval=True)
    stored_property_key = inst.get_with_default(sname,"stored_property_key","AP_cluster")
    cluster_output_file = inst.get_with_default(sname,"cluster_output_file","./AP_cluster.info")
    
    coll = misc.load_collection_with_inst(inst,sname)
    coll.sort(key=lambda struct: struct.struct_id)
    preference = _generate_preference_array(coll,preference)
    
    f = open(dist_mat_path,"r")
    lines = f.read().split("\n")
    while lines[-1]=="":
        lines.pop()
    dist_mat = [[float(x) for x in y.split()] for y in lines]
    
    coll, centers = AP_distance_matrix(coll, dist_mat, affinity_type=affinity_type,
                                       damping=damping, convergence_iter=convergence_iter,
                                       max_iter=max_iter, preference=preference,
                                       stored_property_key=stored_property_key)

    dists = [dist_mat[i][coll.index(centers[coll[i].properties[stored_property_key]])]
                for i in range(len(coll))]

    f = open(cluster_output_file,"w")
    f.write(_print_ap_cluster_information(centers, dists, coll, stored_property_key, dist_mat))
    f.close()


    if inst.has_option(sname,"output_folder"):
        misc.dump_collection_with_inst(inst,sname,coll)

    if inst.has_option(sname,"examplar_output_folder"):
        op = inst.get(sname,"output_folder")
        inst.set(sname,"output_folder",inst.get(sname,"examplar_output_folder"))
        misc.dump_collection_with_inst(inst,sname,centers)
        inst.set(sname,"output_folder",op) 

def affinity_propagation_analyze_preference(inst):
    '''
    Analyze the effect of affinity preference
    '''
    sname = "affinity_propagation_analyze_preference"
    dist_mat_path = inst.get(sname,"dist_mat_input_file")
    affinity_type = inst.get_with_default(sname,"affinity_type",["exponential",1],eval=True)
    damping = inst.get_with_default(sname,"damping",0.5,eval=True)
    convergence_iter = inst.get_with_default(sname,"convergence_iter",15,eval=True)
    max_iter = inst.get_with_default(sname,"max_iter",200,eval=True)
    pref_list = inst.get_eval(sname,"preference_list")
    stored_property_key = inst.get_with_default(sname,"stored_property_key","AP_cluster")
    cluster_output_file = inst.get_with_default(sname,"cluster_output_file",
                                                      "./AP_analyze_preference.info")
    
    coll = misc.load_collection_with_inst(inst,sname)
    coll.sort(key=lambda struct: struct.struct_id)

    f = open(dist_mat_path,"r")
    lines = f.read().split("\n")
    while lines[-1]=="":
        lines.pop()
    dist_mat = [[float(x) for x in y.split()] for y in lines]
    f.close()


    for pref in pref_list:
        preference = _generate_preference_array(coll,pref)
        f = open(cluster_output_file,"a")
    
        coll, centers = AP_distance_matrix(coll, dist_mat, affinity_type=affinity_type,
                                           damping=damping, convergence_iter=convergence_iter,
                                           max_iter=max_iter, preference=preference,
                                           stored_property_key=stored_property_key)

        dists = [dist_mat[i][coll.index(centers[coll[i].properties[stored_property_key]])]
                    for i in range(len(coll))]
        
        labels = np.array([struct.properties[stored_property_key] for struct in coll])
  
        try:
            sil_score = silhouette_score(np.array(dist_mat),labels,metric="precomputed") 
        except:
            sil_score = 1.0
        f.write(" ".join(map(str,pref))+" clusters: %i; mean distance: %f std distance: %f; max distance: %f; silhouette score: %f\n" 
                % (len(centers), np.mean(dists), np.std(dists), max(dists), sil_score))
        f.close()

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

    if inst.has_option(sname,"examplar_output_folder"):
        op = inst.get(sname,"output_folder")
        inst.set(sname,"output_folder",inst.get(sname,"examplar_output_folder"))
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

    if inst.has_option(sname,"examplar_output_folder"):
        op = inst.get(sname,"output_folder")
        inst.set(sname,"output_folder",inst.get(sname,"examplar_output_folder"))
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

    if inst.has_option(sname,"examplar_output_folder"):
        op = inst.get(sname,"output_folder")
        inst.set(sname,"output_folder",inst.get(sname,"examplar_output_folder"))
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
    as well as a list of examplars
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
    as well as a list of examplars
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

#    print(affinity_type)
#    print m, n
    for i in range(m):
        for j in range(n):
#            print i, j
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
 
