'''
Conducts K-mean clustering
'''

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import numpy as np
from utilities import misc

def k_mean_clustering(inst):
    '''
    Conducts simple k-mean clustering and post-clustering selection
    '''
    # TODO: Refact to use the util class, PoolOperation
    sname = "k_mean_clustering"
    info_level = inst.get_info_level(procedure="K_Mean_Clustering")
    structure_dir = inst.get(sname, "structure_dir")
    structure_dir_depth = inst.get_with_default(
            sname, "structure_dir_depth", 0, eval=True)
    structure_suffix = inst.get_with_default(sname, "structure_suffix", "")
    coll = misc.input_pool(structure_dir,
            structure_dir_depth=structure_dir_depth,
            structure_suffix=structure_suffix)

    
    coll = misc.load_collection_with_inst(inst,sname)
    n_clusters = inst.get_eval(sname,"n_clusters")
    #The number of clusters to form as well as the number of centroids to generate.

    max_iter = inst.get_with_default(sname,"max_iter",300,eval=True)
    #Maximum number of iterations of the k-means algorithm for a single run.

    n_init = inst.get_with_default(sname,"n_init",10,eval=True)
    #Number of time the k-means algorithm will be run with different centroid seeds.
    #The final results will be the best output of n_init consecutive runs in terms of inertia.


    init = inst.get_with_default(sname,"init","k-means++")
    #Method for initialization

    algorithm = inst.get_with_default(sname,"algorithm","auto")
    #K-means algorithm to use.

    tol = inst.get_with_default(sname,"tol",0.0001,eval=True)
    #Relative tolerance with regards to inertia to declare convergence

    n_jobs = inst.get_processes_limit(sname)
    # Number of parallel instances to run K-means

    vector_property_name = inst.get(sname,"vector_property_name")
    #Property name of the vector used for k-mean

    cluster_property_name = inst.get_with_default(sname,"cluster_property_name",
                                                        "k_mean_cluster_index")
    #The key under which the cluster will be stored

    distance_property_name = inst.get_with_default(sname,"distance_property_name",None)
    #The key under which the distance from centroid will be stored

    output_info_file = inst.get_with_default(sname,"output_info_file", None)

    output_dir = inst.get_or_none(sname, "output_dir")
    output_format = inst.get_with_default(sname, "output_format", "both")

    if info_level >= 1:
        print("K-Mean clustering called to cluster %i structures into %i groups" 
              % (len(coll),n_clusters))

    #Main Clustering step
    coll, labels, dist_list, centers, n_struct, dists, inertia, sil, output, kmean = \
    conduct_k_mean(coll, vector_property_name,
                   cluster_property_name, 
                   distance_property_name=distance_property_name,
                   n_clusters=n_clusters, max_iter=max_iter, 
                   n_init=n_init, init=init, 
                   algorithm=algorithm, tol=tol, n_jobs=n_jobs)  

    ######Beginning output#####

    if output_info_file!=None:
        f = open(output_info_file,"a")
        f.write(output)
        f.close()
    

    if info_level >= 1:
        print("K-Mean clustering completed")

    if not output_dir is None:
        misc.output_pool(coll, output_dir, output_format)
    
    pcs_enabled = inst.get_boolean(sname,"pcs_enabled")
    if not pcs_enabled:
        return coll
    
    ################## Begins Post-Clustering Selection ##################
    pcs_method = inst.get(sname,"pcs_method")
    n = inst.get_eval(sname,"pcs_number")
    pcs_reserved_variable = inst.get_with_default(sname,"pcs_reserved_variable",None)
    selected = post_clustering_selection(coll, n,
                                         cluster_property_name, pcs_method,
                                         reserved_variable=pcs_reserved_variable)

    pcs_output_dir = inst.get_or_none(sname, "pcs_output_dir")
    pcs_output_format = inst.get_with_default(
            sname, "pcs_output_format", "both")
    if not pcs_output_dir is None:
        misc.output_pool(selected, pcs_output_dir, pcs_output_format)

    if inst.has_option(sname,"pcs_output_info_file"):
        outs = "Diversity information for selected pool from directory %s\n" \
                % inst.get(sname,"structure_dir")
        result = []
        for i in range(len(selected)-1):
            for j in range (i+1,len(selected)):
                result.append(_calculate_distance(selected[i].properties[vector_property_name], 
                                                  selected[j].properties[vector_property_name]))
                outs += "%s %s %f\n" % (selected[i].struct_id, selected[j].struct_id, result[-1])

        outs += "Mean rdf vector difference: %f\n" % numpy.mean(result)
        outs += "Mean rdf vector difference std: %f\n" % numpy.std(result)
        if pcs_method == "property_max" or pcs_method == "property_min":
            property_list = [struct.properties[pcs_reserved_variable] for struct in selected]
            outs += "struct.properties[%s] mean: %f\n" % (pcs_reserved_variable,numpy.mean(property_list))
            outs += "struct.properties[%s] std: %f\n" % (pcs_reserved_variable,numpy.std(property_list))
    outs += "\n"
    write_log.write_log(inst.get(sname,"pcs_output_info_file"),outs)


    return coll, selected


def conduct_k_mean(coll, vector_property_name, 
                   cluster_property_name, distance_property_name=None,
                   n_clusters=8, max_iter=300, n_init=10,
                   init="k-means++", algorithm="auto", 
                   tol=0.0001, n_jobs=1):

    try:
        kmean = KMeans(n_clusters, max_iter=max_iter, n_init=n_init,
                       init=init, algorithm=algorithm, tol=tol,
                       n_jobs=n_jobs)
    except:
        kmean = KMeans(n_clusters, max_iter=max_iter, n_init=n_init,
                       init=init, tol=tol,
                       n_jobs=n_jobs)


    coll.sort(key=lambda struct: struct.struct_id)
    X = [struct.properties[vector_property_name] for struct in coll]

    kmean.fit(X)
    labels = [int(x) for x in kmean.labels_]
    centers = kmean.cluster_centers_
    inertia = kmean.inertia_

    dist_list = []

    n_struct = [0 for x in range(n_clusters)]
    dists = [[] for x in range(n_clusters)]



    for i in range(len(coll)):
    #Process cluster information for each structure
        coll[i].properties[cluster_property_name] = labels[i]
        n_struct[labels[i]] += 1

        dist = np.linalg.norm(np.subtract(coll[i].properties[vector_property_name],
                                          centers[labels[i]]))
        dist_list.append(dist)
        dists[labels[i]].append(dist)

        if distance_property_name!=None:
           coll[i].properties[distance_property_name] = dist

    sil = silhouette_score(np.array(X), np.array(labels), metric="euclidean")

    output = _generate_cluster_output(coll, labels, dist_list, centers, n_struct, dists, inertia, sil)

    return coll, labels, dist_list, centers, n_struct, dists, inertia, sil, output, kmean 


def _generate_cluster_output(coll, labels, dist_list, centers, n_struct, dists, inertia, sil):
    
    st = "Cluster generated: %i\n" % len(centers)
    st += "Inertia of final clustering: %f\n" % inertia
    st += "Silhouette score: %f\n" % sil
    st += "STD of cluster structure number: %f\n" % np.std(n_struct)
    
    st += "Number of structures, mean distance and std of distance of each cluster:\n"
    for i in range(len(centers)):
        st += "%i %i %f %f\n" % (i, n_struct[i], np.mean(dists[i]), np.std(dists[i]))
    st += "\n"

    st += "Cluster centers:\n"
    for i in range(len(centers)):
        st += " ".join(map(str,centers[i])) + "\n"

    st += "\n"

    st += "Cluster assigned and distance to center of each structure:\n"
    for i in range(len(coll)):
        st += "%s %i %f\n" % (coll[i].struct_id, labels[i], dist_list[i])

    return st

def _calculate_distance(v1, v2):
    return numpy.linalg.norm(numpy.subtract(v1, v2))
 

def post_clustering_selection(coll,n,cluster_property_name,method,reserved_variable=None):
    '''
    Selects a number of structures from the collection
    coll: collection
    n: number of structures to be selected
    cluster_property_name: the name (key) used for storing the clustering results
    methods:
    (1) property_min: Takes the min of the property within each cluster
        reserved_varaible should be set to the property name. 
    (2) property_max: Takes the max of the property within each cluster
        reserved_variable should be set to the property name.   
    
    '''
    def select_max(coll,property_name):
        p = [struct.properties[property_name] for struct in coll]
        return p.index(max(p))

    def select_min(coll,property_name):
        p = [struct.properties[property_name] for struct in coll]
        return p.index(min(p))
    
    def get_c(struct):
        return struct.properties[cluster_property_name]

    if n > len(coll):
        raise ValueError("Not enough structures in the collection to select %i structures" % n)

    selected = []
    cluster_indices = []
    cluster_struct_num = []
    struct_in_clusters = []
    original_n = n
    for struct in coll:
    #Puts all strctures into their clusters
        c = get_c(struct)
        if c not in cluster_indices:
            cluster_indices.append(c)
            struct_in_clusters.append([struct])
            cluster_struct_num.append(1)
        else:
            struct_in_clusters[cluster_indices.index(c)].append(struct)
            cluster_struct_num[cluster_indices.index(c)] += 1

    while math.ceil(float(n)/len(cluster_indices)) >= min(cluster_struct_num):
    #If certain clusters do not have enough structures within
    #to meet the average structure per cluster,
    #then the entire cluster will be selected
        c = cluster_struct_num.index(min(cluster_struct_num))
        selected += struct_in_clusters[c]
        n -= cluster_struct_num[c]
        cluster_indices.pop(c)
        struct_in_clusters.pop(c)
        cluster_struct_num.pop(c)
    m = int(math.ceil(float(n)/len(cluster_indices))+0.0001)
    #Now selects m from each cluster leftover
    for i in range(len(cluster_indices)):
        for j in range (m):
            if method == "property_min":
                l = select_min(struct_in_clusters[i],reserved_variable)
            elif method == "property_max":
                l = select_max(struct_in_clusters[i],reserved_variable)
            elif method == "random":
                l = random.randint(0,len(struct_in_clusters[i])-1)
            else:
                raise ValueError("Unsupported method of selection within cluster: %s" % method)
            selected.append(struct_in_clusters[i][l])
            struct_in_clusters[i].pop(l)

    if len(selected) >= original_n:
    #Selecting an average number of structures from each cluster may result in over selection
    #This part removes the additional structure from selected list
        if method == "property_min":
            pool_management.collection_sort(coll,keyword=reserved_variable)
        elif method == "property_max":
            pool_management.collection_sort(coll,keyword=reserved_variable,reverse=True)
        elif method == "random":
            selected = random.sample(selected,original_n)
        else:
            raise ValueError("Unsupported method of selection within Cluster")
        selected = selected[:original_n]

    return selected


