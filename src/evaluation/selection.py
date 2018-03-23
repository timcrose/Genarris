"""
If any part of this module is used for a publication please cite:

X. Li, F. Curtis, T. Rose, C. Schober, A. Vazquez-Mayagoitia, K. Reuter,
H. Oberhofer, and N. Marom "Genarris: Random Generation of Molecular Crystal 
Structures and Fast Screening with a Harris Approximation, ",
J. Chem. Phys., DOI: 10.1063/1.5014038; arXiv 1803.02145 (2018)
"""
'''
Created by Patrick Kilecdi on 01/24/2017

Module for selection based on clusters

'''
from evaluation.evaluation_util import load_pool_operation_keywords, \
        PoolOperation

__author__ = "Xiayue Li, Timothy Rose, Christoph Schober, and Farren Curtis"
__copyright__ = "Copyright 2018, Carnegie Mellon University and "+\
                "Fritz-Haber-Institut der Max-Planck-Gessellschaft"
__credits__ = ["Xiayue Li", "Luca Ghiringhelli", "Farren Curtis", "Tim Rose",
               "Christoph Schober", "Alvaro Vazquez-Mayagoita",
               "Karsten Reuter", "Harald Oberhofer", "Noa Marom"]
__license__ = "BSD-3"
__version__ = "1.0"
__maintainer__ = "Timothy Rose"
__email__ = "trose@andrew.cmu.edu"
__url__ = "http://www.noamarom.com"

def cluster_based_selection(inst):
    '''
    Main module for selection
    '''
    sname = "cluster_based_selection"

    number_selected = inst.get_eval(sname, "number_selected")
    select_max = inst.get_boolean(sname, "select_max")
    preference_property_name = inst.get_or_none(
            sname, "preference_property_name")
    cluster_property_name = inst.get_or_none(sname, "cluster_property_name")

    kwargs = load_pool_operation_keywords(inst, sname)
    kwargs["processes_limit"] = 1
    kwargs["enable_subdir_output"] = False

    executor = PoolOperation(**kwargs)

    return executor.run_operation(_cluster_based_selection,
            name="ClusterBasedSelection",
            args=(number_selected, select_max),
            kwargs={"preference_property_name": preference_property_name,
                    "cluster_property_name": cluster_property_name},
            enable_output=True)

def _cluster_based_selection(coll, number_selected, select_max,
        preference_property_name=None, cluster_property_name=None):
    if not cluster_property_name is None:
        labels = [struct.properties[cluster_property_name]
                  for struct in coll]
    else:
        labels = [0 for struct in coll]

    if not preference_property_name is None:
        preferences = [struct.properties[preference_property_name]
                       for struct in coll]
    else:
        preferences = [random.random() for struct in coll]

    return select_by_clusters(coll, labels, preferences, number_selected,
            select_max=select_max)

def select_by_clusters(coll, labels, preferences, n_selected,
                       select_max = False):
    if (len(coll) != len(labels)) or (len(coll) != len(preferences)):
        raise ValueError("Length of collection does not equal to"
                "length of labels or preferences")

    if n_selected > len(coll):
        raise ValueError("Not enough structures to be selected")

    groups = list(set(labels))
    groups.sort() #Grab and sort all labels
    group_counts = [labels.count(x) for x in groups]

    Coll = []
    for i in range(len(coll)):
        Coll.append((coll[i],labels[i],preferences[i]))

    # Sort collection by: group count, label name, and preference level
    if select_max:
        Coll.sort(key=lambda c: (group_counts[groups.index(c[1])], c[1], -c[2]))
    else:
        Coll.sort(key=lambda c: (group_counts[groups.index(c[1])], c[1], c[2]))

    # Sort groups by group count
    G = groups[:]
    groups.sort(key=lambda g: (group_counts[G.index(g)],g))
    group_counts.sort()

    n_group = len(groups)
    left = n_selected
    selected = []
    count = 0 

    # Select structures first from smaller groups
    for i in range(n_group):
        if group_counts[i] < int(left/(n_group-i)):
            
            #The entire group is selected
            selected += [c[0] for c in Coll[count:count+group_counts[i]]]
            left -= group_counts[i]
        
        else:
            
            selected += [c[0] for c in Coll[count:count+int(left/(n_group-i))]]
            left -= left/(n_group-i)

        count += group_counts[i]

    return selected
    
    

    
