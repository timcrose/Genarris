'''
Created by Patrick Kilecdi on 01/24/2017

Module for selection based on clusters

'''
from utilities import misc
def cluster_based_selection(inst):
    '''
    Main module for selection
    '''
    sname = "cluster_based_selection"

    preference_property = inst.get(sname,"preference_property_key")
    cluster_property = inst.get(sname,"cluster_property_key")
    n_selected = inst.get_eval(sname,"number_selected")
    select_max = inst.get_boolean(sname,"select_max")
     
    coll = misc.load_collection_with_inst(inst,sname)
    if cluster_property!="None":
        labels = [struct.properties[cluster_property] for struct in coll]
    else:
        labels = [0 for struct in coll]

    preferences = [struct.properties[preference_property] for struct in coll]


    selected = select_by_clusters(coll, labels, preferences,
                                 n_selected, select_max)

    if inst.has_option(sname,"output_folder"):
        misc.dump_collection_with_inst(inst,sname,selected)

    return selected


def select_by_clusters(coll, labels, preferences, n_selected,
                       select_max = False):


    if n_selected > len(coll):
        raise ValueError("Not enough structures to be selected")

    groups = list(set(labels))
    groups.sort() #Grab and sort all labels
    group_counts = [labels.count(x) for x in groups]


    Coll = []
    for i in range(len(coll)):
        Coll.append((coll[i],labels[i],preferences[i]))

    if select_max:
        Coll.sort(key=lambda c: (group_counts[groups.index(c[1])], c[1], -c[2]))
    else:
        Coll.sort(key=lambda c: (group_counts[groups.index(c[1])], c[1], c[2]))

    
    G = groups[:]
    groups.sort(key=lambda g: (group_counts[G.index(g)],g))
    group_counts.sort()


    n_group = len(groups)
    left = n_selected
    selected = []
    count = 0 

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
    
    

    
