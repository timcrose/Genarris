# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-

from ibslib.io.io_structures import output_struct_dict
from ibslib.descriptor.find_nearest_structs import find_nearest_structs

import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from sklearn.decomposition import PCA

#struct_dir = 'all_relaxed_structures/all_json_cleaned'
struct_dir = 'original_data/Crossover_Testing/database_json_files'
experimetnal_struct_path = '../Experimental/FUQJIK.json'
output_dir = 'original_data/Crossover_Testing/Crossover_Testing_BrN'
file_format = 'json'
overwrite = True

###############################################################################
# Also want method for identifying best descriptor to use                     #
###############################################################################

# Settings for RDF
descriptor = 'rdf'
kwargs = {
          'atomic_pairs_list': ['Br','N'],
          'atomic_distance_range': [0,15],
          'atomic_distance_increment': 1,
          'smoothing_parameter': 1,
          'descriptor_key': 'rdf'
         }


###############################################################################
# New ibslib implemenetation                                                  #
###############################################################################

struct_dict,results = find_nearest_structs(experimetnal_struct_path,
                                           struct_dir,
                                           descriptor=descriptor,
                                           **kwargs)

output_struct_dict(output_dir,struct_dict,
                   file_format=file_format, overwrite=overwrite)

# Next figure
color = cm.rainbow(np.linspace(0,1,len(struct_dict)))
fig = plt.figure()

rdf_list = copy.deepcopy(results['nearest_descriptor'])
rdf_list.append(results['target_descriptor'])
rdf_array = np.array(rdf_list)

pca_obj = PCA(n_components=2)
pc_vector = pca_obj.fit_transform(rdf_array)

plt.scatter(pc_vector[0:-1,0],pc_vector[0:-1,1], c=color)
plt.scatter(pc_vector[-1,0],pc_vector[-1,1],c='k',s=100,alpha=0.33)

plt.show()

fig = plt.figure()
color = cm.rainbow(np.linspace(0,1,len(results['nearest_descriptor'])))
for i,descriptor_value in enumerate(results['nearest_descriptor']):
    plt.plot(descriptor_value, c=color[i])
    
plt.plot(results['target_descriptor'],c='k')
    
#plt.plot(results['target_descriptor'],c='k')
plt.show()

results_list = []
for i,struct_id in enumerate(results['nearest_id'][0:20]):
    energy = struct_dict[struct_id].get_property('energy_tier1_relaxed')
    results_list.append((struct_id,energy,i))
sorted(results_list, key=lambda x: x[1])