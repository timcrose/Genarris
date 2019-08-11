# -*- coding: utf-8 -*-

from ibslib.io.io_structures import import_structures
#from ibslib.analysis.plotting import plot_many_property_vs_energy, \
#                                     plot_many_property_vs_energy_dict

import sys
sys.path.append('/Users/ibier/Research/Side_Projects/ibslib/ibslib/analysis/')
from plotting import plot_many_property_vs_energy

struct_dir_list = ['/Users/ibier/Research/Results/Hab_Project/FUQJIK/8_mpc/Arjuna_RCD_025_testrun/Initial_Pool/json',
                   '/Users/ibier/Research/Results/Hab_Project/FUQJIK/4_mpc/Property_Based/Initial_Pool/json',
                   '/Users/ibier/Research/Results/Hab_Project/DUXYUQ/Hab_Results/Initial_Results_20180819/init_json_Hab',
                   '/Users/ibier/Research/Results/Hab_Project/BZOXZT/Genarris/combined_IP/initial_Hab_jsons',
                   '/Users/ibier/Research/Results/Hab_Project/GIYHUR/GAtor/IP_37/GIYHUR_IP_37'
                   ]
prop_key_list = ['Habmax','Habmax','Hab','Habmax','Habmax']
energy_key_list = ['energy','energy','energy','energy','energy']
legend_list = ['FUQJIK','FUQJIK01','DUXYUQ','BZOXZT','GIYHUR']
# Order is energy and then property
experimental_data = [(-1212815.055196140,98.78),
                     (-606407.696428319,18.26),
                     (-69842.496704337,10.5),
                     (-43514.036686496,25.39),
                     (-35254.887115816,28.86)]
nmpc_list = [8,4,4,2,2]

kwargs = {
          'marker_list': ['o','^','s','h','D'],
          'color_list': ['tab:blue', 'tab:orange', 'tab:green', 
                         'tab:red',   'tab:pink', 'tab:purple',
                         'tab:cyan', 'tab:gray', 'tab:brown',
                         'tab:olive'],
          'xlabel': '$H_{ab}^{max}$ (meV)',
          'label_size': 32,
          'legend_columns': 2,
          'legend_font_size': 16,
          'figname': 'all_energy_Hab_20190106.pdf',
          'labelpad_x': 5,
          'labelpad_y': 5,
          'experimental_data': experimental_data
        }
#struct_dict_list = []
#for struct_dir in struct_dir_list:
#    struct_dict = import_structures(struct_dir)
#    struct_dict_list.append(struct_dict)
    
struct_dict_list = plot_many_property_vs_energy(struct_dir_list, prop_key_list, 
                             energy_key_list, legend_list, nmpc_list,
                             **kwargs)

#plot_many_property_vs_energy_dict(struct_dict_list, prop_key_list, 
#                             energy_key_list, legend_list, nmpc_list)