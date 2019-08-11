# Python default modules


from ibslib.io import import_structures,output_struct_dict
from ibslib.motif import calc_save_property,calc_and_separate

try: plot_GA_prop
except: from ibslib.analysis.plotting import plot_GA_prop


###############################################################################
# 1. Calculate motifs for a Structure directory and output the Structures     #
#      in a new separate directory with sub-directories for each motif        #
###############################################################################

#struct_dir = "/Users/ibier/Research/Results/Hab_Project/FUQJIK/4_mpc/Energy_Based/cross_025/Hab_Calcs/database_json_files"
#output_dir = "/Users/ibier/Research/Results/Hab_Project/motif_scripts/test_dir"
## Output formats may be "json" or "geometry"
#output_format = 'json'
## Suggested options for motif evaluation
#motif_kwargs = \
#    {
#     'supercell': 3, 
#     'include_negative': True, 
#     'num_mol': 12
#     }
#
#calc_and_separate(struct_dir, output_dir, motif_kwargs, file_format='json',
#                  overwrite=True)

###############################################################################
# 2. Calculate motifs for a Structure directory and output the Structures     #
#     with the motif as a property in the json file                           #
###############################################################################

#struct_dir = "/Users/ibier/Research/Results/Hab_Project/FUQJIK/8_mpc/Arjuna_RCD_025_testrun/database_json_files"
## If output dir is not provided then default behavior is to use struct_dir 
#output_dir = "/Users/ibier/Research/Results/Hab_Project/motif_scripts/FUQJIK/Arjuna_rcd025"
## Output formats may be "json" or "geometry"
#output_format = 'json'
## Suggested options for motif evaluation
#motif_kwargs = \
#    {
#     'supercell': 3, 
#     'include_negative': True, 
#     'num_mol': 12
#     }
#calc_save_property(struct_dir, motif_kwargs, output_dir=output_dir)


###############################################################################
# 3. Motif Plotting                                                           #
#      Much faster if you pre-compute the motifs using 2.                     #
###############################################################################


#struct_dir = "FUQJIK/Arjuna_rcd025"
## If provided, the file will be saved at the provided path
## If len of str is 0 then no file will be saved
#figname = 'FUQJIK/images/Arjuna_rcd025_IP.pdf'
#prop_key = 'Habmax'
#energy_key = 'energy'
## Number of Molecules Per Cell
#nmpc = 4
## True:  Indicating to use motif property. Leave this.
#motif = True
## True:  Motifs are colored by either IP or GAtor
## False: Motifs are colored for best clarity
#color_ip = False
##True:  IP structures are used
##False: IP structures are first separated and included in final plot
#use_ip = True
## True:  Only plots IP
## False: No change to plot
#ip_only = True
#
#plot_GA_prop(struct_dir, prop_key, energy_key, nmpc, 
#             experimental_parameters = [],
#             motif=motif, 
#             color_ip=color_ip, use_ip=use_ip, ip_only=ip_only,
#             figname=figname)
