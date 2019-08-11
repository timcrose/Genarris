# -*- coding: utf-8 -*-

import numpy as np
from matplotlib.pyplot import cm

from ibslib.io import read,write
#from ibslib.analysis.plotting import plot_lattice

struct_dir = "../diversity_analysis/BZOXZT/combined_pool_no_dups"
s = read(struct_dir)

# Show all info about plot_lattice arguments
# For detailed information about arguments, please look at the 
# help output. 
print(help(plot_lattice))

axis_title_fontsize=16

plot_lattice_kwargs = \
    {
      "struct_dict": s,
      "prop": "energy_tier1_relaxed",
      "figname": "",
      "colormap": cm.hot,
      "colormap_lim": (0.1,0.6),
      "relative_energy": True,
      "nmpc": 2,
      "experimental_cell": [],
      "experimental_plot_text_kwargs": 
          {
            "text": "X",
            'family': 'serif',
            'color':  'green',
            'weight': 'normal',
            'size': 28,
          },
      "figsize": (8,6),
      "scatter_kwargs":
         {
            "s": 100, 
         },
     "xlim": [], 
     "ylim": [], 
     "zlim": [],
     "xlabel_kwargs":
         {
           "xlabel": 'a ($\AA$)',
           "fontsize": axis_title_fontsize,
           "labelpad": 10
         },
     "ylabel_kwargs":
         {
           "ylabel": 'b ($\AA$)',
           "fontsize": axis_title_fontsize,
           "labelpad": 10
         },
     "zlabel_kwargs":
         {
           "zlabel": 'c ($\AA$)',
           "fontsize": axis_title_fontsize,
           "labelpad": 10
         },
     "cbar_title_kwargs":
         {
           "ylabel": 'Relative Energy, kJ mol$^{-1}$',
           "fontsize": 16,
           "rotation": 270,
           "labelpad": 30   
         },
     "tick_params_kwargs":
         {
            "axis": "both",
            "labelsize": 16,
            "width": 5,
            "length": 5,
        },
     "cbar_tick_params_kwargs":
         {
           "axis": "both",
           "labelsize": 16,
           "width": 2,
           "length": 5,
         },
     "cbar_ticks": [],
     "xticks": [],
     "yticks": np.arange(4,24,4),
     "zticks": []
    }

plot_lattice(**plot_lattice_kwargs)
