# -*- coding: utf-8 -*-

from ibslib.analysis.plotting import plot_GAtor_Hab

struct_dir = 'database_json_files'
Hab_key = 'Habmax'
energy_key = 'energy'
nmpc = 8
experimental_parameters = ['FUQJIK', -1212815.055196140,
                           98.78]

plot_GAtor_Hab(struct_dir, Hab_key, energy_key, nmpc,
               experimental_parameters=experimental_parameters,
               figname='GAtor_Hab_FUQJIK_8mpc.pdf')