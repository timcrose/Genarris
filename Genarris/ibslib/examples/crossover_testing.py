# -*- coding: utf-8 -*-

from ibslib.io.io_structures import import_structures,\
                                    output_struct, \
                                    output_geo
#from ibslib.crossover.standard_crossover import main as std_cross_main
import numpy as np

import sys
sys.path.append('/Users/ibier/Research/Side_Projects/ibslib/ibslib/crossover')
from standard_crossover import main as std_cross_main

struct_dir = 'structure_raw_pool'
replica = 'test'
nmpc = 8

output_path_struct = 'test_child_struct.json'
output_path_geo = 'test_child_struct.in'

#struct_dict = import_structures(struct_dir)
#struct_id_list = [x for x in struct_dict.keys()]
#struct_choice = np.random.choice(struct_id_list,2)

struct_choice = [struct_id_list[0], struct_id_list[1]]

struct_list = [struct_dict[x] for x in struct_choice]

child_struct = std_cross_main(struct_list,replica,nmpc)

output_geo(output_path_geo, child_struct)


'''
- Can make 4 mpc work easily, but 8 mpc seems to be elusive.
    I wonder why this is. I bet it is a simple fix. Either to how I am calling
    it or a change to the src code. 
- Look at the src to find what's going wrong. 
'''