import os
import json
import sys

'''
- Takes a directory of json files and extracts their geometry
- structure_dir is the directory of json files
- geometry_dir is the output directory for the geometries
- For now, run this script is the directory with both of these subdirectories
'''

structure_dir = "initial_pool"
geometry_dir = "batch_eval"
sys.path.append(structure_dir)


files = [x for x in os.listdir(structure_dir)]
if not os.path.exists(geometry_dir):
        os.makedirs(geometry_dir)

for x in files:
    data = json.load(open(os.path.join(structure_dir, x), 'r'))
    root_name = x.split('.') 	# Splits the name into before and after the period
    output_geometry = open(os.path.join(geometry_dir,root_name[0] + '_geometry'),'a+') # a+ makes read/write 
    output_geometry.write('lattice_vector ')    
    for i in data["properties"]["lattice_vector_a"]:
        output_geometry.write(str(i) + ' ')
    output_geometry.write('\n' + 'lattice_vector ')
    for i in data["properties"]["lattice_vector_b"]:
        output_geometry.write(str(i) + ' ')
    output_geometry.write('\n' + 'lattice_vector ')
    for i in data["properties"]["lattice_vector_c"]:
        output_geometry.write(str(i) + ' ')
    output_geometry.write('\n')
    for i in data["geometry"]:
        output_geometry.write('atom ' + str(i[0]) + ' ' + str(i[1]) + ' ' + str(i[2]) + ' ' + str(i[3]) + '\n')
