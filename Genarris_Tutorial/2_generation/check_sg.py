
import os
import json
import numpy as np

'''
Purpose: 
    Find all unique space groups from either a list of dictionaries produced by the AJ_Extractor or
    Find all unique space groups from a directory of json files
Usage:
    1. Either populate the dict_name_list or json_dir variables using relative paths to the current working directory
    2. Run script which prints the unique sg_list
'''

def main():
    dict_name_list = None
    json_dir = 'structure_raw_pool'
    cwd = os.getcwd()
    sg_list = collect_unique_sg(cwd, dict_name_list, json_dir)
    print(sg_list)

def collect_unique_sg(cwd, dict_name_list=None, json_dir=None):
    '''
    - Accepts either a list of dictionary names or a json directory
    '''
    sg_list = []
    if dict_name_list != None:
        for name in dict_name_list:
            dict_path = os.path.join(cwd, name)
            struct_dict = load_json(dict_path)

            for member in struct_dict:
                sg = struct_dict[member]['space_group']
            if sg not in sg_list:
                sg_list.append(sg)
    elif json_dir != None:
        json_dir_path = os.path.join(os.path.join(cwd, json_dir))
        for struct in os.listdir(json_dir_path):
            struct_path = os.path.join(json_dir_path, struct)
            if '.json' in struct_path:
                struct_dict = load_json(struct_path)
            else:
                continue
            sg = struct_dict['properties']['space_group']
            if sg not in sg_list:
                sg_list.append(sg)
    else: 
        print('Please supply dict_list or json_dir paths')
        return
    
    sg_list = np.sort(sg_list)
    return sg_list

def load_json(json_path):
    with open(json_path) as f:
        struct_dict = json.load(f)
    return struct_dict
           

if __name__ == '__main__':
    main()
