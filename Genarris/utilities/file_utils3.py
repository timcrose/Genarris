
"""
Created on Tue Feb  6 19:16:48 2018

@author: timcrose
"""

import csv, json, os, sys, pickle
from glob import glob
import shutil
from Genarris.utilities import time_utils

def output_from_rank(message_args, rank, mode='a', output_fpath_prefix='output_from_world_rank_'):
    output_fpath = output_fpath_prefix + str(rank)
    with open(output_fpath, mode=mode) as f:
          print(message_args, file=f)

def grep_dir_recursively(search_str, dir_path, read_mode, found_lines, return_list=True, search_from_top_to_bottom=True):
    from file_utils import grep_single_file
    for sub_path in glob(os.path.join(dir_path, '**'), recursive=True):
        if not os.path.isdir(sub_path):
            found_result = grep_single_file(search_str, sub_path, read_mode, found_lines, return_list=return_list, search_from_top_to_bottom=search_from_top_to_bottom)
            if return_list:
                found_lines += found_result
            elif found_result:
                return True
    if return_list:
        return found_lines
    else:
        return False
    
def write_row_to_csv(path, one_dimensional_list, mode='a', delimiter=','):
    if path[-4:] != '.csv':
        raise Exception('path must have .csv extension. path:', path)

    if type(one_dimensional_list) != list:
        raise TypeError('row is not type list, cannot write to csv. type(one_dimensional_list):', \
              type(one_dimensional_list), 'one_dimensional_list:', one_dimensional_list)
    if 'b' == mode[-1]:
        mode = mode[:-1]
    with open(path, mode, newline='') as f:
        csvWriter = csv.writer(f, delimiter=delimiter)
        csvWriter.writerow(one_dimensional_list)
        
def write_rows_to_csv(path, two_Dimensional_list, mode='w', delimiter=','):

    if 'b' == mode[-1]:
        mode = mode[:-1]
    f = open(path, mode, newline='')

    csvWriter = csv.writer(f, delimiter=delimiter)
    for row in two_Dimensional_list:
        if type(row) is not list:
            raise TypeError('row is not type list, cannot write to csv. The type of row is ' + str(type(row)), 'row', row)
        csvWriter.writerow(row)
    f.close()
