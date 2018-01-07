'''
Created on Dec 19, 2013

@author: newhouse
'''
import os

from core.file_handler import cwd, output_file


def local_message(message, replica):
    out_file = os.path.join(cwd, str(replica) + '.out')
    data_file = open(out_file, 'a')
    data_file.write(str(message) + '\n')
    data_file.close()

def error(message, replica=None):
    if replica == None: r = ''
    else: r = str(replica) + ' '
    out_file = os.path.join(cwd, 'error.out')
    data_file = open(out_file, 'a')
    data_file.write(r + str(message) + '\n')
    data_file.close()

def reset_local(replica):
    out_file = os.path.join(cwd, str(replica) + '.out')
    data_file = open(out_file, 'w')
    data_file.write(str('') + '\n')
    data_file.close()
    
def move_to_shared_output(replica):
    local_out_file = os.path.join(cwd, str(replica) + '.out')
    if not os.path.exists(local_out_file): pass
    else: 
        d_file = open(local_out_file, 'r')
        contents_string = d_file.read()
        d_file.close()
        
        data_file = open(output_file, 'a')
        data_file.write('Replica: ' + str(replica) + str(contents_string) + '\n')
        data_file.close()
    reset_local(replica)
