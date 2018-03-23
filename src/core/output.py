"""
If any part of this module is used for a publication please cite:

X. Li, F. Curtis, T. Rose, C. Schober, A. Vazquez-Mayagoitia, K. Reuter,
H. Oberhofer, and N. Marom "Genarris: Random Generation of Molecular Crystal 
Structures and Fast Screening with a Harris Approximation, ",
J. Chem. Phys., DOI: 10.1063/1.5014038; arXiv 1803.02145 (2018)
"""
'''
Created on Dec 19, 2013

@author: newhouse
'''
import os

from core.file_handler import cwd, output_file



__author__ = "Xiayue Li, Timothy Rose, Christoph Schober, and Farren Curtis"
__copyright__ = "Copyright 2018, Carnegie Mellon University and "+\
                "Fritz-Haber-Institut der Max-Planck-Gessellschaft"
__credits__ = ["Xiayue Li", "Luca Ghiringhelli", "Farren Curtis", "Tim Rose",
               "Christoph Schober", "Alvaro Vazquez-Mayagoita",
               "Karsten Reuter", "Harald Oberhofer", "Noa Marom"]
__license__ = "BSD-3"
__version__ = "1.0"
__maintainer__ = "Timothy Rose"
__email__ = "trose@andrew.cmu.edu"
__url__ = "http://www.noamarom.com"

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
