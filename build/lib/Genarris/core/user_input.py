"""
If any part of this module is used for a publication please cite:

X. Li, F. Curtis, T. Rose, C. Schober, A. Vazquez-Mayagoitia, K. Reuter,
H. Oberhofer, and N. Marom "Genarris: Random Generation of Molecular Crystal 
Structures and Fast Screening with a Harris Approximation, ",
J. Chem. Phys., DOI: 10.1063/1.5014038; arXiv 1803.02145 (2018)
"""
'''
Created on Nov 4, 2013

@author: newhouse

This module contains the default user input values.
Values are overridden by textual user input
'''
from ConfigParser import SafeConfigParser
import ast

from Genarris.core.file_handler import default_config, ui_conf

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



DEFAULT_CONFIG_REPLICA = -1

class ListSafeConfigParser(SafeConfigParser):
    '''Inherits SafeConfigParser and provides list parsing with json'''
    
    # TODO: maybe i could use literaleval(super.get()) instead, so to always return lists and ints
    def get_list(self, section, option):
        '''provides list parsing with json'''
        return self.get(section, option).split()
    
    def get_eval(self, section, option):
        return ast.literal_eval(self.get(section, option))
    
def get_config():
    '''
    Reads in default and user defined UI from the filesystem
    '''
    config = ListSafeConfigParser()

    default_config_file = open(default_config, 'r')
    config.readfp(default_config_file)
    default_config_file.close()

    local_config_file = open(ui_conf, 'r')
    config.readfp(local_config_file)
    local_config_file.close()
    return config
