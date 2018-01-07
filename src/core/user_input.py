'''
Created on Nov 4, 2013

@author: newhouse

This module contains the default user input values.
Values are overridden by textual user input
'''
from ConfigParser import SafeConfigParser
import ast

from core.file_handler import default_config, ui_conf


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
