'''
Created on June 17, 2015

'''
import errno
import os
from shutil import rmtree
import time
from hashlib import sha1

#from external_libs.filelock import FileLock

cwd = os.getcwd()
def argument_opt():
        from optparse import OptionParser
        parser=OptionParser()
        parser.add_option('-d','--workdir',action='store',type='string',dest='working_directory',default='ui.conf',help='User input file name (default="ui.conf"')
        return parser.parse_args()


# source directories
src_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
GA_dir = os.path.abspath(os.path.join(src_dir, os.pardir))
res_dir = os.path.join(GA_dir, 'res')

# working directories TODO: make this movable
tmp_dir = os.path.join(cwd, 'tmp/')

# filesystem storage
structure_dir = os.path.join(cwd, 'structures')

# useful files
progress_file = os.path.join(cwd, tmp_dir , 'progress.dat')
default_config = os.path.join(GA_dir, 'res', 'default.conf')
(options,argv)=argument_opt() #Retrieve the user_input from command line
ui_conf = os.path.join(cwd, options.user_input)
replica_file = os.path.join(tmp_dir, 'replica_index.dat')
output_file = os.path.join(cwd, 'output.out')

# constants
INITIAL_POOL_REFID = -1

def mkdir_p(path):
    '''
    makes full directory path
    '''
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise
        
def mkdir_p_clean(path):
    '''
    removes directory and recreates it
    '''
    if os.path.exists(path): rmtree(path)
    mkdir_p(path)
    
def my_import(name, package=''):
    '''
    dynamically (at runtime) imports modules portentially specified in the UI
    taken from http://function.name/in/Python/__import__ 
    '''
    name = package + '.' + name
    mod = __import__(name)
    components = name.split('.')
    for comp in components[1:]:
        mod = getattr(mod, comp)
    return mod 

def write_data(filepath, filename, contents):
    if not filename == None: filepath = os.path.join(filepath, filename)
    d_file = open(filepath, 'w')
    d_file.write(str(contents))
    d_file.close()

def read_data(filepath, filename=None):
    if filename is not None: full_filepath = os.path.join(filepath, filename)
    else: full_filepath = filepath
    d_file = open(full_filepath, 'r')
    contents_string = d_file.read()
    d_file.close()
    return contents_string
    
def get_progress():
    with FileLock(progress_file):
        p_file = open(progress_file, 'r')
        progress_string = p_file.readline()
        p_file.close()
    return progress_string

def set_progress(progress_string):
    with FileLock(progress_file):
        p_file = open(progress_file, 'w')
        p_file.write(progress_string)
        p_file.close()
    return True

def get_index(index_path):
    if not os.path.exists(index_path):
        data_file = open(index_path, 'w')
        data_file.write(str(0))
        data_file.close()
        
    with FileLock(index_path):
        index_file = open(index_path, 'r')
        index = int(index_file.readline())
        index_file.close()
    return index

def get_and_increment_index(index_path):
    '''
    retrieves the next valid index and increments the index in the shared file
    uses FileLock to avoid conflicts 
    '''
    
    # create if does not exist
    if not os.path.exists(index_path):
        mkdir_p(os.path.abspath(os.path.join(index_path, os.pardir)))
        data_file = open(index_path, 'w')
        data_file.write(str(0))
        data_file.close()
        
    with FileLock(index_path):
        index_file = open(index_path, 'r')
        index = int(index_file.readline())
        index_file.close()
    
        data_file = open(index_path, 'w')
        data_file.write(str(index + 1))
        data_file.close()
    return index

def get_random_index(seed=None):
    LENGTH_OF_INDEX = 10
    return sha1(repr(time.time())+str(seed)).hexdigest()[:LENGTH_OF_INDEX]

def print_to_file(message):
    with FileLock(output_file):
        data_file = open(output_file, 'a')
        data_file.write(str(message) + '\n')
        data_file.close()

if __name__ == '__main__':
    print cwd
