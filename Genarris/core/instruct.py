"""
If any part of this module is used for a publication please cite:

X. Li, F. Curtis, T. Rose, C. Schober, A. Vazquez-Mayagoitia, K. Reuter,
H. Oberhofer, and N. Marom "Genarris: Random Generation of Molecular Crystal 
Structures and Fast Screening with a Harris Approximation, ",
J. Chem. Phys., DOI: 10.1063/1.5014038; arXiv 1803.02145 (2018)
"""
'''
Created by Patrick Kilecdi (Xiayue Li) on 09/29/2015
This module stores the class Instruct,
which provides the instruction for structure generation
'''

try:
    from ConfigParser import SafeConfigParser
except ImportError:
    from configparser import SafeConfigParser
import ast, os
from Genarris.utilities import list_utils
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
from Genarris.utilities import file_utils


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

class Instruct(SafeConfigParser):
    '''
    This is the class that records the instruction for structural generation
    All of the information is set in the dictionary self.I
    '''
    def __init__(self):
        '''
        In initialization, sets the default (recommended) values for certain parameters
        '''
        SafeConfigParser.__init__(self)
        self.add_section("structure_generation")
        self.add_section("Genarris_master")
        self.add_section("FHI-aims")
        self.update_keywords([["Genarris_master","processes_limit","1"],["FHI-aims","execute_style","safe_subprocess"],["FHI-aims","runjob_modes","4"],["FHI-aims","runjob_thres","4"],["FHI-aims","multiple_launches","1"],["Genarris_master","info_level","1"]])

    def check_keywords(self,clist=[]):
        '''
        Checks whether or not the keys in the clist has been properly assigned to self.I (exists and non-void)
        Returns True if yes, the missing key if not
        Default clist includes: cell_vol, 
        '''
        for item in clist:
            if not self.has_option(item[0],item[1]):
                return item
        return True

    def get_boolean(self,section,option):
        '''
        Check and see if the section and option is specified. If so, use it if valid, otherwise raise an error. Return False by default.
        '''
        if self.has_option(section,option):
            value = self.get(section,option)
                
            if value == 'True' or value == 'true' or value == 'TRUE':
                return True
            elif value == 'False' or value == 'false' or value == 'FALSE':
                return False
            else:
                raise ValueError("Optional boolean flag has to be set to True, true, TRUE, False, false, or FALSE if present")

        # default to False
        return False

    def get_eval(self,section,option):
        try:
            return ast.literal_eval(self.get(section,option))
        except:
            return eval(self.get(section,option))

    def get_with_default(self,section,option,default,eval=False):
        self.set_default(section,option,str(default))
        if not eval:
            return self.get(section,option)
        else:
            return self.get_eval(section,option)

    def get_or_none(self, section, option, eval=False):
        if self.has_option(section, option):
            if eval:
                return self.get_eval(section, option)
            else:
                return self.get(section, option)
        else:
            return None

    def get_keywords(self,klist,eval=False):
        '''
        Retrieve a list of keywords from the configuration
        klist=[[necessary_keyword_section,necessary_keyword_option],[optional_keyword_section,optional_keyword,default_value]]
        get_keywords will not check for the existence of the necessary keyword before (will cause error if not exist)
        if eval=True, then will return call ast to do literal_eval
        '''
        result=[]
        for item in klist:
            if len(item)==2:
                if eval:
                    result.append(self.get_eval(item[0],item[1]))
                else:
                    result.append(self.get(item[0],item[1]))
            elif len(item)==3:
                try:
                    if eval:
                        result.append(self.get_eval(item[0],item[1]))
                    else:
                        result.append(self.get(item[0],item[1]))
                except ValueError:
                    raise ValueError("Malformed string when attempting to literal_eval an optional keyword")
                except:
                    result.append(item[2])
            else:
                raise ValueError("Unsupported entry for klist")
        return result

    def get_keywords_single_section(self,section,klist,eval=False):
        '''
        Retrive keywords from the same section
        The section in klist is omitted
        if any element is a single string, then considred necessary keyword
        if any element a list of two values, then considred [optional_keyword_name,default_value]
        '''
        result=[]
        for item in klist:
            if isinstance(item,str):
                if eval:
                    result.append(self.get_eval(section,item))
                else:
                    result.append(self.get(section,item))
            elif isinstance(item,list) and len(item)==2:
                try:
                    if eval:
                        result.append(self.get_eval(section,item[0]))
                    else:
                        result.append(self.get(section,item[0]))
                except ValueError:
                    raise ValueError("Malformed string when attempting to literal_eval an optional keyword")
                except:
                    result.append(item[1])
            else:
                raise ValueError("Unsupported entry for klist")
        return result
    
    def get_kv_single_section(self, section, keyword_list, eval=False):
        result = {}
        for key in keyword_list:
            if self.has_option(section, key):
                if eval:
                    result[key] = self.get_eval(section, key)
                else:
                    result[key] = self.get(section, key)
        return result



    def load_instruction_from_file(self,path):
        '''
        Reads from a path the instruction
        And update the keywords in self.I
        Interpretation of the configuration file is done with SafeConfigParser
        '''
        f = open(path,"r")
        self.readfp(f)
        f.close()

    
    def set_keywords_single_section(self,section,klist,clear_first=False):
        '''
        Receive a section and a list of [option, value]
        if clear_first, will remove and add the section beforehands
        '''
        if clear_first:
            self.remove_section(section)
            self.add_section(section)
        for option, value in klist:
            self.set(section,option,value)

    def set_default(self,section,option,value):
        '''
        Assign the value to [section,option] if the particular option does not exist
        '''
        if not self.has_section(section):
            self.add_section(section)
        if not self.has_option(section,option):
            self.set(section,option,value)
            return True
        return False
    
    def transfer_keywords(self,src_section,dst_section,klist):
        '''
        Transfer namesake options from one section to another
        '''
        for option in klist:
            if self.has_option(src_section,option):
                self.set(dst_section,option,self.get(src_section,option))

    def update_keywords(self,klist):
        '''
        Updates the configuration with klist
        '''
        for item in klist:
            self.set(item[0],item[1],str(item[2]))


    #Below is a list of functions to retrieve specific information
    # from an instruct object
    def procedure_exist(self,procedure):
        '''
        Checks whether or not a procedure exist
        '''
        return procedure in self.get_eval("Genarris_master","procedures")

    def has_procedure(self,procedure):
        '''
        Same as procedure_exist
        '''
        return procedure in self.get_eval("Genarris_master","procedures")   

    def set_procedures(self, procedures):
        if not type(procedures) is list:
            procedures = [procedures]
        self.set("Genarris_master", "procedures", str(procedures))

    def get_info_level(self,procedure=None):
        if procedure!=None and not self.has_procedure(procedure):
            return self.get_eval("Genarris_master","info_level")-1
        return self.get_eval("Genarris_master","info_level")
        
    def set_info_level(self,info_level):
        self.set("Genarris_master","info_level",str(info_level))

    def get_list(self,section,option):
        '''
        Get a list option that is potentially a string
        '''
        l = self.get(section,option)
        try:
            l = eval(l)
        except:
            pass
        if not isinstance(l,list):
            return [l]
        else:
            return l

    def get_master_working_dir(self):
        return self.get("Genarris_master","working_dir")
    
    def get_master_err_path(self):
        return self.get("Genarris_master","master_err_path")

    def get_tmp_dir(self, sname=None):
        if sname is None or not self.has_option(sname, "tmp_dir"):
            return os.path.abspath(self.get("Genarris_master","tmp_dir"))
        else:
            return os.path.abspath(self.get(sname, "tmp_dir"))
    
    def get_nodelist(self):
        return self.get_eval("parallel_settings","nodelist")

    def get_master_node(self):
        return self.get("Genarris_master","master_node")
    
    def get_master_processes_limit(self):
        return self.get("Genarris_master","processes_limit")

    def get_processes_limit(self, section):
        return self.get_with_default(section,"processes_limit",\
            self.get_master_processes_limit(),eval=True)

    def grant_permission(self,path):
        '''
        Grants permission to the given path by reading the Genarris_master.group_permission
        '''
        if self.has_option("Genarris_master","group_permission"):
            value = self.get("Genarris_master","group_permission")
            if value != "TRUE": 
                raise ValueError("Genarris_master.group_permission present but not set to TRUE; Make sure to enable a feature by setting the parameter to TRUE, all caps; else, omit the option")
            if os.path.isfile(path):
                os.system("chmod g=u "+path)
            elif os.path.isdir(path):
                os.system("chmod -R g=u "+path)
            else:
                raise ValueError("Not a file or directory: "+path)

    def make_path_absolute(self, section, option):
        if self.has_option(section, option):
            path = self.get(section, option)
            self.set(section, option, os.path.abspath(path))

    def write_to_file(self, file_path):
        f = open(file_path, "w")
        self.write(f)
        f.close()

    def __deepcopy__(self,memo):
        '''
        Due to the inability to deepcopy a configuration file
        Will generate a temporary config file and read it back in
        '''
        config_string = StringIO()
        self.write(config_string)
        config_string.seek(0)
        copied = Instruct()
        copied.readfp(config_string)
        return copied

    def get_inferred(self, sname, sname_list, options, default=None, eval=False, type_=str, required=True):
        '''
        sname: str
            Current section name
        sname_list: list of str or str
            list of section names to check. The first section that has option defined
            will get used.
        options: list of str or str
            list of option names that correspond to the list of snames given
        default: any
            default value for option
        eval: bool
            True if you want to evaluate the expression with python's eval function
        type_: any
            type of the desired output.
        required: bool
            True: will throw error if option DNE in sname_list
            False: will return "none_value" if option DNE in sname_list

        return: option value

        Purpose: The option might have been specified in a different (usually previous
            in a workflow of sections) section so we'd like to be able to grab that
            option value so the user doesn't have to specify redundant options in
            all sections. This function iterates over a list of section names and 
            returns the option value that exists in the first of these sections.
        '''
        if type(sname_list) is str:
            sname_list = [sname_list]
        if type(options) is str:
            options = [options] * len(sname_list)
        if len(options) != len(sname_list):
            raise Exception('Got different number of options and section names')
        found_sname = None
        for sname, option in zip(sname_list, options):
            if self.has_option(sname, option):
                if (type_ == 'dir' and os.path.isdir(self.get(sname, option))) or \
                    (type_ == 'file' and os.path.isfile(self.get(sname, option))) or \
                    (type_ != 'dir' and type_ != 'file'):

                    found_sname = sname
                    break
        if found_sname is None and (required or sname in sname_list):
            raise Exception('could not find options ' + ', '.join(options) + ' in any of the '+
            'following sections: ' + ', '.join(sname_list), 'Or, path didnt exist')
        if default is not None:
            value = self.get_with_default(found_sname,option,default,eval=eval)
        elif type_ is str or type_ == 'dir' or type_ == 'file':
            value = self.get(found_sname, option)
        elif type_ is bool:
            value = self.get_boolean(found_sname, option)
        elif type_ is list:
            value = self.get_list(found_sname, option)
        elif type_ is None:
            value = self.get_or_none(found_sname, option, eval=eval)
        elif eval:
            value = self.get_eval(found_sname, option)
        else:
            value = 'none_value'
        if required and value == 'none_value':
            raise Exception('Could not find options', options, 'in', sname_list)
        elif sname in sname_list and value == 'none_value':
            raise Exception('Could not find options', options, 'in', sname_list)
        
        if required and \
            (type_ == 'dir' and not os.path.isdir(value)) or \
            (type_ == 'file' and not os.path.isfile(value)):

            raise Exception('type_', type_, 'in', options, 'in', sname_list, 'DNE')

        elif sname in sname_list and \
            (type_ == 'dir' and not os.path.isdir(value)) or \
            (type_ == 'file' and not os.path.isfile(value)):

            raise Exception('type_', type_, 'in', options, 'in', sname_list, 'DNE')            

        return value

def get_random_index(seed=None):
        '''
        Outputs a random index
        '''
        from hashlib import sha1
        import time
        LENGTH_OF_INDEX = 10
        return sha1(repr(time.time())+str(seed)).hexdigest()[:LENGTH_OF_INDEX]

def get_procedure_name_from_section_name(sname):
    '''
    sname: str
        section name

    Return: str
        the procedure name from a given sname

    Purpose: provide the procedure name from a given section name
    Notes: Requires the format of procedure names to be the capital letter
        of every portion of sname that is delimited by an underscore. The exception is
        the portion "fhi" which should be "FHI" in its procedure name.
    '''
    return '_'.join(['FHI' if portion == 'fhi' else portion.capitalize() for portion in sname.split('_')])


def get_last_active_procedure_name(inst, sname, iteration=0):
    '''
    inst: Instruct
    sname: str
        current section name
    iteration: int
        The number of times this section has showed up previously in the active procedure list
        in the conf file. This allows sections to be done multiple times such as AP.
    
    Return: str
        The procedure name coming before the current procedure as specified in the
        active procedure list. If this procedure is first, then return 'none'. The procedure
        name is made all lowercased so it is a section name.
    '''
    procedure_name = get_procedure_name_from_section_name(sname)
    procedures = inst.get_eval('Genarris_master', 'procedures')
    
    procedure_idx = list_utils.indices(procedures, procedure_name)[iteration]
    if procedure_idx == 0:
        return 'none'
    else:
        return procedures[procedure_idx - 1].lower()

    
def get_molecule_path(inst, sname):
    if inst.has_option(sname, 'molecule_path'):
        molecule_path = inst.get(sname, 'molecule_path')
        return molecule_path
    elif inst.has_section('relax_single_molecule'):
        molecule_path_list = file_utils.find(os.path.abspath(inst.get('relax_single_molecule', 'aims_output_dir')), 'geometry.in.next_step')
        if len(molecule_path_list) == 1:
            molecule_path = molecule_path_list[0]
            return molecule_path
        else:
            last_section = get_last_active_procedure_name(inst, sname)
            if last_section == 'relax_single_molecule':
                raise Exception('Was not able to infer molecule_path. We do not want to use the unrelaxed geometry here.')
    last_section = get_last_active_procedure_name(inst, sname)
    sname_list = [sname, last_section, 'estimate_unit_cell_volume', 'harris_single_molecule_prep', 'pygenarris_structure_generation', 'structure_generation_batch', 'harris_approximation_batch']
    molecule_path = inst.get_inferred(sname, sname_list, ['molecule_path'] * len(sname_list), type_='file')
    return molecule_path