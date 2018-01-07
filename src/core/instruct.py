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

try:
	from StringIO import StringIO
except ImportError:
	from io import StringIO


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
		self.update_keywords([["structure_generation","bravais_system",-1],["structure_generation","p_tolerance",0.5],["structure_generation","angle_range",[60,120]],\
		["structure_generation","space_group","-1"],["structure_generation","lattice_variance",[0.6,2]],\
		["Genarris_master","processes_limit","1"],["FHI-aims","execute_style","safe_subprocess"],["FHI-aims","runjob_modes","4"],["FHI-aims","runjob_thres","4"],["FHI-aims","multiple_launches","1"],["Genarris_master","info_level","1"]])

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
		Check and see if the section and option is specified
		if specified, has to set to "TRUE", else an error will be raised
		'''
		if self.has_option(section,option):
			if self.get(section,option)!="TRUE":
				raise ValueError("Optional boolean flag has to be set to TRUE if present")
			return True
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


	#Below is a list of functions to retrieve specific information from an instruct object
	
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


def get_random_index(seed=None):
        '''
        Outputs a random index
        '''
        from hashlib import sha1
        import time
        LENGTH_OF_INDEX = 10
        return sha1(repr(time.time())+str(seed)).hexdigest()[:LENGTH_OF_INDEX]

