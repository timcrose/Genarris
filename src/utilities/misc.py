'''
Miscellaneous utilities functions
'''

from core import structure
from utilities import write_log
from external_libs.filelock import FileLock
from external_libs import matrix_op
import random
import copy
import numpy
import numpy as np
import os, shutil, time

def random_rotation_matrix():
	'''
	Generates and returns a random rotation matrix
	'''
	u=random.uniform(0,1)
	v=random.uniform(0,1)
	theta=360*u
	phi=numpy.rad2deg(numpy.arccos(2*v-1))
	ovec=[numpy.sin(numpy.deg2rad(phi))*numpy.cos(numpy.deg2rad(theta)),numpy.sin(numpy.deg2rad(phi))*numpy.sin(numpy.deg2rad(theta)),numpy.cos(numpy.deg2rad(phi))]
	rot=random.random()*360
	c=numpy.cos(numpy.deg2rad(rot));s=numpy.sin(numpy.deg2rad(rot))
	x,y,z=ovec
	return ([[x*x*(1-c)+c,x*y*(1-c)-z*s,x*z*(1-c)+y*s],
		[x*y*(1-c)+z*s,y*y*(1-c)+c,y*z*(1-c)-x*s],
		[x*z*(1-c)-y*s,y*z*(1-c)+x*s,z*z*(1-c)+c]])

def get_random_index(seed=None):
	'''
	Outputs a random index
	'''
	from hashlib import sha1
	import time
	LENGTH_OF_INDEX = 10
	return sha1((repr(time.time())+str(seed)).encode()).hexdigest()[:LENGTH_OF_INDEX]

def struct_id_assign(struct,struct_id_scheme,struct_index=0):
	'''
	Assigns struct_id according to the specified struct_id_scheme
	'''
	name = ""
	for item in struct_id_scheme:
		if item == "struct_index":
			name += str(struct_index)
		elif item == "random_index":
			name += get_random_index()
		elif item in struct.properties:
			name += str(struct.properties[item])
		else:
			name += item
	struct.struct_id = name
	return name


def frac_coor_adjust(vec,create_duplicate=True):
	'''
	Adjust the fractional coordinates to all have values between 0 to 0..9999
	'''
	if create_duplicate:
		vec = copy.deepcopy(vec)
	for j in range (3):
		if (vec[j] < -0.0001):
			vec[j] += int(-vec[j]+1)
		elif (vec[j] > 0.99999):
			vec[j] -= int(vec[j] + 0.00001)
	return vec

def molecule_harris_info_prep(trans,orien):
	'''
	Generate the molecule_harris_info character of the structure given by trans and orien
	'''
	if round(numpy.linalg.det(orien)) == -1:

		is_enantiomer = True
		orien = numpy.dot(orien,[[1,0,0],[0,1,0],[0,0,-1]])
	else:
		is_enantiomer = False
	angles =  list(numpy.rad2deg(matrix_op.euler_from_matrix(orien,"rzyz")))
	for k in range (3):
		if angles[k]<0:
			angles[k] += 360
	return angles+list(trans)+[is_enantiomer]


def load_structure_with_inst(inst,section,option="structure_format",structure_path=None):
	'''
	Read the structure from a specific structure path given
	The function will try to read the format as specified by inst.get(section,option)
	If the option does not exist, will try both json and geometry parsing
	Raise ValueError if the option exists but is not geometry or json
	Return structure if successfully loaded; return False if not
	'''
	if structure_path == None:
		inst.get(section,"structure_path")
	f = open(structure_path,"r")
	st = f.read()
	f.close()
	struct = structure.Structure()
	if not inst.has_option(section,option):
		try:
			struct.loads(st)
		except:
			try:
				struct.build_geo_whole_atom_format(st)
			except:
				return False
	else:
		format = inst.get(section,option)
		if format == "geometry":
			try:
				struct.build_geo_whole_atom_format(st)
				struct.properties["original_structure_directory"] = structure_path
			except:
				return False
		elif format == "json":
			try:
				struct.loads(st)
				struct.properties["original_structure_directory"] = structure_path
			except:
				return False
                else:
                    return False
	if struct.struct_id == None:
		if inst.has_option(section,"structure_suffix"):
			struct.struct_id = os.path.split(structure_path)[-1][:-len(inst.get(section,"structure_suffix"))]
		else:
			struct.struct_id = os.path.split(structure_path)[-1]
	
	return struct

def load_collection_with_inst(inst,section):
	'''
	Loads a collection of structure from inst according to the parameters set in the given section
	Relavant parameters: structure_dir, structure_dir_depth, structure_suffix, structure_format
	'''
	struct_files = get_structure_list(inst,section,False)
	return [load_structure_with_inst(inst,section,structure_path=path) for path in struct_files]

def dump_collection_with_inst(inst,section,coll):
	for struct in coll:
		dump_structure_with_inst(inst,section,struct)	

def dump_structure_with_inst(inst,section,struct):
	'''
	Reads the options in section and output the structure
	'''
	if not inst.has_option(section,"output_format"):
		write_log.write_master_err(inst,"dump_structure_with_inst missing output_format in section "+section)
		return False
#		raise ValueError("dump_structure_with_inst missing output_format in section "+section)
	if not inst.has_option(section,"output_file") and not inst.has_option(section,"output_folder"):
		write_log.write_master_err(inst,"dump_structure_with_inst missing both output_file and output_folder in section "+section)
		return False
 #               raise ValueError("dump_structure_with_inst missing both output_file and output_folder in section "+section)

	output_formats = inst.get_eval(section,"output_format")
	filename = get_random_index()
	for output_format in output_formats:
		if output_format != "json" and output_format != "geometry":
			write_log.write_master_err_and_raise(inst,"dump_structure_with_inst unknown output_format"+output_format)

	#Now determines the path to print the structure
		if inst.has_option(section,"output_file"):
			output_path = inst.get(section,"output_file")
		else:
			if struct.struct_id!=None:
				name = str(struct.struct_id)
			else:
				name = filename
			if inst.has_option(section,"output_suffix"):
				name += inst.get_eval(section,"output_suffix")[output_formats.index(output_format)]
			elif output_format == "json":
				name += ".json"
			elif output_format == "geometry":
				name += ".in"

			output_path = os.path.join(inst.get(section,"output_folder"), name )
		
		safe_make_dir(os.path.dirname(output_path))
		f = open(output_path,"w")
		if output_format == "geometry":
			f.write(struct.get_geometry_atom_format())
		elif output_format == "json":
			f.write(struct.dumps())
		f.close()
		inst.grant_permission(output_path)		

		if inst.get_info_level()>=2:
			write_log.write_master_log(inst,"1 structure printed to "+output_path+" in "+output_format+" format")

def dump_structure(struct, output_folder, output_formats, output_suffixs=None):

    for output_format in output_formats:
        if output_format != "json" and output_format != "geometry":
	    raise ValueError("dump_structure unknown output_format: "+output_format)


        name = struct.struct_id
        if output_suffixs!=None:
	    name += output_suffixs[output_formats.index(output_format)]
        elif output_format == "json":
            name += ".json"
        elif output_format == "geometry":
	    name += ".in"

	output_path = os.path.join(output_folder, name)
		
        safe_make_dir(os.path.dirname(output_path))
        f = open(output_path,"w")
        if output_format == "geometry":
	    f.write(struct.get_geometry_atom_format())
        elif output_format == "json":
	    f.write(struct.dumps())
	f.close()


def safe_copy_folder(src,dst):
	'''
	Safely copy the src folder to the destination directory
	i.e., if src = "/home/myfolder", and dst = "/home/shared/"
	will copy it to "/home/shared/myfolder"
	Will not copy if src is already in dst
	'''
	src_abs = os.path.abspath(src)
	dst_abs = os.path.abspath(dst)
	src_dir = os.path.dirname(src_abs)
	src_name = os.path.basename(src)
	if src_dir!=dst_abs:
		shutil.copytree(src_abs,os.path.join(dst_abs,src_name))

def safe_make_dir(directory):
	'''
	Silently make dir
	'''
	try:
		os.makedirs(directory)
	except:
		pass

def safe_rmdir(directory):
    try:
        os.rmdir(directory)
    except:
        pass

def safe_remove_files(file_list):
    for f in file_list:
        try:
            os.remove(f)
        except:
            pass

def list_directory_file(paths,suffix="",depth=0):
	'''
	Walks the directory until the given depth and give all the files and directories
	If path is a single element string, then only list one directory
	If path is a list, then list all of the folders in the directory
	'''
	if isinstance(paths,str):
		paths = [paths]

	for i in range (len(paths)):
		paths[i] = os.path.abspath(paths[i])

	results=[]
	for i in range (depth+1):
		if len(paths) == 0:
			break

		new_path = []
		for path in paths:
			results += [os.path.join(path,name) for name in os.listdir(path)\
			if os.path.isfile(os.path.join(path,name)) and (suffix=="" or name[len(name)-len(suffix):]==suffix)]
			new_path += [os.path.join(path,name) for name in os.listdir(path) if os.path.isdir(os.path.join(path,name))]
		
		paths = new_path[:]

	return results

def get_structure_list (inst,section,prioritize_path=True):
	'''
	Takes the instruction under the specified section and returns a list of structure paths
	If structure_path is present under the section, then returns a list of a single single structure path
	If structure_dir is used, then returns a list of paths
	'''
	if inst.has_option(section,"structure_path") and (prioritize_path or not inst.has_option(section,"structure_dir")):
		return [inst.get(section,"structure_path")]

	elif inst.has_option(section,"structure_dir"):
		inst.set_default(section,"structure_suffix","")
		inst.set_default(section,"structure_dir_depth","0")
		return list_directory_file(inst.get(section,"structure_dir"),\
		inst.get(section,"structure_suffix"),inst.get_eval(section,"structure_dir_depth"))
	
	else:
		write_log.write_master_err_and_raise(inst,"Section %s missing both necessary path or directory to structure files" % section)

def get_and_lock_file(paths,suffix="",depth=0,message="locked"):
	'''
	Randomly selects an available file under the paths
	Locks the file, and returns the file name
	'''
	list_of_files = list_directory_file(paths,suffix,depth)
	while len(list_of_files)>0:
		chosen = int(random.uniform(0,len(list_of_files)))
		name = list_of_files[chosen]
		if not os.path.isfile(name+".lock"):
			f = open(name+".lock","w")
			f.write(message)
			f.close()
			return name
		else:
			list_of_files.pop(chosen)
	return False

def retrieve_and_truncate_last_line(file_path):
	'''
	Returns the last non-empty line of the file
	Returns False if the file is empty
	'''
	f = open(file_path,"r+")
	f.seek(0,os.SEEK_END)
	pos = f.tell()-1
	result = ""
	while pos>=0:
		char = f.read(1)
		if char!="\n":
			result = char+result
		elif result!="":
			break
		pos -= 1
		if pos>=0:
			f.seek(pos,os.SEEK_SET)
	f.seek(pos+1,os.SEEK_SET)
	f.truncate()
	if result=="":
		return False
	else:
		return result

def retrieve_integer_and_subtract(file_path,subtract=1):
	'''
	Returns the integer stored in the file
	And subtract it by subtract
	'''
	dirname = os.path.dirname(file_path)
	filename = file_path[len(dirname)+1:]
	with FileLock(filename,dirname,600):
		f = open(file_path,"r")
		result = int(f.read())
		f.close()
		f = open(file_path,"w")
		f.write(str(result-subtract))
		f.close()
	return result

def unit_vector(vector):
	""" Returns the unit vector of the vector.  """
	return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
	'''
	Returns the angle between v1 and v2
	'''
	v1_u = unit_vector(v1)
	v2_u = unit_vector(v2)
	return np.rad2deg(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)))

def half_gaussian_sampling_lower(center,std,lowerbound=None):
	'''
	Gaussian sample that returns a value smaller than the center and greater than lowerbound
	'''
	if lowerbound!=None and lowerbound >= center:
		raise ValueError("half_gaussian_sampling_lower lowerbound greater or equal to center")

	if std == 0:
		return center

	k = np.random.normal(loc=center,scale=std)
	while k>=center or (lowerbound!=None and k<lowerbound):
		k = np.random.normal(loc=center,scale=std)

	return k

def half_gaussian_sampling_upper(center,std,upperbound=None):
	'''
	Gaussian sample that returns a value smaller than the center and greater than lowerbound
	'''
	if upperbound!=None and upperbound <= center:
		raise ValueError("half_gaussian_sampling_lower lowerbound greater or equal to center")

	if std == 0:
		return center

	k = np.random.normal(loc=center,scale=std)
	while k<=center or (upperbound!=None and k>upperbound):
		k = np.random.normal(loc=center,scale=std)

	return k


def write_active(fdir):
    '''
    Creates or overwrite a active.info file
    Recording current time stamp
    '''
    f=open(os.path.join(fdir,"active.info"),"w")
    f.write(str(time.time()))
    f.close()

def read_active(folder):
    '''
    Reads the active.info in the folder and returns the time
    '''
    try:
        f=open(os.path.join(folder,"active.info"))
        last_time=float(f.read())
        return last_time
    except:
        return False

	
