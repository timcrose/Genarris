"""
If any part of this module is used for a publication please cite:

X. Li, F. Curtis, T. Rose, C. Schober, A. Vazquez-Mayagoitia, K. Reuter,
H. Oberhofer, and N. Marom "Genarris: Random Generation of Molecular Crystal 
Structures and Fast Screening with a Harris Approximation, ",
J. Chem. Phys., DOI: 10.1063/1.5014038; arXiv 1803.02145 (2018)
"""
'''
Created on Jun 26, 2015
@author: Patrick Kilecdi
'''
import structure, structure_handling
from hashlib import sha1
import shutil
import os, time, multiprocessing, subprocess
#import bgqtools
import random
import numpy as np
import inspect
import utilities, copy
from external_libs.filelock import FileLock
from utilities import file_handling

__author__ = "Xiayue Li, Timothy Rose, Christoph Schober, and Farren Curtis"
__copyright__ = "Copyright 2018, Carnegie Mellon University and "+\
                "Fritz-Haber-Institut der Max-Planck-Gessellschaft"
__credits__ = ["Xiayue Li", "Luca Ghiringhelli", "Farren Curtis", "Tim Rose",
               "Christoph Schober", "Alvaro Vazquez-Mayagoita",
               "Karsten Reuter", "Harald Oberhofer", "Noa Marom"]
__license__ = "BSD-3"
__version__ = "180324"
__maintainer__ = "Timothy Rose"
__email__ = "trose@andrew.cmu.edu"
__url__ = "http://www.noamarom.com"



cwd=os.getcwd()
#print "This is working directory "+cwd
#bin="/projects/GAtor/avazquez/soft/aims.150205.scalapack.mpi.x"
bin="/lustre/project/nmarom/bin/aims.150205.scalapack.mpi.x"
tmp=os.path.join(cwd,"tmp/")
#print "This is temp directory "+tmp
relaxation="aims_light"
#folder_in_name="btm_4mpc_chiral_ga_2nd_1000"
#folder_in_name="tier_1_1-1117"
#folder_in_name="MBD_test_on_mira"
folder_in_name="tier_2_combined"
input_path=os.path.join(cwd,"inputs",folder_in_name)
other_path="/lustre/project/nmarom/pool_management_0.2/inputs/tier_2_combined_new_control"
#bk_path_success=os.path.join(cwd,"bk",folder_in_name); bk_enabled = False
#bk_path_failed=os.path.join(cwd,"bk",folder_in_name+"_failed") 
#bk_path_hung=os.path.join(cwd,"bk",folder_in_name+"_hung")
#bk_path_killed=os.path.join(cwd,"bk",folder_in_name+"_killed") 
#If a aims job fails though it has started outputing, the file will be backed up to the _failed path
#when the job gets hung, whether at launching or in the middle, it is back-uped to hung
#When remove_lock updates a structure from a folder that a job finished left behind, the bk folder will go to bk_path_killed
output_path=os.path.join(cwd,"outputs","tier_2_combined")
#To begin, place the pool in the input_path, and prepare the output_path directory

def main():
	'''
	Since the things to be done with done with this pool management can be so diverse,
	there wouldn't be a set routine and user interface. 
	Simple calls of the functions should help effectively maange a large pool
	'''
	coll=[]
	make_directory(tmp)
	make_directory(output_path)
#	if bk_enabled:
#		make_directory(bk_path)
	write_log("",time_print=False)
	write_log("",time_print=False)
	write_log("",time_print=False)
	write_log("new pool management task initiating at "+time.strftime("%Y-%m-%d %H:%M:%S"))
	#Routine command for structure pulling
	coll+=pull_structures(input_path,type="flat",suffix=".json",parse="loads",id_assign=None)
	for struct in coll:
		try:
			if struct.properties["energy_tier_2"]!=None:
				struct.properties["energy"]=float(struct.properties["energy_tier_2"])
		except:
			pass
	coll+=pull_structures(other_path,type="flat",suffix=".json",parse="loads",id_assign=None)

#	coll+=pull_structures(input_path,type="flat",suffix=".next_step",parse="build_geo_whole_atom_format",id_assign="all",id_assign_type="folder_name",update_enabled=False,energy_extract="aims.out")
	
	#Optional cleaning of the output path
#	clean_directory(output_path)

	#Routine 1: sorting a collection according to a keyword in properties
#	collection_sort(coll,keyword="energy",reverse=False)
#	collection_print(coll,output_path,type="file",naming_scheme="list_index",suffix=".data",parse="dumps") 
	#In order to preserve the ranking, you might want to name the structures using the list_index scheme

	
	#Routine 2: doing batch operation on a collection
	#Look in structure.handling for a list of currently featured operations
	#Let me know if you have a feature that you would like to add
#	coll+=pull_structures(input_path,type="file",suffix=".json",parse="json")
#	import structure_handling
#	collection_operation(coll,structure_handling.cell_lower_triangular,create_duplicate=False)

	#Routine 3: starting collection batch relax; first call
#	collection_relax_cypress(coll=coll,remove_lock=True)
#	collection_relax_bgq("cetus",coll=coll,remove_lock=False)

	#Routine 4: calling an duplicate of the relaxation scheme
#	collection_relax_cypress(coll=None,remove_lock=False)
#	collection_relax_bgq("cetus",coll=None,remove_lock=False)

	#Routine 5: restarting a relaxation scheme after failure or walltime runout
#	collection_relax_cypress(coll=None, remove_lock=True)
#	collection_relax_bgq("Cetus",nodes_per_replica=128,coll=None,remove_lock=True)

        #------------------Remark on relaxation------------------
        #To do batch relaxation on Cypress, first submit a job that uses routine 3
        #and then submit multiple jobs with routine 4
        #If the job is not finished in time, do Routine 5


	#Routine 6: printing certain features of the structure along side
#	import get_properties
#	collection_extract(coll,[get_properties.get_id,get_properties.get_energy,get_properties.get_volume],outf=folder_in_name+".info")

#        import molecular_descriptors
#        import get_properties
#        collection_extract_multiprocessing(coll,[get_properties.get_energy,[molecular_descriptors.coulombic_matrix,{"super_cell":[-1,1]}]],outf="cm_from_"+folder_in_name+".info",number_of_multi=20)
	
	#Routine 7: does diversity analysis on the pool
#	import molecular_descriptors
#	import get_properties
#	import structure_comparison
#	print collection_diversity_pairwise_random(coll,[get_properties.get_energy_difference,[molecular_descriptors.coulombic_matrix_comparison,{"super_cell":False}],[structure_comparison.is_duplicate,[2,10],{'create_duplicate':True,"extended_check":False}]],report="diversity_btm_sym_2_whole_random.out")
#        print collection_diversity_pairwise_random(coll[:1000],[get_properties.get_energy_difference,[molecular_descriptors.coulombic_matrix_comparison,{"super_cell":False}],[structure_comparison.is_duplicate,[2,15],{'create_duplicate':True,"extended_check":False}]],report="diversity_btm_sym_2_top_1000_random.out")
#        print collection_diversity_pairwise_random(coll[:100],[get_properties.get_energy_difference,[molecular_descriptors.coulombic_matrix_comparison,{"super_cell":False}],[structure_comparison.is_duplicate,[2,15],{'create_duplicate':True,"extended_check":False}]],report="diversity_btm_sym_2_top_100_random.out")
#	print collection_diversity_pairwise_full(coll,[get_properties.get_energy_difference,[structure_comparison.is_duplicate,[4,15],{'create_duplicate':True,"extended_check":False}]],report="diversity_"+folder_in_name+".out")
#	print collection_diversity_pairwise_full(coll,[[structure_comparison.is_duplicate,[4,15],{'create_duplicate':True,"extended_check":False}]],report="diversity_"+folder_in_name+".out")

	#Routine 8: manipulates propreties of the structures in the collection
#	property_set(coll,"mixed_ranking","list_index")
#	property_set(coll,"harris_ranking","struct_id")
#	property_set(coll,"energy_relaxed_geometry","property","energy")

#	clean_directory(output_path)
#	for i in range(len(coll)):
#		if coll[i].properties["energy_tier1"]<-100000:
#			coll[i].properties["energy_per_mol"]=coll[i].properties["energy_tier1"]/4
#		else:
#			coll[i].properties["energy_per_mol"]=coll[i].properties["energy_tier1"]/2
#	collection_sort(coll,keyword="energy_per_mol",reverse=False)
#			sum_4=i
#			break
#	print "This is how many 4mpc "+str(sum_4)
#	for i in range (len(coll)):
#		if coll[i].properties["energy"]>-100000:
#			break
#	print "this is i", i
#	coll=collection_remove_duplicate(coll[:i],4,15,energy_window=0.05,is_duplicate_tolerance=50,extended_check=False,report="tier_2_complete_duplicate_4mpc.info",report_style="pairs")

#	coll=collection_remove_duplicate(coll[1000:2000],4,15,energy_window=0.05,is_duplicate_tolerance=1,extended_check=False,report="4mpc_chiral_2nd_1000_duplicate.info")

	
	#Routine command for collection outputs
	#-------------------------Warning!-------------------------------
	#Warning: if only relaxation is called, no need to print the collection
	#Relaxation outputs new structure to the output_path as it goes
	collection_print(coll,output_path,type="file",naming_scheme="struct_id",suffix=".json",parse="dumps",zero_fill=4)
#	collection_print(coll[200:400],output_path+"2nd_200",type="file",naming_scheme="struct_id",suffix=".json",parse="dumps",zero_fill=4)
#	collection_print(coll[400:600],output_path+"3rd_200",type="file",naming_scheme="struct_id",suffix=".json",parse="dumps",zero_fill=4)
#	collection_print(coll[600:],output_path+"4th_200",type="file",naming_scheme="struct_id",suffix=".json",parse="dumps",zero_fill=4)
#	collection_print(coll,output_path,type="file",naming_scheme="struct_id",suffix=".in",parse="get_geometry_atom_format",zero_fill=4)
#	collection_print(coll],output_path,type="file",naming_scheme="list_index",suffix=".json",parse="dumps")
	
def pull_structures(input_path,type="flat",suffix=None,parse="loads",log=True,id_assign=None, id_assign_type="file_name", update_enabled=False, energy_extract=None, geometry_extract=None):
	'''
	type="flat" reads in all files in the current input path and that in the folders
	Essentially flattens the input_path and conduct Type="file"
	type="file" reads in only the files in the input_path
	type="folder" reads in only the files in the folders

	suffix=None allows all types of input suffix
	suffix=".json" would require the file name to have the suffix of .json

	parse: the attribute name of the Structure class that would build a sruct from a string

	id_assign will enable assigning struct.id
	id_assign="all" assigns all id
	id_assign="if_absent" assigns only when the struct_id is not already assigned

	id_assign_type="file_name": assigns the file from which the structure is pulled as the struct.id
	id_assign_type="folder_name": assigns the folder name as the struct.id

	update_enabled=True will disable pulling in of structures with .more suffix, and yet try to update any none .more structure with a .more having the same name
	The .more file must be a json readable string

	energy_extract="aims.out": Only useful currently for the folder, attempts to find the file aims.out, and extract the energy and update the structure energy
	geometry_extract="geometry.in.next_step": Only useful currently for the folder mode, attempts to find geometry.in.next_step and update the geometry of the structure; if geometry_extract is set as any other string than that, it will be understood to be an aims output

	Returns a list of structures
	'''
	coll=[]
	if id_assign!=None and id_assign!="all" and id_assign!="if_absent":
		raise ValueError(id_assign+": unknown type of id_assign for pull_structures")
	if type=="flat" or type=="file":
		list_of_files = [name for name in os.listdir(input_path) \
				if os.path.isfile(os.path.join(input_path,name))]
		for name in list_of_files:
			if suffix==None or name[len(name)-len(suffix):]==suffix\
			and not (update_enabled and name[len(name)-5:]==".more")\
			and not name==energy_extract and not name==geometry_extract:
				f=open(os.path.join(input_path,name))
				input_string=f.read()
#				print input_string
				new_structure=structure.Structure()
				method = getattr(new_structure,parse)
				method(input_string)

				if id_assign=="all" or (id_assign=="if_absent" and new_structure.struct_id==None):
					if id_assign_type=="file_name":
						new_structure.struct_id=name[:len(name)-len(suffix)]
					elif id_assign_type=="folder_name":
						new_structure.struct_id=os.path.basename(input_path)
					else:
						raise ValueError("Unknown id_assign_typ")
				if update_enabled and os.path.isfile(os.path.join(input_path,name+".more")):
					f=open(os.path.join(input_path,name+".more"))
					input_string=f.read()
					old_structure=structure.Structure()
					old_structure.loads(input_string)
					structure_handling.struct_properties_updates(new_structure,old_structure)
				if energy_extract!=None and os.path.isfile(os.path.join(input_path,energy_extract)): #The structures being read in will have their energy updated
					new_structure.properties["energy"]=file_handling.extract_energy_aims(os.path.join(input_path,energy_extract))


				coll.append(new_structure)
				print "New structure collected! length of coll %i" % len(coll)

	if type=="flat" or type=="folder":
		list_of_folders = [name for name in os.listdir(input_path) \
				if os.path.isdir(os.path.join(input_path,name))]
		for name in list_of_folders:
			coll+=pull_structures(os.path.join(input_path,name),type="flat",suffix=suffix,parse=parse,log=False,id_assign=id_assign,id_assign_type=id_assign_type,update_enabled=update_enabled,energy_extract=energy_extract,geometry_extract=geometry_extract)

	if log: #This is the outermost layer
		write_log("pull_structures has been called to retrieve structure from %s \n%i structures have been retrieved" % (input_path,len(coll)))	
	return coll

def collection_diversity_pairwise_full (coll,list_of_functions,report=""):
	'''
	Conducts a diversity check on the collection
	coll: the collection being checked
	list_of_functions: list of functions that has the first two arguments as two structures
	with arguments attached to the function
	and outputs a single scalar value
	
	report: directory of a pair-wise report on the diversity
	'''
	write_log("collection_diversity_pairwise_full is called")
	write_log("collection length is %i" % len(coll))
	for i in range (len(list_of_functions)):
		try:
			list_of_functions[i]=list(list_of_functions[i])
		except:
			list_of_functions[i]=[list_of_functions[i]]

		if not inspect.isfunction(list_of_functions[i][0]):
			raise ValueError("Non-function received")
		if len(list_of_functions[i])==1:
			list_of_functions[i]=[list_of_functions[i][0],[],{}]
		elif len(list_of_functions[i])==2:
			if isinstance(list_of_functions[i][1],dict):
				list_of_functions[i]=[list_of_functions[i][0],[],list_of_functions[i][1]]
			else:
				list_of_functions[i].append({})
	
	diff=[[[] for i in range (len(coll))] for j in range (len(coll))]
	for i in range (len(coll)):
		for j in range(i+1,len(coll)):
			for info in list_of_functions:
				diff[i][j].append(info[0](coll[i],coll[j],*info[1],**info[2]))
#			print diff[i][j]
			if report!="":
				f=open(report,"a")
#				outs="%i %i " % (i,j)
				outs=coll[i].struct_id+" "+coll[j].struct_id+" "
				for value in diff[i][j]:
					outs+=str(value)+" "
				outs+="\n"
				f.write(outs)
				f.close()
				
	result=[0 for i in range (len(list_of_functions))]
	outs=""
	for i in range (len(coll)-1):
		for j in range (i+1,len(coll)):
			outs+=str(i)+" "+str(j)+" "
			for k in range (len(list_of_functions)):
				result[k]+=diff[i][j][k]
				outs+=str(diff[i][j][k])+" "
			outs+="\n"
	for k in range (len(list_of_functions)):
		result[k]/=len(coll)**2
	write_log("This is the result of diversity check "+str(result))
	return result
	
def collection_diversity_pairwise_random(coll,list_of_functions,sampling=2000,report=""):
	'''
	Conducts a diversity check on the collection
	coll: the collection being checked
	list_of_functions: list of functions that has the first two arguments as two structures
	with arguments attached to the function
	and outputs a single scalar value
	
	sampling: number of random pairs selected from the collection
	
	report: directory of a pair-wise report on the diversity
	'''
	write_log("collection_diversity_pairwise_random is called")
	write_log("collection length is %i, sampling is %i" % (len(coll),sampling))
	for i in range (len(list_of_functions)):
		try:
			list_of_functions[i]=list(list_of_functions[i])
		except:
			list_of_functions[i]=[list_of_functions[i]]

		if not inspect.isfunction(list_of_functions[i][0]):
			raise ValueError("Non-function received")
		if len(list_of_functions[i])==1:
			list_of_functions[i]=[list_of_functions[i][0],[],{}]
		elif len(list_of_functions[i])==2:
			if isinstance(list_of_functions[i][1],dict):
				list_of_functions[i]=[list_of_functions[i][0],[],list_of_functions[i][1]]
			else:
				list_of_functions[i].append({})
	
	result=[0 for i in range (len(list_of_functions))]
	for i in range (sampling):
		s1=int(random.uniform(0,len(coll)))
		s2=int(random.uniform(0,len(coll)))
		outs="%i %i %i " % (i,s1,s2)
		for info in list_of_functions:
			diff = info[0](coll[s1],coll[s2],*info[1],**info[2])
			result[list_of_functions.index(info)]+=diff
			outs+=str(diff)+" "
		outs+="\n"
		if report!="":
			f=open(report,"a")
			f.write(outs)
			f.close()
			
	for k in range (len(list_of_functions)):
		result[k]/=sampling
	write_log("This is the result of diversity check "+str(result))
	return result


def collection_extract(coll,list_of_functions,outf=""):
	'''
	Produce a list of lists of properties extracted from the structures 
	
	coll: collection to be extracted

	list_of_functions: functions that take struct as an argument and returns the properties
	see get_properties.py for a list of featured functions

	outf: optional; the file name where the information will be appended to the file
	'''
	result=[]
	for struct in coll:
		result.append([function(struct) for function in list_of_functions])
	if outf!='':
		f=open(outf,"a")
		out_string=''
		for llist in result:
			new_line=''
			for item in llist:
				new_line+=str(item)+' '
			out_string+=new_line+"\n"
		f.write(out_string)
		f.close()
	write_log("Collection_extract is called. Results outputed to "+outf)
	return result

def collection_extract_multiprocessing(coll,list_of_functions,outf="",number_of_multi=1):
        '''
        Takes in a collection and list_of_functions, and returns a list of the result
        if multi>1, multiprocessing will be used to call struct_extract
        '''
        def structure_extract(struct):
                return [info[0](struct,*info[1],**info[2]) for info in list_of_functions]

        if number_of_multi<1:
                raise ValueError("multi<1")

	write_log("collection_extract_multiprocessing is called. number_of_multi=%i" % (number_of_multi))
        list_of_functions=utilities.normalize_list_of_functions(list_of_functions)

        if number_of_multi==1:
                result=[]
                for struct in coll:
                        result.append(structure_extract(struct))
        else:
                pool=multiprocessing.Pool(processes=number_of_multi)
                result=pool.map(structure_extract,coll)

        if outf!="":
                f=open(outf,"w")
		for info in result:
	                f.write(str(info)+"\n")
                f.close()
        return result
	

def collection_operation(coll,function,args=[],create_duplicate=True):
	'''
	Carries operation on an entire collection

	coll: the collection that the operation will be done on

	function: the operation that will be carried out on the collection
	should take struct() as the first argument, and return a struct()

	args: optional; list or tuple of constant of functions
	Warning: make sure the functions in args take one structure as an argument
	And returns the desired argument for the main function of this operation

	create_duplicate: whether or not the coll will be overwritten
	'''
	if create_duplicate:
		coll=copy.deepcopy(coll)
	if isinstance(args,tuple):
		args=list[args]
	for i in range (len(coll)):

		current_args=[]
		for arg in args:
			if inspect.isfunction(arg):
				current_args.append(arg(coll[i]))
			else:
				current_args.append(arg)

		coll[i]=function(coll[i],*current_args)
	
	return coll

def collection_print(coll,output,type="file",naming_scheme="random",suffix="",parse="dumps",zero_fill=0):
	'''
	Outputs the structures in coll to a given path
	type="file": outputs indidual file with the structure information
	type="folder": creates folders and place individual file under different folder
	naming_scheme="random": gets a random index from get_random_index
	naming_scheme="index": name the structures according to its index in the list
	naming_scheme="struct_id": name the structures according to its struct_id;
	if none, get_random_index will be called to assign one
	suffix=".json": adds ".json" to the file name
	'''
	make_directory(output)
	for struct in coll:
		method=getattr(struct,parse)
		out_str=method()
		if naming_scheme=="random":
			name=get_random_index()
		elif naming_scheme=="list_index":
			name=str(coll.index(struct))
		elif naming_scheme=="struct_id":
			try:
				name=struct.struct_id
				if name=='' or name==None:
					name=get_random_index()
			except:
				name=get_random_index()
		else:
			raise NameError("Unknown type of naming_scheme")
		if type=="file":
			f=open(os.path.join(output,name.zfill(zero_fill)+suffix),"w")
			f.write(out_str)
			f.close()
		elif type=="folder":
			make_directory(os.path.join(output,name))
			f=open(os.path.join(output,name.zfill(zero_fill),name.zfill(zero_fill)+suffix),"w")
			f.write(out_str)
			f.close()
		else:
			raise NameError("Unknown output type")		
#	write_log("%i structures have been printed to directory: %s" % (len(coll),output))

def collection_relax_bgq(machine,nodes_per_replica=32,coll=None,remove_lock=False):
	'''
	Relaxes the collection on Mira or Cetus
	machine="cetus" or "Cetus" will boot blocks of 128 to the smallest
	machine="Mira" will boot blocks of 512 
	if coll!=None, will print the collection to the temporary file
	Boot the blocks and multiprocess single_relax
	If remove_lock is enabled, will remove the lock files in the tmp directory
	'''
	if not os.path.isfile(os.path.join(cwd,"control.in")):
		raise RuntimeError("control.in not found in the working directory")
		return
	if coll!=None:
		collection_print(coll,tmp,naming_scheme="struct_id",parse="get_geometry_atom_format")
		collection_print(coll,tmp,naming_scheme="struct_id",suffix=".more",parse="dumps")
		write_log("%i structures added to tmp for relaxation" % len(coll))
	if machine=="cetus":
		machine=="Cetus"
	if machine=="mira":
		machine=="Mira"
	
	partsize, partition, job_id = bgqtools.get_cobalt_info()
	if (partsize<128 and machine=="Cetus") or (partsize<512 and machine=="Mira"):
		raise RuntimeError("Please submit a job with more nodes")
		return
	
	if machine=="Cetus":
		block_size=min(partsize,128)
	else:
		block_size=min(partsize,512)

	blocks = bgqtools.get_bootable_blocks(partition,partsize,block_size)
	bgqtools.boot_blocks(blocks)
	corners = bgqtools.block_corner_iter(blocks,nodes_per_replica)

	if remove_lock:
		scavenge_tmp()


	replica_name=[] #Sets up the replica names with the corner information encoded
	while True:
		try:
			corner=next(corners)
			replica_name.append([corner[0]+'%'+corner[1]+'%'+corner[2],nodes_per_replica,machine])
		except:
			break

	number_of_multi=len(replica_name)
	pool=multiprocessing.Pool(processes=number_of_multi)
	pool.map(relax_serial,replica_name)
#	relax_serial(replica_name[0])
	
def collection_relax_cypress(nodes_per_replica=2,coll=None,remove_lock=False):
        '''
        Relaxes the collection on Cypress
        Only runs one process; requires external multiple submission to achieve true parallelism
	
	Setting remove_lock=True not only removes the .lock files, but also tries to take out the geometry.in.next_step file 
        '''
        if not os.path.isfile(os.path.join(cwd,"control.in")):
                raise RuntimeError("control.in not found in the working directory")
                return
        if coll!=None:
                collection_print(coll,tmp,naming_scheme="struct_id",parse="get_geometry_atom_format")
		collection_print(coll,tmp,naming_scheme="struct_id",suffix=".more",parse="dumps")		
        if remove_lock:
		scavenge_tmp()

        replica_info=[get_random_index(),nodes_per_replica,"Cypress"]
        relax_serial(replica_info)

def collection_remove_duplicate(coll, nmpc, napm, energy_window=0.1, is_duplicate_tolerance=1, extended_check=False, keep_rule="energy_low", create_duplicate=True,report="",report_style="duplicates"):
	'''
	Removes the duplicate within the collection by using structure_comparison.is_duplicate()
	energy_window: only structures with energy difference smaller than that will be compared; set to a large value to compare all
	is_duplicate_tolerance: if is_duplicate_residual<tolerance, then determined as duplicate
	extended_check: whether or not the is_duplicate calling will involve setting up super_cell
	keep_rule: rule to decide which one in the duplicate pair to keep
	create_duplicate: whether to do a deepcopy on coll at the beginning
	report: the file to report this duplicate check
	report_stype="duplicates": only prints out the duplicates to be removed
	report_stype="pairs": prints out the pairs of duplicates
	'''
	def write_report(ind):
		if report!="":
			f=open(report,"a")
			if report_style=="duplicate":
				f.write(coll[ind[1]].struct_id+" ")
			elif report_style=="pairs":
				f.write(coll[ind[0]].struct_id+" "+coll[ind[1]].struct_id+"\n")
			f.close()

	import structure_comparison
	if create_duplicate:
		coll = copy.deepcopy(coll)
	write_log("collection_remove_duplicate is called to work on a collection of length %i, nmpc=%i, napm=%i, energy_window=%f, is_duplicate_tolerance=%f, extended_check=%s, keep_rule=%s" % (len(coll),nmpc,napm,energy_window,is_duplicate_tolerance,str(extended_check),keep_rule))
	left = 0
	while left < len(coll)-1:
		right = left+1
		remove = False #remove = True will mean that the current left should be thrown out
		while right<len(coll)-1:
			print coll[right].struct_id, coll[right].properties["energy"]
			print coll[left].struct_id, coll[left].properties["energy"]
			if abs(coll[right].properties["energy"]-coll[left].properties["energy"])>energy_window:
				right+=1
				continue
			resi = structure_comparison.is_duplicate(coll[left],coll[right],nmpc,napm,create_duplicate=True,extended_check=extended_check)
			if resi<is_duplicate_tolerance:
				if keep_rule=="energy_low":
					if coll[left].properties["energy"]<coll[right].properties["energy"]:
						write_report([left,right])
						coll.pop(right)
						right -= 1 #Move it back one so that when moved forward 1 later, remains at the same spot
					else:
						remove = True
						break
				elif keep_rule=="index_low":
					write_report([left,right])
					coll.pop(right)
					right -= 1
				else:
					raise ValueError("Unknown type of keep_rule")
			right += 1
		if remove:
			write_report([right,left]) #Reversed because the left one is the one to be removed
			coll.pop(left)
		else:
			left += 1
	write_log("Final length of collection: "+str(len(coll)))
	return coll

def collection_sort(coll,keyword="energy",reverse=False):
	'''
	Quicksorts the structure collection by the key word passed in
	'''
	coll.sort(key=lambda struct: getattr(struct,"properties")[keyword],reverse=reverse)

def collection_struct_id_assign(coll,struct_id_scheme=["random_index"]):
	'''
	Generates struct_id for the batch of structure using the struct_id_scheme
	'''
	for struct in coll:
		name = ""
		for item in struct_id_scheme:
			if item == "list_index":
				name += str(coll.index(struct)).zfill(len(str(len(coll))))
			elif item == "random_index":
				name += get_random_index()
			elif item in struct.properties:
				name += str(struct.properties[item])
			else:
				name += item
		struct.struct_id = name

def relax_serial(replica_info):
	'''
	Runs and operate the Relax() class
	'''
	replica=replica_info[0]
	nodes_per_replica=replica_info[1]
	machine=replica_info[2]

        if not os.path.isfile(os.path.join(cwd,"control.in")):
                raise RuntimeError("control.in not found in the working directory")
                return
	if not os.path.isdir(os.path.join(cwd,"tmp")):
		raise RuntimeError("tmp file not created! It should hold the structures to be relaxed")
		return
	
	time.sleep(random.uniform(0,5)) #Avoids too many replicas beginning at the same time

	process=Relax(replica,machine,tmp,cwd,bin,nodes_per_replica)
	count=0
	while True:

		found=process.setup()

		if found==False: #No more structures found in the tmp file
			write_log("Process %s cannot find anymore structure to relax in tmp, exiting" % (replica))
			remove_directory(os.path.join(tmp,replica))
			break

		write_log("Process %s beginning fhi-aims relaxation on structure %s" % (replica,process.filename))
		process.execute()

		if process.is_successful():

			write_log("Process %s now extracting geometry from structure %s" % (replica,process.filename))
			relaxed_struct=process.extract()
			relaxed_struct.struct_id=process.filename

#			try:
			f=open(os.path.join(tmp,process.filename+".more"),"r")
			write_log("Process %s found a .more file for structure %s" %(replica,process.filename))
			old_struct=structure.Structure()
			old_struct.loads(f.read())
			structure_handling.struct_properties_updates(relaxed_struct,old_struct)
			f.close()
			os.remove(os.path.join(tmp,process.filename+".more"))
#			except:
#				write_log("Process %s fails to update the new structure %s with original information" % (replica,process.filename))

			collection_print([relaxed_struct],output_path,type="file",naming_scheme="struct_id",suffix=".json",parse="dumps")

			if bk_enabled:
				bk_folder(tmp,replica,bk_path_success,"const",process.filename)

			process.cleanup()
			write_log("Process %s successfully relaxed structure %s, outputed to directory %s" % (replica,process.filename,output_path))
			count=0

		else:
			write_log("Relaxation failure for process %s with structure %s" % (replica,process.filename))
			#Currently, relaxation failure may be caused by job termination. 
			#Will back up the folder again to tmp file, to be picked up by the next scavenger scheme
			#And the .lock file is not removed until then
			if bk_enabled:
				bk_folder(tmp,replica,bk_path_failed,"random")
				bk_folder(tmp,replica,tmp,"random") 
			process.filename=None
			count+=1
			if count>10:
#				raise RuntimeError("Repeated relaxation failure; check relaxation module")
				write_log("Repeated relaxation failure for replica "+replica)
				break

def property_set(coll,name,type="const",value=None):
	'''
	Sets the struct.propreties[name} according to the type
	'''
	write_log("property_set is called to set the %s property, according to type=%s" % (name,type))
	if type=="const":
		write_log("The constant property set is "+str(value))
	for struct in coll:
		if type=="const":
			struct.properties[name]=value
		elif type=="list_index":
			struct.properties[name]=coll.index(struct)
		elif type=="struct_id":
			struct.properties[name]=struct.struct_id
		elif type=="property":
			struct.properties[name]=struct.properties[value]
		else:
			raise ValueError("Unknown type of assignment in property_set")
	return coll

def scavenge_tmp():
	'''
	Called if remove_lock is set to true
	Searches for valid relaxation folders in tmp (ones with restart.info)
	Updates the geometry outside the folder with geometry.in.next_step
	Copies the folder to bk with a random name if bk_enabled set to True
	Removes the folder
	'''
	time.sleep(30)
	list_of_folders = [name for name in os.listdir(tmp) \
			if os.path.isdir(os.path.join(tmp,name))]	
	for folder in list_of_folders:
		if os.path.exists(os.path.join(tmp,folder,"restart.info")):
			write_log("scavenge_tmp found folder %s containing restart.info" % folder)
			last_active_time=read_active(os.path.join(tmp,folder))
			if last_active_time==False:
				write_log("Unable to find active.info in the folder. Moving on")
				continue
			current_time=time.time()
			if current_time-last_active_time<30:
				write_log("Last update of active.info within 30 seconds. Moving on")
				continue

			f=open(os.path.join(tmp,folder,"restart.info"),"r")
			filename=f.read()
			f.close()

			if os.path.exists(os.path.join(tmp,folder,"geometry.in.next_step")):
				restart_structure=structure.Structure()
				try:
					restart_structure.build_geo_from_atom_file(os.path.join(tmp,folder,"geometry.in.next_step"))
				except: #Bad restart
					try:
						write_log("geometry.in.next_step potentially corrupted")
						shutil.rmtree(os.path.join(tmp,folder))
					except:
						write_log("Unable to remove the old folder ; only going to remove active info")
						if os.path.exists(os.path.join(tmp,folder,"active.info")):
							os.remove(os.path.join(tmp,folder,"active.info"))
					continue
					
				old_structure=structure.Structure()
				old_structure.build_geo_from_atom_file(os.path.join(tmp,filename))
				if len(old_structure.geometry)==len(restart_structure.geometry):
					write_log("geometry.in.next_step found in the folder %s ; updating to tmp file" % (folder))
					f=open(os.path.join(tmp,filename),"w")
					f.write(restart_structure.get_geometry_atom_format())
					f.close()
					os.system("chmod g=u "+os.path.join(tmp,filename))
				else:
					write_log("geometry.in.next_step found in the folder %s , but compromised; geometry not updated")
			if bk_enabled:
				bk_folder(tmp,folder,bk_path_killed,"random")
			try:
				shutil.rmtree(os.path.join(tmp,folder))
			except:
				write_log("Unable to remove the scavenged folder ; only removing active.info")
				if os.path.exists(os.path.join(tmp,folder,"active.info")):
					os.remove(os.path.join(tmp,folder,"active.info"))

			try:
				os.remove(os.path.join(tmp,filename+".lock"))
			except:
				write_log("WARNING! Unable to remove "+filename+".lock")
				pass

			write_log("Folder %s cleaned and .lock file removed" % folder)

def scavenge_folder(folder,check_active=True,active_wait=30):
	'''
	Scavenges a particular folder in tmp folder
	
	if check_active set to true, will compare the time in active.info with the current time
	if the time difference is greater than active_wait, then the folder will be scavenged
	returns True if the function calling this function should move to consider removing the folder
	
	WARNING! This function is not group-running friendly (account permission not reset)
	'''
	write_log("scavenge_folder scavenging folder %s." % folder)
	if not os.path.exists(os.path.join(tmp,folder,"restart.info")):
		write_log("Folder %s does not have a restart.info file. Moving on." % folder)
	if check_active:
		last_active_time=read_active(os.path.join(tmp,folder))
		if last_active_time==False:
			write_log("check_active set to True, but unable to find active.info in the folder. Moving on")
			return False
		current_time=time.time()
		if current_time-last_active_time<active_wait:
			write_log("Last update of active.info within %i seconds. Moving on" % (active_wait))
			return False

	f=open(os.path.join(tmp,folder,"restart.info"),"r")
	filename=f.read()
	if not os.path.isfile(os.path.join(tmp,filename)):
		write_log("File indicated by restart.info in folder %s does not exist in the tmp directory. Moving on." % folder)
		return False
		
	f.close()

	if os.path.exists(os.path.join(tmp,folder,"geometry.in.next_step")):
		restart_structure=structure.Structure()
		try:
			restart_structure.build_geo_from_atom_file(os.path.join(tmp,folder,"geometry.in.next_step"))
		except: #Bad restart
			write_log("geometry.in.next_step file in folder %s corrupted. Moving on")
			try:
				os.remove(os.path.join(tmp,filename+".lock"))
			except:
				pass
			return True
	
		old_structure=structure.Structure()
		old_structure.build_geo_from_atom_file(os.path.join(tmp,filename))
		if len(old_structure.geometry)==len(restart_structure.geometry):
			write_log("geometry.in.next_step found in the folder %s ; updating to tmp file" % (folder))
			f=open(os.path.join(tmp,filename),"w")
			f.write(restart_structure.get_geometry_atom_format())
			f.close()
		else:
			write_log("geometry.in.next_step found in the folder %s , but compromised; geometry not updated")


class Relax():
	'''
	Conducts relaxation with a given replica name
	Read from the tmp directory for structures to be relaxed
	Make sure a control.in file is located in the working directory
	'''
	def __init__(self,replica,mode,tmp,cwd,bin,nodes):
		'''
		replica is the replica name. The aims job will run under tmp/replica
		mode tells Relax() what type of job submission technique should be used
		tmp provides the tmp directory
		cwd provides the directory where the control.in can be found
		bin provides the path to aims binary
		'''
 
		self.replica=replica
		self.mode=mode
		self.filename=""
		self.tmp=tmp
		self.cwd=cwd
		self.bin=bin
		self.working_dir=os.path.join(tmp,replica)
		self.nodes=nodes

		if not os.path.isfile(os.path.join(self.cwd,"control.in")):
			raise RuntimeError("control.in not found in the working directory")
			return
		if not os.path.isdir(tmp):
			raise RuntimeError("tmp directory not found")
			return

	def setup(self):
		clean_directory(os.path.join(self.tmp,self.replica+'/'))
		success=self.get_structure()
		if success==False: #No more structures in the tmp directory
			return False
		shutil.copy(os.path.join(self.cwd,"control.in"),os.path.join(self.tmp,self.replica))
		shutil.copyfile(os.path.join(self.tmp,self.filename),os.path.join(self.tmp,self.replica,"geometry.in"))
		f=open(os.path.join(self.tmp,self.replica,"restart.info"),"w")
		f.write(self.filename)
		f.close()
		write_active(self.working_dir)
		self.set_permission() #Set group permission same as user

	def execute(self):
		'''
		Calls the aims relaxation job according the mode being specified
		'''
		if self.mode=="Cetus" or self.mode=="cetus" or self.mode=="Mira" or self.mode=="mira" or self.mode=="bgq":
			block_size=self.nodes
			#Will run it with modes=4 and thre=4
			modes=4; thres=4
			try:
				l=self.replica.index("%")
				block=self.replica[0:l]
				rest=self.replica[l+1:]
				l=rest.index("%")
				corner=rest[:l]
				shape=rest[l+1:]
				arglist=["runjob","--np",str(modes*block_size),"-p",str(modes),"--envs","OMP_NUM_THREADS="+str(thres),"--verbose","INFO","--block",block,"--corner",corner,"--shape",shape,"--cwd",self.working_dir,"--exe",self.bin]
#				command='runjob --np %i -p %i --envs OMP_NUM_THREADS=%i  --verbose=INFO --block %s --corner %s --shape %s --cwd %s : %s > %s &' % (modes*block_size,modes,thres,block,corner,shape,self.working_dir,self.bin,os.path.join(self.working_dir,'aims.out'))
			except: #Only has a block name
				arglist=["runjob","--np",str(modes*block_size),"-p",str(modes),"--envs","OMP_NUM_THREADS="+str(thres),"--verbose","INFO","--block",block,"--corner",corner,"--shape",shape,"--cwd",self.working_dir,"--exe",self.bin]
#				command='runjob --np %i -p %i --envs OMP_NUM_THREADS=%i  --verbose=INFO --block %s --cwd %s : %s > %s &' % (modes*block_size,modes,thres,self.replica,self.working_dir,self.bin,os.path.join(self.working_dir,'aims.out'))
#			os.system(command)
		elif self.mode=="Cypress" or self.mode=="cypress":
			arglist=["mpirun","-wdir",self.working_dir,self.bin]

		aimsout=os.path.join(self.working_dir,"aims.out")
		for i in range (10):
			outfile=open(aimsout,"w")
			get_execute_clearance(request_folder=self.working_dir)
#			print "this is arglist", arglist
			p=subprocess.Popen(arglist,stdout=outfile)
			time.sleep(1)
			try:
				status=p.poll()
			except: #OSError Errno3
				write_log("Replica %s encountering node failure ; unable to spawn process anymore" % self.replica)
				time.sleep(86400)
				return False
			time_limit=60
			for j in range (time_limit): #Allow 60 seconds for aims to start outputting
				if (p.poll()!=None) or (os.stat(aimsout).st_size>512):
					break
				write_active(self.working_dir)
				self.set_permission()
				time.sleep(1)
			if (os.stat(aimsout).st_size>512):
				write_log("Replica %s successfully launched aims job. Beginning outputing" % self.replica)
				break
			outfile.close()
			write_log("Replica %s failed to launch aims job." % self.replica)
			try:
				p.send_signal(2)
			except:
				write_log("Replica %s process no longer exists ; kill failure ; possible nodes failure" % self.replica)
				time.sleep(86400)
			active_sleep(60,self.working_dir)

			if i==9:
				write_log("Replica %s repeatedly fails to launch aims job. Exiting." % self.replica)
				return False
#				raise RuntimeError("Repeated failure to launch aims job in execute()")

		counter=0; last=os.stat(aimsout).st_size
		while counter<60 and p.poll()==None: #The output file needs to update at least once in every 5 minutes
			write_active(self.working_dir)
			self.set_permission()
			time.sleep(10) #This value should not be greater than 20 for active.info update
			if os.stat(aimsout).st_size>last:
				last=os.stat(aimsout).st_size
				counter=0
			else:
				counter+=1
		if counter==60:
			try:
				write_log("Replica %s aims job hung for some reason." % self.replica)
				p.send_signal(2)
				active_sleep(60,self.working_dir)
			except:
				write_log("Replica %s process unable to kill ; possible nodes failure ; now sleeping" % self.replica)
				time.sleep(86400)
				pass
		outfile.close()
		write_log("Replica %s aims job exited with status %s" % (self.replica,str(p.poll())))
		


	def extract(self):
	        '''
        	Reads the output files from relaxation and specifies properties of the structure (e.g. energy)
	        Returns: Structure or False
        	'''
	        # creates empty structure
        	self.result_struct = structure.Structure()
	        #try: self.result_struct.set_lattice_vectors(self.input_structure.get_lattice_vectors())
	
        	# reads output file to extract energy
	        energy = self.extract_energy()
        	print "energy", energy
	        if energy == False: return False  # did not converge
        	self.result_struct.set_property('energy', energy)

	        # sets the result structure's lattice vectors based on output file
        	lats= self.extract_lats()
	
	        latA = [float(lats[0][0]), float(lats[0][1]), float(lats[0][2])]
        	latB = [float(lats[1][0]), float(lats[1][1]), float(lats[1][2])]
	        latC = [float(lats[2][0]), float(lats[2][1]), float(lats[2][2])]
        	self.result_struct.set_property('lattice_vector_a', latA)
	        self.result_struct.set_property('lattice_vector_b', latB)
        	self.result_struct.set_property('lattice_vector_c', latC)

	        latA = np.asarray(latA)
        	latB = np.asarray(latB)
	        latC = np.asarray(latC)

        	temp_vol = np.dot(np.cross(latA, latB), latC)
	        alpha = self.angle(latB, latC)
        	beta = self.angle(latA, latC)
	        gamma = self.angle(latA, latB)
        	a = np.linalg.norm(latA)
	        b = np.linalg.norm(latB)
        	c = np.linalg.norm(latC)
	        self.result_struct.set_property('cell_vol', temp_vol)
        	self.result_struct.set_property('alpha',alpha)
	        self.result_struct.set_property('beta', beta)
        	self.result_struct.set_property('gamma', gamma)
	        self.result_struct.set_property('a',a)
        	self.result_struct.set_property('b', b)
	        self.result_struct.set_property('c', c)


        	print "latA",latA
	        print "latB", latB
        	print "latc", latC
	        print "alpha", alpha
	        print "beta", beta
	        print "gamma", gamma

        	# sets the result structure's geometry based on output file
	        extract_success = self.extract_geometry()

        	return self.result_struct

	def cleanup(self):
		'''
		Remove the structure file and .lock file
		'''
		try:
			os.remove(os.path.join(self.tmp,self.filename))
		except:
			pass
		try:
			os.remove(os.path.join(self.tmp,self.filename+'.lock'))
		except:
			pass
	
	def release(self):
		'''
		Only remove .lock file
		'''
		try:
			os.remove(os.path.join(self.tmp,self.filename+'.lock'))
		except:
			pass
		self.filename=None

	def extract_energy(self):
        	'''
	        reads the resulting energy from the FHI-aims log
        	Returns: float if successful, False if unsuccessful
	        '''
        	aims_out = open(os.path.join(self.working_dir, 'aims.out'))
	        while True:
        	    line = aims_out.readline()
		    if not line: return False  # energy not converged
	            if '  | Total energy corrected        :' in line:
	                tokens = line.split()
        	        energy = float(tokens[5])  # converts from SI string to float
                	return energy

	def extract_lats(self):
        	print "extract lats..."
	        aims_out = open(os.path.join(self.working_dir, 'aims.out'))
        	while True:
                	line = aims_out.readline()
	                if 'Final atomic structure' in line: break
        	aims_out.readline()  # skip one line
	        atom_string = ''
        	while True:
                	line = aims_out.readline()
	                if 'atom' in line: break
        	        else: atom_string += line

	        latout = atom_string.split()
        	lat_A = [latout[1], latout[2], latout[3]]
	        lat_B = [latout[5], latout[6], latout[7]]
        	lat_C = [latout[9], latout[10], latout[11]]
	        lats = [lat_A, lat_B, lat_C]

        	return lats

	def angle(self,v1,v2):
	        numdot = np.dot(v1,v2)
        	anglerad = np.arccos(numdot/(self.leng(v1)*self.leng(v2)))
	        angledeg = anglerad*180/np.pi
        	return (angledeg)

	def leng(self,v):
	        length = np.linalg.norm(v)
        	return length

	def extract_geometry(self):
	        '''
        	Reads the FHI-aims output file and builds the structure's geometry
	        Returns: True if successful
	        '''
        	aims_out = open(os.path.join(self.working_dir, 'aims.out'))
	        while True:
        	    line = aims_out.readline()
	            if not line: return False  # not converged
        	    if 'Final atomic structure' in line: break
	        aims_out.readline()  # skip one line
        	aims_out.readline()  # skip one line
	        aims_out.readline()  # skip one line
        	aims_out.readline()  # skip one line
	        atom_string = ''
        	while True:
	            line = aims_out.readline()
        	    if 'Fractional coordinates:' in line: break
	            else: atom_string += line
	        self.result_struct.build_geo_whole_atom_format(atom_string)
        	return True

	def is_successful(self):
	        '''
        	checks if relaxation/optimization was successful
	        '''
        	aims_out = os.path.join(self.working_dir, 'aims.out')
		aims_out = open(aims_out,"r")
	        while True:
        	    line = aims_out.readline()
	            if line == '':
        	        break
	            if 'Have a nice day' in line:
			aims_out.close()
        	        return True
		aims_out.close()
	        return False


	def get_structure(self):
		'''
		Finds a structure currently under the tmp directory 
		Places a lock on it for use
		'''
		list_of_files = [name for name in os.listdir(self.tmp) \
                                if os.path.isfile(os.path.join(self.tmp,name)) and not os.path.isfile(os.path.join(self.tmp,name+'.lock')) and name[len(name)-5:]!='.lock' and name[len(name)-5:]!=".more"]
		found=False
		while len(list_of_files)>0:
			chosen=int(random.uniform(0,len(list_of_files)))
			name=list_of_files[chosen]
			if not os.path.isfile(os.path.join(self.tmp,name+".lock")):
				f=open(os.path.join(self.tmp,name+".lock"),"w")
				f.write("locked")
				f.close()
				os.system("chmod g=u "+os.path.join(self.tmp,name+".lock"))
				found=True
				self.filename=name
				break
			else:
				list_of_files=list_of_files[:chosen]+list_of_files[chosen+1:]
		return found

	def set_permission(self):
		'''
		Allow other group users to access this folder
		'''
		os.system("chmod -R g=u "+self.working_dir)


def get_random_index(seed=None):
    LENGTH_OF_INDEX = 10
    return sha1(repr(time.time())+str(seed)).hexdigest()[:LENGTH_OF_INDEX]

def clean_directory(path):
	'''
	Clean a directory by removing it and recreating it
	'''
	if "*" in path:
		raise RuntimeError("Clean directory is made to clean a directory with a * symbol in it; will not perform for safety")
		return
	if os.path.exists(path):
		shutil.rmtree(path)
	make_directory(path)

def make_directory(path):
	'''
	Creates a directory
	'''
	try:
		os.makedirs(path)
		os.system("chmod -R g=u "+path)
	except:
		pass

def remove_directory(path):
	'''
	Removes a directory
	'''
	if os.path.exists(path):
		shutil.rmtree(path)

def write_log(message,file="Logfile",time_print=True):
	'''
	Writes to the logfile
	'''
	f=open(file,"a")
	if time_print:
		message=time.strftime("%Y-%m-%d %H:%M:%S")+" "+message
	f.write(message+"\n")
	f.close()
	os.system("chmod g=u "+file)

def write_active(folder):
	'''
	Write to a active.info file the current time
	'''
	f=open(os.path.join(folder,"active.info"),"w")
	f.write(str(time.time()))
	f.close()
	os.system("chmod g=u "+os.path.join(folder,"active.info"))

def read_active(folder):
	'''
	Reads the active.info in the folder and returns the time
	'''
	try:
		f=open(os.path.join(folder,"active.info"))
		last_time=float(f.read())
		return last_time	
	except:
		write_log("Unable to read active.info")
		return False

def active_sleep(total_time,folder,write_active_interval=10):
	'''
	Sleep for the total time but continues to write_active in the folder
	'''
	times = int(total_time/write_active_interval)
	for i in range (times):
		write_active(folder)
		time.sleep(write_active_interval)
	write_active(folder)
	try:
		time.sleep(total_time-write_active_interval*times)
		write_active(folder)
	except:
		pass
	

def bk_folder(fdir,folder,bk_path,naming_scheme="original",nname=get_random_index()):
	'''
	Backs up a folder in the bk_folder
	'''
	if naming_scheme=="original":
		nname=folder
	elif naming_scheme=="random":
		nname=get_random_index(seed=folder)
	elif naming_scheme=="const":
		nname=nname
	else:
		raise ValueError("Unknown naming_scheme in bk_folder")

	try:
		shutil.copytree(os.path.join(fdir,folder),os.path.join(bk_path,nname))
		os.system("chmod -R g=u "+os.path.join(bk_path,nname))
		write_log("Folder %s backed-up to folder %s" % (os.path.join(fdir,folder),os.path.join(bk_path,nname)))
	except:
		write_log("bk failure when copying folder %s to folder %s" % (os.path.join(fdir,folder),os.path.join(bk_path,nname)))

def get_execute_clearance(request_folder="",buffer=3,max_wait=1000):
	'''
	Reads the execute.info in the working directory and gets clearance for executing commands such as runjob and mpirun
	'''
	for i in range (max_wait):
		with FileLock("execute.info",60):
			if not os.path.exists(os.path.join(cwd,"execute.info")):
				data_file=open(os.path.join(cwd,"execute.info"),"w")
				data_file.write(str(time.time()))
				data_file.close()
				os.system("chmod g=u "+os.path.join(cwd,"execute.info"))
				return True

			data_file=open(os.path.join(cwd,"execute.info"))
			last_time=float(data_file.read())
			data_file.close()
			current_time=time.time()
			if (current_time-last_time>buffer) or (i==max_wait):
                                data_file=open(os.path.join(cwd,"execute.info"),"w")
                                data_file.write(str(time.time()))
                                data_file.close()
				os.system("chmod g=u "+os.path.join(cwd,"execute.info"))
				if i==1000 and current_time-last_time<=buffer:
					write_log("Error in obtaining clearance for execution. Proceed with caution!")
					return False
                                return True
		print "buffering taking place in get_execute_clearance"
		if request_folder!="":
			write_active(request_folder)
		time.sleep(buffer)

if __name__ == '__main__':
	main()


