'''
Created by Patrick Kilecdi

Carries out Harris Approximation calculation on a list of structuree
'''

from core import structure, structure_handling, structure_comparison
from external_libs import matrix_op
#from external_libs.aimsrotate import rotate
#from aimsutils.rotate import rotate2
from aimsrotate import rotate
import shutil, os, copy, numpy
from evaluation import run_aims, utility
from utilities import write_log, misc, parallel_run
import multiprocessing
import time

def harris_approximation_batch(inst):
	'''
	Pulls in the structure from the working directory
	And calls harris_approximation_single to evaluate each structure
	'''
	info_level = inst.get_info_level()
	sname = "harris_approximation_batch"

	if info_level>=3:
		write_log.write_master_log(inst,"Harris approximation reporting",
		additional_item=[[sname,"structure_dir"]])

	if inst.has_option("harris_approximation_batch","working_dir"):
		working_dir = inst.get(sname,"working_dir")
	else:
		working_dir = os.path.join(
		inst.get("Genarris_master","working_dir"),"tmp")
	working_dir = os.path.abspath(working_dir)
	if not os.path.isdir(working_dir):
		os.makedirs(working_dir)

	structure_dir,structure_suffix,structure_format,exec_control_path = \
	inst.get_keywords_single_section(sname,["structure_dir",
	["structure_suffix",None],"structure_format","exec_control_path"])

	structure_dir_depth = inst.get_keywords_single_section(
	sname,[["structure_dir_depth",0]],eval=True)[0]

	if inst.has_option("harris_approximation_batch","processes_limit"):
		processes_limit = inst.get_eval(sname,"processes_limit")
	else:
		processes_limit = inst.get_eval("Genarris_master","processes_limit")
	
	if not os.path.isdir(structure_dir):
		write_log.write_master_err_and_raise(inst,
		"Harris approximation batch structure directory not found")
	
	structure_dir_copy = inst.get_boolean(sname,"structure_dir_copy")
	if structure_dir_copy: 
	#Will copy a version of the structure directory 
	#b/c once the structure is evaluated, it is deleted from the directory
		new_structure_dir = os.path.join(working_dir,
		"tmp_"+os.path.basename(structure_dir))
		try:
			shutil.copytree(structure_dir,new_structure_dir)
		except:
			pass
		if info_level>=1:
			write_log.write_master_log(inst,
			"Structures copied from folder %s to folder %s" 
			% (structure_dir,new_structure_dir))

		#Replica don't need to copy again
		inst.remove_option(sname,"structure_dir_copy") 
		inst.set(sname,"structure_dir",new_structure_dir)
		structure_dir = new_structure_dir

	list_of_struct_path = misc.list_directory_file(structure_dir,
	suffix=structure_suffix,depth=structure_dir_depth)
	spread_processes_enabled = inst.get_boolean(sname,
				"spread_processes_across_nodes")

	if spread_processes_enabled:
		nodelist = inst.get_nodelist()
	if spread_processes_enabled and len(nodelist)>1: #If there are more than 1 nodes, will launch additional Python instances on the other nodes 
		
		nodes_and_insts = []
		sample_inst = copy.deepcopy(inst)
		sample_inst.set("Genarris_master",
				"procedures",
				str(["Harris_Approximation_Batch"]))

		sample_inst.set(sname,"processes_limit",
				str(int(processes_limit/len(nodelist))))

		processes_limit = int(processes_limit/len(nodelist)) 
		#This will influence the rest of the running on this node

		sample_inst.remove_option(sname,"spread_processes_across_nodes") 
		#Remove the option so that spawned processes 
		#will not continue to spawn more ssh

		inst.remove_option(sname,"spread_processes_across_nodes") 
		#Remove the option so that spawned multiprocesses 
		#will not continue to spawn more ssh

		nodelist.remove(inst.get_master_node())
		for node in nodelist:
			new_inst = copy.deepcopy(sample_inst)
			new_inst.set("Genarris_master","master_node",node)
			nodes_and_insts.append([node,new_inst])
		if info_level>=4:
			write_log.write_master_log(inst,"This is nodelist "
						   + str(nodelist))

		mods = inst.get_with_default("parallel_settings",
					     "modules_launch",[],True)
		addc = inst.get_with_default("parallel_settings",
					     "additional_commands","''",True)
		subp = parallel_run.ssh_spawn_python_subprocess(nodes_and_insts, 
								mods,addc)
		if info_level >= 1:
			write_log.write_master_log(inst,sname+" ssh spawned"+
			" additional python subprocesses at "+str(nodelist))
	
	if processes_limit > 1:	#Will launch multiprocessing on this node
		new_inst = copy.deepcopy(inst)
		new_inst.set("harris_approximation_batch","processes_limit","1")
		inst_list = []
		for i in range (processes_limit-1):
			inst_list.append(copy.deepcopy(new_inst))

		mulp = multiprocessing.Pool(processes=processes_limit-1)
		mulp.map_async(harris_approximation_batch,inst_list)
		if info_level >= 1:
			write_log.write_master_log(inst,sname + 
			" multiprocessed additional %i processes at node %s " 
			% (processes_limit-1,inst.get_master_node()))
	
	sample_inst = copy.deepcopy(inst)
	sample_inst.remove_section("harris_approximation_single")
	sample_inst.add_section("harris_approximation_single")
	sample_inst.transfer_keywords("harris_approximation_batch",
				      "harris_approximation_single",
				       ["exec_control_path",
					"structure_output_dir",
					"energy_name",
					"molecule_harris_folders",
					"molecule_name",
					"enantiomer_name",
					"molecule_path",
					"molecule_format",
					"rotations_title",
					"structure_format",
					"full_reverse_napm",
					"full_reverse_prep_control",
					"match_molecule_tolerance",
					"mpirun_processes",
					"update_poll_interval",
					"update_poll_times"])

	if inst.has_option(sname,"clean_up"):
		if inst.get(sname,"clean_up")=="5":
			sample_inst.set("harris_approximation_single","clean_up","4")
		else:
			sample_inst.set("harris_approximation_single","clean_up",
					inst.get(sname,"clean_up"))
	
	#Begin linear, continuous evaluation of structures
	struct_path = misc.get_and_lock_file(structure_dir,
					     suffix=structure_suffix,
					     depth=structure_dir_depth)

	list_of_folders = []
	while struct_path!=False:
		if info_level>=2:
			 write_log.write_master_log(inst,sname + " beginning "+
						    "to evaluate structure " + 
						    struct_path)
		start_time = time.time()
			
		struct_name = os.path.basename(struct_path)
		struct_name = struct_name[:len(struct_name)-len(structure_suffix)]
		
		new_inst = copy.deepcopy(sample_inst)
		single_working_dir = os.path.join(working_dir,
					"harris_" + struct_name +
					"_" + misc.get_random_index())
		#Adding a random index to avoid clashing

		new_inst.set_keywords_single_section("harris_approximation_single",
		[["structure_path",struct_path],["working_dir",single_working_dir]])
		list_of_folders.append(single_working_dir)

		if inst.has_option(sname,"structure_output_dir"):
			new_p = os.path.join(inst.get(sname,"structure_output_dir"),
					     struct_name+".json")
			new_inst.set("harris_approximation_single",
				     "structure_output_path",new_p)

		if harris_approximation_single(new_inst):
			try: #Just in case a structure got double booked
				os.remove(struct_path)
			except:
				pass
			try:
				os.remove(struct_path+".lock")
			except:
				pass
		
		if info_level>=2:
			write_log.write_master_log(inst,sname+" completed with"+
			" structure %s. Time taken: %f seconds" 
			% (struct_path,time.time()-start_time))


		struct_path = misc.get_and_lock_file(structure_dir,
						     suffix=structure_suffix,
						     depth=structure_dir_depth)


	#Clean-up
	if inst.has_option("harris_approximation_batch","clean_up") \
	and inst.get("harris_approximation_batch","clean_up")=="5":
		for folder in list_of_folders:
			shutil.rmtree(folder)

	if info_level >= 2:
		write_log.write_master_log(inst,
		"harris_approximation_batch calculation completed.")
	
	if processes_limit > 1: #Joining the multiprocessing processe
		if info_level >= 1:
			write_log.write_master_log(inst,
			"harris_approximation_batch beginning to join multiprocesses")
		mulp.close()
		mulp.join()
		if info_level >= 1:
			write_log.write_master_log(inst,
			"harris_approximation_batch multiprocesses joined")

	if spread_processes_enabled and (len(nodelist)>1 \
	or nodelist[0]!=inst.get_master_node()):
		if info_level >= 1:
			write_log.write_master_log(inst,\
			"harris_approximation_batch beginning to join ssh processes")
		for sp in subp:
			sp.wait()
		if info_level >= 1:
			write_log.write_master_log(inst,\
			"harris_approximation_batch ssh subprocesses joined")


def harris_approximation_single(inst):
	'''
	Prepares and runs one instance of Harris Approximation
	'''
	info_level = inst.get_info_level()
	if info_level>=2 or (info_level>=1 and inst.procedure_exist("Harris_Approximation_Single")):
		write_log.write_master_log(inst,"Harris Approximation Single module launched",additional_item=[["harris_approximation_single","working_dir"],["harris_approximation_single","structure_path"]])

	structure_path = inst.get("harris_approximation_single","structure_path")
	working_dir, rotations_title, molecule_name, enantiomer_name = \
	inst.get_keywords_single_section("harris_approximation_single",[["working_dir",os.path.dirname(structure_path)],["rotations_title","molecule"],"molecule_name","enantiomer_name"])
	working_dir = os.path.abspath(working_dir)
	if not os.path.isdir(working_dir):
		os.makedirs(working_dir)

	
	struct = misc.load_structure_with_inst(inst,"harris_approximation_single","structure_format",structure_path)
	if struct == False:
		write_log.write_master_err(inst,"harris_approximation_single load structure failure. Structure path: "+structure_path)
		return False

	folder_copy = True
	full_reverse = False
	if not ("molecule_harris_info" in struct.properties and isinstance(struct.properties["molecule_harris_info"],list)):
	#Prepare the molecule_harris_info for rotations.in if the structure does not already have it
		success = False
		match_molecule_tolerance = inst.get_keywords([["harris_approximation_single","match_molecule_tolerance",1]],eval=True)[0]
		#First attempt to map the standard molecule
		if inst.has_option("harris_approximation_single","molecule_path"):
			molecule = misc.load_structure_with_inst(inst,"harris_approximation_single","molecule_format",inst.get("harris_approximation_single","molecule_path"))
			success = reverse_harris_info_prep(struct,molecule,inst.get_master_err_path(),match_molecule_tolerance)

			if success and (info_level >= 2 or (info_level>=1 and inst.procedure_exist("Harris_Approximation_Single"))):
				write_log.write_master_log(inst,"Successfully maps standard molecule to the structure: "+structure_path)

			if not success and not inst.has_option("harris_approximation_single","full_reverse_napm"):
				write_log.write_master_err(inst,"harris_approximation_single preparing structure's molecule_harris_info failure. Structure path: "+structure_path)
				return False

		#Then attempt to conduct full reverse scheme
		if not success and inst.has_option("harris_approximation_single","full_reverse_napm"):
			#Begin by taking the first molecule from the structure and attempts to map the rest of it by that molecule
			full_reverse = True
			napm = inst.get_eval("harris_approximation_single","full_reverse_napm")
			molecule = copy.deepcopy(struct)
			molecule.geometry = molecule.geometry[:napm]
			if "lattice_vector_a" in molecule.properties: #Strike out the lattice vectors to make single molecule
				del molecule.properties["lattice_vector_a"]
			molecule_com = structure_handling.cm_calculation(molecule,range(len(molecule.geometry)))
			structure_handling.cell_translation(molecule,[-x for x in molecule_com],False)

			success = reverse_harris_info_prep(struct,molecule,inst.get_master_err_path(),match_molecule_tolerance)
			if not success:
				write_log.write_master_err(inst,"harris_approximation_single failed to conduct full reverse scheme on structure: "+structure_path)
				return False
			if info_level >= 2 or (info_level>=1 and inst.procedure_exist("Harris_Approximation_Single")):
				write_log.write_master_log(inst,"Successfully conducted the mapping part of full reverse scheme on structure: "+structure_path)


			f = open(os.path.join(working_dir,"single_molecule.in"),"w")
			f.write(molecule.get_geometry_atom_format())
			f.close()

			inst.grant_permission(working_dir)
			#Now begins to prepare (by default) the single molecule calculations necessary for Harris rotation

			inst_1 = copy.deepcopy(inst)
			inst_1.remove_section("harris_single_molecule_prep")
			inst_1.add_section("harris_single_molecule_prep")
			inst_1.set("harris_single_molecule_prep","molecule_path",os.path.join(working_dir,"single_molecule.in"))
			inst_1.set("harris_single_molecule_prep","molecule_format","geometry")
			inst_1.set("harris_single_molecule_prep","control_path",inst.get("harris_approximation_single","full_reverse_prep_control"))
			inst_1.set("harris_single_molecule_prep","output_dir",os.path.join(working_dir,molecule_name))
			inst_1.set("harris_single_molecule_prep","enantiomer_output_dir",os.path.join(working_dir,enantiomer_name))
			inst_1.transfer_keywords("harris_approximation_single","harris_single_molecule_prep",["mpirun_processes","update_poll_interval","update_poll_times"])
			harris_single_molecule_prep(inst_1)
			folder_copy = False #No need to copy the folder into the current working directory
			inst.grant_permission(working_dir)
			if info_level >= 2 or (info_level>=1 and inst.procedure_exist("Harris_Approximation_Single")):
				write_log.write_master_log(inst,"Successfully conducted first molecule preparation for structure: "+structure_path)
		else:
			write_log.write_master_log(inst,"Harris Approximation failed on structure: %s; not enough information provided" % structure_path)
			write_log.write_master_err(inst,"Harris Approximation failed on structure: %s; not enough information provided" % structure_path)
			return False


	rotations_write(struct,molecule_name,enantiomer_name,molecule_name=rotations_title,fname=os.path.join(working_dir,"rotations.in"))
	if not full_reverse: 
		#If full reverse scheme is called, molecule harris folder would already be located in the respective temporary calculation folder
		molecule_harris_folders = inst.get_eval("harris_approximation_single","molecule_harris_folders")

		if folder_copy: #No need to copy if the full reverse scheme was used
			for folder in molecule_harris_folders:
			#Moving all the necessary folders containing single molecule information to the working_directory
				if not os.path.isdir(os.path.abspath(folder)):
					write_log.write_master_err_and_raise(inst,"molecule_harris_folder missing. Path = "+os.path.abspath(folder))
				misc.safe_copy_folder(folder,working_dir)

	inst.grant_permission(working_dir)
	if info_level>=2 or (info_level>=1 and inst.procedure_exist("Harris_Approximation_Single")):
		write_log.write_master_log(inst,"Harris_Approximation_Single working directory preparation completes. About to launch HA core",additional_item=[["harris_approximation_single","working_dir"],["harris_approximation_single","structure_path"]])

	#Begin launching harris_single_run to carry out the core calculation
	inst = copy.deepcopy(inst)
	inst.remove_section("harris_single_run")
	inst.add_section("harris_single_run")
	inst.set("harris_single_run","working_dir",working_dir)
	inst.set("harris_single_run","control_path",inst.get("harris_approximation_single","exec_control_path"))
	inst.transfer_keywords("harris_approximation_single","harris_single_run",["mpirun_processes","update_poll_interval","update_poll_times"])
	harris_single_run(inst)

	if info_level>=2 or (info_level>=1 and inst.procedure_exist("Harris_Approximation_Single")):
		write_log.write_master_log(inst,"Harris_Approximation_Single module HA core running complete. Beginning postprocessing.",additional_item=[["harris_approximation_single","working_dir"],["harris_approximation_single","structure_path"]])



	if inst.has_option("harris_approximation_single","structure_output_path"): #Begins extracting energy if a structure output is required
		harris_energy = utility.aims_extract_energy(os.path.join(os.path.join(working_dir,"harris"),"aims.out"))
		if harris_energy == False:
			write_log.write_master_err(inst,"Harris Approximation failed to yield an energy.",additional_item=[["harris_approximation_single","working_dir"],["harris_approximation_single","structure_path"]])
		energy_name = inst.get_keywords([["harris_approximation_single","energy_name","energy"]])[0]
		struct.set_property(energy_name,harris_energy)

		structure_output_path = inst.get("harris_approximation_single","structure_output_path")
		try:
			os.makedirs(os.path.dirname(structure_output_path))
		except:
			pass
		f = open(structure_output_path,"w")
		f.write(struct.dumps())
		f.close()
		inst.grant_permission(structure_output_path)

	#Now initiate clean up of harris folders

	clean_up = inst.get_keywords([["harris_approximation_single","clean_up",1]],eval=True)[0]
	if clean_up not in [0,1,2,3,4]:
		write_log.write_master_err(inst,"Unknown clean up type for harris_approximation_single")
	else:
		if clean_up >= 1 and not full_reverse: #Removes the molecule folders
			for folder in molecule_harris_folders:
				if os.path.dirname(os.path.abspath(folder))!=working_dir:
					shutil.rmtree(os.path.join(working_dir,os.path.basename(folder)),ignore_errors=True)
		harris_dir = os.path.join(working_dir,"harris")
		if clean_up == 2 and os.path.isfile(os.path.join(harris_dir,"restart.new000")): #Additional removes restart.neww000
			os.remove(os.path.join(harris_dir,"restart.new000"))
		if clean_up == 3:
		#Removes all but aims.out, geometry.in and control.in
			files = os.listdir(harris_dir)
			for name in files:
				if name!="aims.out" and name!="geometry.in" and name!="control.in":
					os.remove(os.path.join(harris_dir,name))
		if clean_up == 4: #Removes the entire folder
			shutil.rmtree(harris_dir,ignore_errors=True)

	return True	
	
	
def harris_single_run(inst):
	'''
	Runs a single instance of Harris Approximation
	Note: this function merely carries out the calling of aimsrotate and aims binary to evaluate the rotated geometry
	For load-balancing reasons, this function is kept at a relatively low position in hierarchy
	Other functions are to call this to finish up the Harris Approximation
	'''
	working_dir, control_path = inst.get_keywords_single_section("harris_single_run",["working_dir","control_path"])
	current_dir = os.getcwd()
	os.chdir(working_dir)

	#First create the harris folder with the rotated geometry
	rotated_dir = os.path.join(working_dir,"harris")
#	r = rotate2.Rotations()
#	r.load_rotations("rotations.in")
#	r.write_restartfile("harris")
	R = rotate.RotateMixed(working_dir+"/")
	R.write_files(folder=rotated_dir)

	os.chdir(current_dir)


	inst.grant_permission(working_dir)


	# Now calls aims_single_run to run in the folder created
	if os.path.abspath(os.path.join(rotated_dir,"control.in"))!=os.path.abspath(control_path):
		shutil.copyfile(control_path,os.path.join(rotated_dir,"control.in"))
	inst = copy.deepcopy(inst)
	inst.remove_section("aims_single_run")
	inst.add_section("aims_single_run")
	inst.set("aims_single_run","working_dir",rotated_dir)
	inst.transfer_keywords("harris_single_run","aims_single_run",["mpirun_processes","update_poll_times","update_poll_interval"])
	inst.set("aims_single_run","mpirun_hosts","['"+inst.get("Genarris_master","master_node")+"']")
	run_aims.aims_single_run(inst)	
	
def harris_single_molecule_prep(inst):
	'''
	Prepares /a single molecule to be used for Harris approximation	
	'''
	info_level = inst.get("Genarris_master","info_level")
	sname = "harris_single_molecule_prep"
	if info_level>=2 or (info_level>=1 and 
	inst.procedure_exist("Harris_Single_Molecule_Prep")):
		message = ("Harris single molecule preparation "+
		"beginning for original molecule.")

		write_log.write_master_log(inst,message,
		additional_item=[[sname,"output_dir"]])

	molecule_path, control_path = \
	inst.get_keywords_single_section(sname,
					["molecule_path",
					"control_path"])

	#First move the molecule to output_dir
	struct = misc.load_structure_with_inst(inst, sname,
					      "molecule_format",molecule_path)

	if "lattice_vector_a" in struct.properties: 
	#So that lattice vectors are not printed out in geometry
		del struct.properties["lattice_vector_a"]
	
	#Now places the molecule's COM at the origin
	molecule_com = structure_handling.cm_calculation(struct,
					range(len(struct.geometry)))

	structure_handling.cell_translation(struct,
					    [-x for x in molecule_com],
					    False)

	run_single_molecule_prep(inst,struct)
	if info_level>=2 or (info_level>=1 and 
	inst.procedure_exist("Harris_Single_Molecule_Prep")):
		message = ("Harris single molecule preparation "+
		"completed for original molecule.")

		write_log.write_master_log(inst,message,
		additional_item=[[sname,"output_dir"]])

	if inst.has_option(sname,"enantiomer_output_dir"):
		if info_level>=2 or (info_level>=1 and 
		inst.procedure_exist("Harris_Single_Molecule_Prep")):
			message = ("Harris single molecule preparation "+
			"moving on to molecule's enantiomer.")

			write_log.write_master_log(inst,message,
			additional_item=[[sname,"enantiomer_output_dir"]])

		structure_handling.cell_reflection_z(struct,False)
		run_single_molecule_prep(inst,struct,True)	
		if info_level>=2 or (info_level>=1 and 
		inst.procedure_exist("Harris_Single_Molecule_Prep")):
			message = ("Harris single molecule preparation "+
			"completed for molecule's enantiomer")
			
			write_log.write_master_log(inst,message,
			additional_item=[[sname,"enantiomer_output_dir"]])

def run_single_molecule_prep(inst,struct,enantiomer_run=False):
	'''
	Sends aims single run with the correct inst to carry out the single molecule run
	'''
	if enantiomer_run==False:
		working_dir = inst.get("harris_single_molecule_prep","output_dir")
	else:
		working_dir = inst.get("harris_single_molecule_prep","enantiomer_output_dir")
	try:
		os.makedirs(working_dir)
	except:
		pass
	f = open(os.path.join(working_dir,"geometry.in"),"w")
	f.write(struct.get_geometry_atom_format())
	f.close()
	if inst.get("harris_single_molecule_prep","control_path")!=os.path.join(working_dir,"control.in"):
		shutil.copyfile(inst.get("harris_single_molecule_prep","control_path"),os.path.join(working_dir,"control.in"))

	#Now calls run_aims.aims_single_run to run the job
	inst = copy.deepcopy(inst)
	if inst.has_section("aims_single_run"):
		inst.remove_section("aims_single_run")
	inst.add_section("aims_single_run")
	inst.set("aims_single_run","working_dir",working_dir)
	inst.transfer_keywords("harris_single_molecule_prep","aims_single_run",["update_poll_interval","update_poll_times","mpirun_processes","additional_arguments"])

	run_aims.aims_single_run(inst)

def reverse_harris_approximation(inst):
	'''
	Finds the rotational angle of the other molecules in the structure from the first one
	'''
	structure_path = inst.get("reverse_harris_approximation","structure_path")
	if inst.has_option("reverse_harris_approximation","working_dir"):
		working_dir = inst.get("reverse_harris_approximation","working_dir")
	else:
		working_dir = os.path.dirname(structure_path)

	original_dir = os.path.join(working_dir,"molecule_original")
	enantiomer_dir = os.path.join(working_dir,"molecule_enantiomer")
	
	nmpc = int(inst.get("reverse_harris_approximation","NMPC"))
	napm = int(inst.get("reverse_harris_approximation","NAPM"))
	prep_control_path = inst.get("reverse_harris_approximation","prep_control_path")
	exec_control_path = inst.get("reverse_harris_approximation","exec_control_path")

	
	struct = misc.load_structure_with_inst(inst,"reverse_harris_approximation","structure_format",structure_path)
	if struct == False:
		write_log.write_master_err(inst,"reverse_harris_approximation unable to pull in the structure "+structure_path)
		return False

	if nmpc*napm!=len(struct.geometry):
		write_log.write_master_err(inst,"reverse_harris_approximation: nmpc and napm doees not agree with the length of the structure "+structure_path)
	
	if not os.path.exists(working_dir):
		os.makedirs(working_dir)


	first_molecule = copy.deepcopy(struct)
	first_molecule.geometry = first_molecule.geometry[:napm]
	if "lattice_vector_a" in first_molecule.properties: #Strike out the lattice vectors to make single molecule
		del first_molecule.properties["lattice_vector_a"]
	molecule_com = structure_handling.cm_calculation(first_molecule,range(len(first_molecule.geometry)))
	structure_handling.cell_translation(first_molecule,[-x for x in molecule_com],False)

	#Prepare the harris rotational information here
	if not reverse_harris_info_prep(struct,first_molecule,inst.get("Genarris_master","master_err_path")):
		write_log.write_master_err(inst,"reverse_harris_approximation failed to prepare reverse_harris_information for structure"+structure_path)
		return False		

	f = open(os.path.join(working_dir,"single_molecule.in"),"w")
	f.write(first_molecule.get_geometry_atom_format())
	f.close()

	inst.grant_permission(working_dir)
	#Now begins to prepare (by default) the single molecule calculations

	inst_1 = copy.deepcopy(inst)
	inst_1.remove_section("harris_single_molecule_prep")
	inst_1.add_section("harris_single_molecule_prep")
	inst_1.set("harris_single_molecule_prep","molecule_path",os.path.join(working_dir,"single_molecule.in"))
	inst_1.set("harris_single_molecule_prep","molecule_format","geometry")
	inst_1.set("harris_single_molecule_prep","control_path",prep_control_path)
	inst_1.set("harris_single_molecule_prep","output_dir",os.path.join(working_dir,"molecule_original"))
	inst_1.set("harris_single_molecule_prep","enantiomer_output_dir",os.path.join(working_dir,"molecule_enantiomer"))
        inst_1.set("harris_single_molecule_prep","update_poll_interval", "1")
        inst_1.set("harris_single_molecule_prep","update_poll_times", "120")
	harris_single_molecule_prep(inst_1)

	inst.grant_permission(working_dir)
	# Now begins to call aimsrotate to rotate the file and run the harris approximation
	
	rotations_write(struct,"molecule_original","molecule_enantiomer",fname=os.path.join(working_dir,"rotations.in"))
	inst_1 = copy.deepcopy(inst)
	inst_1.remove_section("harris_single_run")
	inst_1.add_section("harris_single_run")
	inst_1.set("harris_single_run","working_dir",working_dir)
	inst_1.set("harris_single_run","control_path",exec_control_path)
        inst_1.set("harris_single_run","update_poll_interval", "1")
        inst_1.set("harris_single_run","update_poll_times", "120")
	harris_single_run(inst_1)

def reverse_harris_approximation_batch(inst):
	'''
	Set up folders and do batch reverse_harris_approximation
	'''
	if inst.has_option("reverse_harris_approximation_batch","working_dir"):
		working_dir = inst.get("reverse_harris_approximation_batch","working_dir")
	else:
		working_dir = inst.get("Genarris_master","working_dir")
	structure_suffix, structure_format, prep_control_path, exec_control_path,nmpc,napm =\
	inst.get_keywords_single_section("reverse_harris_approximation_batch",["structure_suffix",["structure_format",None],"prep_control_path","exec_control_path","NMPC","NAPM"])
	
	#Now calls reverse_harris_approximation to handle each individual structure
	list_of_files = [name for name in os.listdir(working_dir) if name[len(name)-len(structure_suffix):]==structure_suffix]
#	current_dir = os.getcwd()
	for name in list_of_files:
		inst_1 = copy.deepcopy(inst)
		inst_1.remove_section("reverse_harris_approximation")
		inst_1.add_section("reverse_harris_approximation")
		inst_1.set_keywords_single_section("reverse_harris_approximation",[["working_dir",os.path.join(working_dir,"reverse_harris_"+name[:len(name)-len(structure_suffix)])],\
		["structure_path",os.path.join(working_dir,name)],["prep_control_path",prep_control_path],["exec_control_path",exec_control_path],["NMPC",nmpc],["NAPM",napm]],clear_first=True)
		if structure_format!=None:
			inst_1.set("reverse_harris_approximation","structure_format",structure_format)
		reverse_harris_approximation(inst_1)

def rotations_write(struct,original_name,enantiomer_name,molecule_name="molecule",fname=""):
	'''
	Taking a structure with the following information
	struct.properties["molecule_harris_info"][#molecule]=[Tx,Ty,Tz,angle_1,angle_2,angle_3,is_enantiomer]
	Outputs a proper rotations.in file
	'''
	result = "title %s\n" % molecule_name
	if "lattice_vector_a" in struct.properties: #Will report error if lattice_vector_a exists but not the rest
		result+="lattice_vector %f %f %f\n" % (struct.properties["lattice_vector_a"][0],struct.properties["lattice_vector_a"][1],struct.properties["lattice_vector_a"][2])
		result+="lattice_vector %f %f %f\n" % (struct.properties["lattice_vector_b"][0],struct.properties["lattice_vector_b"][1],struct.properties["lattice_vector_b"][2])
		result+="lattice_vector %f %f %f\n" % (struct.properties["lattice_vector_c"][0],struct.properties["lattice_vector_c"][1],struct.properties["lattice_vector_c"][2])
	for info in struct.properties["molecule_harris_info"]:
		if info[6]:
			choice = enantiomer_name
		else:
			choice = original_name
		result += "%f %f %f %f %f %f %s\n" % (info[0],info[1],info[2],info[3],info[4],info[5],choice)

	if fname == '':
		try:
			print(result)
		except:
			pass
		return result
	else:
		fout = open(fname,'w')
		fout.write(result)
		fout.close()

def reverse_harris_info_prep(struct,molecule,error_path=None, match_molecule_tolerance=1):
	'''
	Prepares the struct.properties["molecule_harris_info"]
	'''
	if len(struct.geometry)%len(molecule.geometry)!=0:
		raise ValueError("structure's atom sum is not an interger multiple of the molecule's atom sum; len(struct.geometry)=%i and len(molecule.geometry)=%i" % (len(struct.geometry),len(molecule.geometry)))
	struct.properties["molecule_harris_info"]=[]
	for i in range(len(struct.geometry)/len(molecule.geometry)):
		new_molec = copy.deepcopy(struct)
		new_molec.geometry = new_molec.geometry[i*len(molecule.geometry):(i+1)*len(molecule.geometry)]
		struct.properties["molecule_harris_info"].append(reverse_harris_match_molecule(new_molec,molecule,error_path,match_molecule_tolerance))
	if False in struct.properties["molecule_harris_info"]:
		return False
	return True

def reverse_harris_match_molecule(s1,s2,error_path=None,match_molecule_tolerance=1):
	'''
	s1 is the molecule to be mapped
	s2 is the standard molecule, whose COM lies at the origin
	Mapping is successful if the residual is smaller than match_molecule_tolerance
	See core.structure_comparison.residual_single_molecule for details about residual calculation	
	Outputs translation and Euler angles in rzyz form if mapping successful
	WARNING: it is crucial for the origin of s2 to lie at origin for this to work
	'''
	if len(s1.geometry)!=len(s2.geometry):
		raise ValueError("Two structures do not have the same number of atoms")
	s2 = copy.deepcopy(s2)	

	match_molecule_length_requirement=0.1 #Length required for difference vectors for calculating the rotation axis
	#Difference vectors are drawn from one atom in s1 to the same atom in s2
	match_molecule_cross_tolerance=0.001 #Length required for the rotation axis before normalization

	standard_back_up = copy.deepcopy(s2)
	napm = len(s1.geometry)

	cm1=structure_handling.cm_calculation(s1,range(napm))
	cm2=structure_handling.cm_calculation(s2,range(napm))
	trans=[cm2[j]-cm1[j] for j in range (3)]
	structure_handling.cell_translation(s1,trans,False) #Fixes the cm of two molecules to be at the same place
	resi = structure_comparison.residual_single_molecule(s1,s2,napm,0,0)
	min_resi = resi #min_resi keeps track of minimum resi for output in case mapping failed

	trans = [-x for x in trans]
	if resi<match_molecule_tolerance:
	        return [0,0,0]+trans+[False] #No rotational part

	lll=0
	while lll<2: #The second time around, the molecule will be mapped by enantiomer instead
		lll+=1
		vec1=vec2=None
		i=0
		rotation_axis=[0,0,0]
		chosen=None
		while (vec2==None) and (i<napm): 
            	#find two vectors of sufficient length that can help find the rotational axis
			diff=[s1.geometry[i][j]-s2.geometry[i][j] for j in range (3)]
			leng=numpy.linalg.norm(diff)
			if leng<match_molecule_length_requirement:
		                #AN atom is invariant under the rotation. 
		                #The rotational axis will be defined by the vector from the cm to the atom
		                #unless the atom sits on the center
				rotation_axis=[s1.geometry[i][j]-cm1[j] for j in range (3)]
				if numpy.linalg.norm(rotation_axis)>match_molecule_cross_tolerance:
					break
			if leng>match_molecule_length_requirement:
				if vec1==None:
					vec1=diff
					chosen=i
					i+=1
					continue
				if numpy.linalg.norm(numpy.cross(vec1,diff))>match_molecule_cross_tolerance:
					vec2=diff
					break
			i+=1
		if numpy.linalg.norm(rotation_axis)<match_molecule_cross_tolerance: 
		#If no atom is invariant under rotation, then the rotation axis will be the cross product of the difference vectors found above
			if vec2==None:
				vec2=[0,0,0]
			rotation_axis=numpy.cross(vec1,vec2)
		rl=numpy.linalg.norm(rotation_axis)
		if rl>match_molecule_cross_tolerance:
			for j in range (3):
				rotation_axis[j]/=rl #Normalize the rotation vector
			if chosen==None:
				chosen=0
				while (chosen<napm) and (numpy.linalg.norm(numpy.cross(rotation_axis,[s1.geometry[chosen][j] for j in range (3)]))<match_molecule_cross_tolerance):
					chosen+=1
			v1=[s1.geometry[chosen][j] for j in range (3)]
			v2=[s2.geometry[chosen][j] for j in range (3)] #COM already alligned at the origin
			m1=numpy.dot(v1,rotation_axis)
			v1=[v1[j]-m1*rotation_axis[j] for j in range (3)]
			#Find the perpendicular displacement vector from the atom to the rotational axis
			#i.e. the radius
			m2=numpy.dot(v2,rotation_axis)
			v2=[v2[j]-m2*rotation_axis[j] for j in range (3)]
			l1=numpy.linalg.norm(v1); l2=numpy.linalg.norm(v2)
			direction=numpy.cross(v1,v2)
			rad=numpy.arcsin(numpy.linalg.norm(direction)/(l1*l2))
			if numpy.sign(numpy.dot(v1,v2))==-1:
				#Greater than 90 degree rotation
				rad=numpy.pi-rad
			if numpy.linalg.norm(direction)>match_molecule_cross_tolerance:
				#If the two atoms are 180 apart, then direction will be a 0 vector
				rad*=numpy.sign(numpy.dot(direction,rotation_axis))
			#If the current direction does not point in the same direction as the rotational axis,
			#then the angle of rotation needs to be flipped

			#Now rotate the structure with the found rotation and test residual
			structure_handling.cell_rotation(s2,vec=rotation_axis,origin=[0,0,0],rad=-rad,create_duplicate=False)
			resi=structure_comparison.residual_single_molecule(s1,s2,napm,0,0)
			min_resi = min(resi,min_resi)

			if resi<match_molecule_tolerance:
				#Prepare Euler angles for output
				x,y,z = rotation_axis
				c = numpy.cos(-rad)
				s = numpy.sin(-rad)
				mat =  [[x*x*(1-c)+c,x*y*(1-c)-z*s,x*z*(1-c)+y*s],
					[x*y*(1-c)+z*s,y*y*(1-c)+c,y*z*(1-c)-x*s],
					[x*z*(1-c)-y*s,y*z*(1-c)+x*s,z*z*(1-c)+c]]
				angles = list(numpy.rad2deg(matrix_op.euler_from_matrix(mat,"rzyz")))
				if lll == 1:
					return angles+trans+[False] #Not enantiomer yet
				else:
					return angles+trans+[True]

		#Flips the standard molecule now to map instead the enantiomer onto s1
		s2 = copy.deepcopy(standard_back_up)
		structure_handling.cell_reflection_z(s2,False)
		resi=structure_comparison.residual_single_molecule(s1,s2,napm,0,0)
		min_resi = min(resi,min_resi)

		if resi<match_molecule_tolerance: #Only a mirror reflection
			return [0,0,0]+trans+[True]
	if error_path!=None:
		write_log.write_log(error_path,"reverse_harris_match_molecule: Unable to match molecule, minimum residual="+str(min_resi))
	return False


