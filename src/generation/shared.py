"""
If any part of this module is used for a publication please cite:

X. Li, F. Curtis, T. Rose, C. Schober, A. Vazquez-Mayagoitia, K. Reuter,
H. Oberhofer, and N. Marom "Genarris: Random Generation of Molecular Crystal 
Structures and Fast Screening with a Harris Approximation, ",
J. Chem. Phys., DOI: 10.1063/1.5014038; arXiv 1803.02145 (2018)
"""
'''
Created on June 17, 2015
Author: Patrick Kilecdi
'''
from core import structure, structure_handling, instruct
from generation import sgroup
import random
import numpy
from utilities import misc, write_log
import copy

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



def single_structure_generation(inst,molecule=None):
	'''
	Generates a structure according to the instruction 
	specified by the section structure_generation
	Optionally reads in a molecule, 
	which should be a structure.Structure() type, 
	to avoid re-reading molecule's geometry file
	'''
	def fill_molecule_by_wyckoff_list(struct,current_wyc):
		'''
		Recursively fills in molecule
		current_wyc is the index of wyckoff_list
		'''
		for i in range (fill_limit):

			struct_new = place_molecule_space_group(struct,
			molecule,space_group,wyckoff_list[current_wyc]) 
			#Note: this function will not alter the original structure

			if closeness_check(struct_new,inst,len(molecule.geometry)):
				if current_wyc == len(wyckoff_list)-1: 
				#Has finished filling all the Wyckoff position
					return struct_new
				else: #Move on to the next
					result = fill_molecule_by_wyckoff_list(
					struct_new,current_wyc+1)
					if result!=False: #Successful fill
						return result

		return False #Unsuccessful fill
			
	inst = copy.deepcopy(inst)
	sname = "structure_generation"

	if molecule==None: #Need to load molecule
		molecule = misc.load_structure_with_inst(inst,sname,
		"molecule_format",inst.get(sname,"molecule_path"))

	master_limit, sg_limit, lv_limit, wyc_limit, fill_limit =\
	inst.get_keywords_single_section("structure_generation_counters",\
					[["master_counter",1024],
					["space_group_counter",64],
					["lattice_vector_counter",16],
					["wyckoff_list_counter",4],
					["fill_molecule_counter",5]],True)
	#e.g., if lattice_vector_counter is set to 8,
	#then a new set of lattice_vector will be generated
	#If 8 attempts of filling the unit cell has failed
	#Due to failed closeness check

	info_level = inst.get_info_level()
	master_counter = 0
	sg_counter = sg_limit #For initializing
	lv_counter = lv_limit
	wyc_counter = wyc_limit 

	# ------------------ Initializing Space Group -------------------

	space_group, space_group_allowed, wyckoff_list, nmpc = \
	inst.get_keywords_single_section("structure_generation",
					[["space_group",-1],
					["space_group_allowed",[]],
					["wyckoff_list",[]],"NMPC"],True)

	molecule_name, enantiomer_name = \
	inst.get_keywords_single_section("structure_generation",
					[["molecule_name","molecule_original"],
					["enantiomer_name","molecule_enantiomer"]],
					False)

	if space_group!=-1: #Space group pre-specified
		sg_limit += master_limit + 1 #Prevent space group update
		sg = sgroup.Sgroup(space_group)
		sg.wyckoff_preparation(nmpc)

	if wyckoff_list!=[]: #List of Wyckoff positions pre-specified
		wyc_limit += master_limit + 1

		if space_group!=-1 and nmpc!=sg.wyckoff_counter(wyckoff_list):
			write_log.write_master_err(inst,"Space group specified"+
			" but specified Wyckoff number unsuitable for NMPC")
			raise ValueError("Space group specified but specified"+
					" Wyckoff number unsuitable for NMPC")

		if space_group == -1:
			if space_group_allowed==[]: 
			#No allowed space group specified
			#Automatically select space groups that 
			#coupled with the given wyckoff position, 
			#gives the right amount of nmpc; 
				space_group_allowed = \
				sgroup.allowed_sg_wyckoff_list(nmpc,wyckoff_list)

				if len(space_group_allowed)==0:
					message = ("wyckoff_list wpecified as %s, " +
					"nmpc set as %i; " +
					"cannot find suitable space group")\
					% (str(wyckoff_list),nmpc)

					write_log.write_master_err(inst,message)
					raise ValueError(message)

			else: 
			#Will prune away space groups allowed 
			#but does not match nmpc with the given wyckoff_list
				sg_allowed = space_group_allowed[:]
				for i in sg_allowed:
					sg = sgroup.Sgroup(i)
					if nmpc!=sg.wyckoff_counter(wyckoff_list):
						space_group_allowed.pop(
						space_group_allowed.index(i))

				if len(space_group_allowed)==0:
					message = ("With the Wyckoff list %s, "+
					"none of the allowed space group %s " +
					"adds up to nmpc %i") \
					% (str(wyckoff_list),str(sg_allowed),nmpc)
					write_log.write_master_err(inst,message)
					raise ValueError(message)

	elif space_group==-1:
	#No Wyckoff position list specified
		if space_group_allowed == []:
			space_group_allowed = sgroup.allowed_sg_nmpc(nmpc)
		else:
			check_list = sgroup.allowed_sg_nmpc(nmpc)
			space_group_allowed = [x for x in space_group_allowed 
						if x in check_list]

		if len(space_group_allowed)==0:
			message = "No allowed space group can give the nmpc"
			write_log.write_master_err(inst,message)
			raise ValueError(message)
	
	if inst.has_option(sname,"chiral_or_racemic"):
		sym = inst.get(sname,"chiral_or_racemic")
		if sym == "0":
			space_group_allowed = \
			sgroup.select_chiral_sg(space_group_allowed)

		elif sym == "1":
			space_group_allowed = \
			sgroup.select_racemic_sg(space_group_allowed)
		else:
			message = ('structure_generation.chiral_or_racemic ' +
			'should be either 0 (for chiral space groups) '+
			'or 1 (for racemic space groups)')
			write_log.write_master_err_and_raise(inst,message)
			

	if info_level >= 2:
		message = ("Upon space group initializing, " + 
			"nmpc = %i, wyckoff_list = %s, space_group = %i, " +
			"space_group_allowed = %s") % (nmpc,str(wyckoff_list),
			space_group,str(space_group_allowed))

		write_log.write_master_log(inst,message)

	#------------- Beginning Main Structure Generation -------------

	#At this point, space_group_allowed non-empty if space_group=-1
	#All the space groups in space_group allowed will make sense
	#If space_group = -1, sg will be the correct space group
	success = False
	while master_counter < master_limit:
		master_counter += 1
		sg_counter += 1
		lv_counter += 1
		wyc_counter += 1
		if sg_counter > sg_limit: #Update space group
			space_group = space_group_allowed[\
			int(random.uniform(0,len(space_group_allowed)))]

			sg = sgroup.Sgroup(space_group)
			sg.wyckoff_preparation(nmpc)
			sg_counter = 0
			lv_counter = lv_limit+1 #Need to update lattice vectors
			inst.set(sname,"space_group",str(space_group))
			if wyc_limit < master_limit: #Wyckoff list not predefined
			#Now need to update wyckoff list
				wyc_counter = wyc_limit+1
		
		if lv_counter > lv_limit: #Update lattice vectors
			struct = structure.Structure()
			lattice_vector_generation(struct,inst)
			if info_level >= 3:
				message = ("New lattice vectors " +
				"for space group %i generated as %s, "
				"with angles alpha=%f, beta=%f, gamma=%f") % \
				(space_group,
				str([struct.properties["lattice_vector_a"],
				struct.properties["lattice_vector_b"],
				struct.properties["lattice_vector_c"]]),
				struct.properties["alpha"],
				struct.properties["beta"],
				struct.properties["gamma"])

				write_log.write_master_log(inst,message)
			lv_counter = 0 
		
		if wyc_counter > wyc_limit:
			wyckoff_list = sg.wyckoff_selection(nmpc)
			wyc_counter = 0
			if info_level >= 3:
				message = ("New Wyckoff list selected as %s " + 
				"for space group %i and nmpc %i" )% \
				(str(wyckoff_list),space_group,nmpc)

				write_log.write_master_log(inst,message)
	
		success = fill_molecule_by_wyckoff_list(struct,0)
		if success!=False:
			break

	if success!=False:
		if info_level>=2:
			message = ("Successful generation of structure " + 
			"in space group %i, with wyckoff list as %s") % \
			(space_group,str(wyckoff_list))

			write_log.write_master_log(inst,message)

		success.properties["space_group"] = space_group
		success.properties["wyckoff_list"] = wyckoff_list
		success.properties["NMPC"] = nmpc
		success.properties["molecule_name"] = molecule_name
		success.properties["enantiomer_name"] = enantiomer_name
		if inst.get_boolean(sname,"orthogonalize_unit_cell"):
			success = structure_handling.cell_modification(
			success,nmpc,molecule.get_n_atoms(),True)	
		return success
	else:
		if info_level>=2:
			message = "single_molecule_generation unsuccessful"
			write_log.write_master_log(inst,message,
			additional_item = [["structure_generation","space_group"],
					["structure_generation","wyckoff_list"]])
		return False	
		
				
def place_molecule_space_group (struct,molecule,sgn,wycn,create_duplicate=True):
	'''
	Randomly places a molecule into the cell according to the symmetry requirement of the Wyckoff position
	Both structure and molecule will should be of structure.Structure() class
	'''
	if create_duplicate:
		struct = copy.deepcopy(struct)
	sg = sgroup.Sgroup(sgn)
	lattice_vectors = numpy.transpose([struct.properties["lattice_vector_a"],
					   struct.properties["lattice_vector_b"],
					   struct.properties["lattice_vector_c"]])
	latinv = numpy.linalg.inv(lattice_vectors)

	for i in range (10):
		fst_frac = sg.wycgen(wycn) 
		#This generates the fractional coordinates 
		#that fits the corresponding Wyckoff number

		fst_orien = misc.random_rotation_matrix()
	
		list_of_fracs = [] #List of translation vectors
		list_of_orien = [] #List of potential rotational matrices
	
		for ll in range (0, len(sg.btran)): 
		#Loop over Bravais centering
			for l in range (0, len(sg.op)): 
			#Loop over point group operations
				new_frac = numpy.dot(sg.op[l],fst_frac)
				new_frac = numpy.add(new_frac,sg.trans[l]) 
				#Add the Translation required by the operation
				new_frac = numpy.add(new_frac,sg.btran[ll]) 
				#Add the translation from the Bravais centering
	
				misc.frac_coor_adjust(new_frac,False)
				space_group_rotation = numpy.dot(lattice_vectors,
				numpy.dot(sg.op[l],latinv)) 
				#Fractional rotation to absolute rotation

				new_orien = numpy.dot(space_group_rotation,fst_orien)
			
				c = True
				for i in range (len(list_of_fracs)): 
				#If a special Wyckoff selection, 
				#after the symmetry operation is applied,
				#Two or more molecules will sit on the same site.
				#Here groups them into sites for later selection
					diff = numpy.linalg.norm(
					numpy.subtract(list_of_fracs[i],new_frac))

					if diff < 0.00001: #Same site
						c = False
						list_of_orien[i].append(new_orien) 
						break
				if c:
					list_of_fracs.append(new_frac)
					list_of_orien.append([new_orien])
		if len(list_of_fracs)==sg.wmult[wycn]*sg.bmult:
			break

	if len(list_of_fracs)!=sg.wmult[wycn]*sg.bmult:
		message = ("Repeated failure to match number of molecules "+
		"to that of Wyckoff position; check space group module for: "
		"space group=%i, wycn=%i" % (sg.type,wycn))
		raise RuntimeError(message) #!!!

	final_orien = []
	for candidates in list_of_orien: 
	#Pick one of the several sitting on the same site 
	#if a special Wyckoff position is selected
		chosen = int(random.uniform(0,len(candidates)))
		final_orien.append(candidates[chosen])

	final_trans = []
	for frac in list_of_fracs:
		#Convert to absolute coordinates
		final_trans.append(numpy.dot(lattice_vectors,frac))

	for i in range (len(final_trans)):
		#Now add the molecules into the structure
		k = "molecule_harris_info"
		n = misc.molecule_harris_info_prep(final_trans[i], 
						   final_orien[i])
		if k in struct.properties:
			struct.properties[k].append(n)
		else:
			struct.properties[k] = [n]

		new_mole = structure_handling.cell_transform_mat(\
		molecule,final_orien[i])

		structure_handling.cell_translation(new_mole,
						    final_trans[i],False)

		struct.geometry = numpy.concatenate((struct.geometry,
						     new_mole.geometry))
	
	return struct
	

def closeness_check (struct,inst,napm):
	'''
	Calls the respective functions in structure_handling to conduct the closeness check on the structure
	The existence of the closeness criteria (e.g., COM_dist) in the structure_generation section enables the calling
	'''
	sname = "structure_generation"
	original_struct = struct
	struct = structure_handling.cell_modification(struct,len(struct.geometry)/napm,napm)
	COM_dist, atom_dist, inverse_dist, specific_radius_proportion, specific_radius_std, specific_radius_lower_bound = inst.get_keywords([[sname,"COM_dist",None],[sname,"atom_dist",None],[sname,"inverse_dist",None],[sname,"specific_radius_proportion",None],[sname,"specific_radius_std",0],[sname,"specific_radius_lower_bound",None]],True)

	pair_dist = inst.get_with_default(sname,"pair_dist",None)

	if COM_dist!=None:
		result = structure_handling.COM_distance_check(struct,nmpc=len(struct.geometry)/napm,lowerbound=COM_dist)
		if not result:
#			print "Failed COM test"
			return False
	if atom_dist!=None:
		result = structure_handling.atom_distance_check_1(struct,nmpc=len(struct.geometry)/napm,lowerbound=atom_dist)
			
		if not result:
#			print "Failed atom_dist test"
			if inst.get_info_level()>=3:
				write_log.write_master_log(inst,"Structure failed intermolecular atomic distance check")

			return False

	if specific_radius_proportion!=None:
		sr = misc.half_gaussian_sampling_lower(specific_radius_proportion, specific_radius_std, specific_radius_lower_bound)
		#This allows a fuzzy sr to be used, according to a Gaussian distribution. If specific_radius_std=0, then only samples the center

		
		cr = inst.get_with_default(sname,"sr_custom_radius",{},eval=True)
		crp = inst.get_with_default(sname,"sr_custom_radius_pairs",{},eval=True)
			
		result = structure_handling.specific_radius_check(struct,
								  nmpc=len(struct.geometry)/napm,
								  proportion=sr,
								  custom_radius=cr,
								  custom_radius_pairs=crp)

		if not result:
			if inst.get_info_level()>=3:
				write_log.write_master_log(inst,"Structure failed specific radius check; sr="+str(sr))
			return False

		try:
			original_struct.properties["sr"].append(result)
		except:
			original_struct.properties["sr"] = [result]


	return True

def lattice_vector_generation(struct,inst):
    '''
        Reads in the Bravais type of the lattice
        And create a set of Bravais lattice forming a cell with the constrained volume
        Calls diacheck to perform the diagonal rule test on the new cell
        btype= 1-->Triclinic, 2-->monoclinic, 3-->orthorhombic,
        4-->tetragonal, 5-->Cubic
	'''
    sname = "structure_generation"
    check=inst.check_keywords([["structure_generation","unit_cell_volume"]])
    if not (check==True):
        raise ValueError("Missing key in Instruct: "+check)

	#Retrive information from inst
    standard_volume,vrange,vstd,vupper,p_tolerance,angle_range=\
    inst.get_keywords_single_section("structure_generation",["unit_cell_volume",["volume_ratio_range",[1,1]],["unit_cell_volume_std",0],["unit_cell_volume_upper_bound",None],["p_tolerance",0.25],["angle_range",[30,150]],],True)
    #dcheck = inst.get_boolean("structure_generation","diagonal_rule_enabled")
    btype,sgp,ax_variance,by_variance,cz_variance = \
	inst.get_keywords_single_section("structure_generation",[["bravais_system",-1],["space_group",-1],["ax_variance",[0.5,2]],["by_variance",[0.5,2]],["cz_variance",[0.5,2]]],True)

    if btype == -1:
        btype = int(random.uniform(1,6))

    if sgp!=-1 and sgp!=0:
        sg=sgroup.Sgroup(sgp)
        btype = sg.blt


    sequence_interp={0:[0,1,2],1:[0,2,1],2:[1,0,2],3:[1,2,0],4:[2,0,1],5:[2,1,0]} #In order for equivalence between the three lattice vectors, a number is randomized, which is then interpreted by this value into actual 

    struct.properties['bravais_system']=btype
#	total_volume=standard_volume*random.uniform(vrange[0],vrange[1])
    total_volume = misc.half_gaussian_sampling_upper(standard_volume,vstd,upperbound=vupper)
    cube_root = total_volume**(1.0/3)
    angles=["alpha","beta","gamma"]
    lengths=["a","b","c"]
    vectors=["lattice_vector_a","lattice_vector_b","lattice_vector_c"]
    for i in range (3):
        struct.properties[vectors[i]]=[0,0,0]
    principle_variance=[ax_variance,by_variance,cz_variance]
	#The above three values allows streamlined interacting with struct.properties
	
    counter=0
    while True:
        counter+=1
        if counter==1000:
            raise RuntimeError("Repeated failure to generate Bravais lattice; check variance parameters")
        struct.properties['alpha']=struct.properties['beta']=struct.properties['gamma']=90
        if btype==1:
            for j in range (3):	
                struct.properties[angles[j]]=random.uniform(angle_range[0],angle_range[1])
            if struct.properties["alpha"]+struct.properties["beta"]+struct.properties["gamma"]>355: #Cannot be flat on the plane
                continue
            if struct.properties["alpha"]>struct.properties["beta"]+struct.properties["gamma"]-5 or \
			   struct.properties["gamma"]>struct.properties["beta"]+struct.properties["alpha"]-5 or \
			   struct.properties["beta"]>struct.properties["gamma"]+struct.properties["alpha"]-5: 
			   #This will force the structure to be a plane
                continue

                p=calc_p(struct.properties["alpha"],struct.properties["beta"],struct.properties["gamma"])
                if p<p_tolerance:
                    continue
        elif btype==2:
            struct.properties[angles[1]]=random.uniform(angle_range[0],angle_range[1]) #Forcing beta to be the non-90 angle
            p=palc_p(struct.properties["alpha"],struct.properties["beta"],struct.properties["gamma"])
            if p<p_tolerance:
                continue


        if btype==1 or btype==2 or btype==3: #Triclinic, Monoclinic and Orthorhombic cells have a!=b!=c
            sequence = int(random.uniform(0,6))
            v1, v2, v3 = sequence_interp[sequence]
            struct.properties[vectors[v1]][v1] = random.uniform(principle_variance[v1][0],principle_variance[v1][1])*cube_root
            vol_remain = total_volume / struct.properties[vectors[v1]][v1]
            lower = max(principle_variance[v2][0]*cube_root,vol_remain/(cube_root*principle_variance[v3][1]))
            #Has to leave a certain range for v3 to fall into the required range
            upper = min(principle_variance[v2][1]*cube_root,vol_remain/(cube_root*principle_variance[v3][0]))
            if lower>upper:
                continue
            struct.properties[vectors[v2]][v2] = random.uniform(lower,upper)
            struct.properties[vectors[v3]][v3] = vol_remain / struct.properties[vectors[v2]][v2]

        elif btype==4: #Tetragonal cell has a=a!=c
            lower = max(cube_root*ax_variance[0],cube_root*by_variance[0],(total_volume/(cube_root*cz_variance[1]))**0.5)
            upper = min(cube_root*ax_variance[1],cube_root*by_variance[1],(total_volume/(cube_root*cz_variance[0]))**0.5)
            if lower>upper:
                raise RuntimeError("Principle component variance unsuitable for tetragonal structure generation")
            struct.properties["lattice_vector_a"][0]=struct.properties["lattice_vector_b"][1]=random.uniform(lower,upper)
            struct.properties["lattice_vector_c"][2]=total_volume/(struct.properties["lattice_vector_a"][0]**2)

        elif btype==5: #Cubic structure the same for all
            for i in range (3):
                struct.properties[vectors[i]][i] = cube_root
        else:
            raise ValueError("Unsupported Bravais Lattice Type number")
        alpha,beta,gamma = numpy.deg2rad(struct.properties["alpha"]), numpy.deg2rad(struct.properties["beta"]), numpy.deg2rad(struct.properties["gamma"])

        struct.properties["a"]=struct.properties["lattice_vector_a"][0]
        struct.properties["b"]=struct.properties["lattice_vector_b"][1]/numpy.sin(gamma)
        struct.properties["lattice_vector_b"][0] = struct.properties["b"]*numpy.cos(gamma)
        struct.properties["c"]=struct.properties["lattice_vector_c"][2]*struct.properties["lattice_vector_b"][1]/(numpy.sin(beta)**2*struct.properties["lattice_vector_b"][1]**2-(struct.properties["b"]*numpy.cos(alpha)-struct.properties["lattice_vector_b"][0]*numpy.cos(beta))**2)**0.5
        struct.properties["lattice_vector_c"][0]=struct.properties["c"]*numpy.cos(beta)
        struct.properties["lattice_vector_c"][1]=(struct.properties["c"]*struct.properties["b"]*numpy.cos(alpha)-struct.properties["lattice_vector_c"][0]*struct.properties["lattice_vector_b"][0])/struct.properties["lattice_vector_b"][1]
        struct.properties["unit_cell_volume"] = numpy.linalg.det([struct.properties["lattice_vector_a"],struct.properties["lattice_vector_b"],struct.properties["lattice_vector_c"]])


#		if not dcheck or structure_handling.cell_diagonal_rule(struct)==True:
			#If the user wishes the cell to be automatically orthogonal
			#The cell_diagonal_rule check should be enabled
			#WARNING: This will disable exploring the alternative space group settings
        break
    return struct
 	
"""
def lattice_vector_generation_obselete(struct,inst):
        '''
        Reads in the Bravais type of the lattice
        And create a set of Bravais lattice forming a cell with the constrained volume
        Calls diacheck to perform the diagonal rule test on the new cell
        btype= 1-->Triclinic, 2-->monoclinic, 3-->orthorhombic,
        4-->tetragonal, 5-->Cubic
	'''
	check=inst.check_keywords([["structure_generation","unit_cell_volume"],["structure_generation","p_tolerance"],["structure_generation","bravais_system"]])
	if not (check==True):
		raise ValueError("Missing key in Instruct: "+check)

	#Retrive information from inst
	if int(inst.get("structure_generation","bravais_system"))==-1:
		btype=int(random.uniform(1,6))
	else:
		btype=int(inst.get("structure_generation","bravais_system"))
	standard_volume,vrange,vstd,vupper,p_tolerance,lattice_variance,angle_range,dcheck=\
	inst.get_keywords_single_section("structure_generation",["unit_cell_volume",["volume_ratio_range",[1,1]],["unit_cell_volume_std",0],["unit_cell_volume_upper_bound",None],"p_tolerance","lattice_variance","angle_range",["diagonal_rule_enabled",True]],True)

#	print "This is btype", btype

	sequence_interp={0:[0,1,2],1:[0,2,1],2:[1,0,2],3:[1,2,0],4:[2,0,1],5:[2,1,0]} #In order for equivalence between the three lattice vectors, a number is randomized, which is then interpreted by this value into actual 

        struct.properties['bravais_system']=btype
#	total_volume=standard_volume*random.uniform(vrange[0],vrange[1])
	total_volume = misc.half_gaussian_sampling_upper(standard_volume,vstd,upperbound=vupper)
	print("In lattice generation, this is standard_volume = %f,vstd = %f ,upperbound = %f" % (standard_volume,vstd,upperbound))
	angles=["alpha","beta","gamma"]
	lengths=["a","b","c"]
	vectors=["lattice_vector_a","lattice_vector_b","lattice_vector_c"]
	#The above three values allows streamlined interacting with struct.properties
        if btype==1:
            while True:
		for j in range (3):
			struct.properties[angles[j]]=random.uniform(angle_range[0],angle_range[1])
                p=calc_p(struct.properties["alpha"],struct.properties["beta"],struct.properties["gamma"])
                if p<p_tolerance:
                    continue
                abc=total_volume/p
                sequence=int(random.uniform(0,6))
		#Generate a random sequence for the three lattice vector lengths to be calculated
                for i in range (2):
			struct.properties[lengths[sequence_interp[sequence][i]]]=abc**(1/(3-i+0.0))*random.uniform(lattice_variance[0],lattice_variance[1])
                        abc/=struct.properties[lengths[sequence_interp[sequence][i]]]
                if abc<(total_volume/p)**(1/(3+0.0))*lattice_variance[0]:
		#The last vector length is too small
			continue
		struct.properties[lengths[sequence_interp[sequence][2]]]=abc
                structure_handling.set_lattice_vectors(struct)
		if not dcheck or structure_handling.cell_diagonal_rule(struct)==True:
			break
        if btype==2: #Always requires beta to be non-90
            while True:
#                picked=int(random.uniform(0,3))
                struct.properties['alpha']=struct.properties['beta']=struct.properties['gamma']=90
		struct.properties[angles[1]]=random.uniform(angle_range[0],angle_range[1]) #Forcing beta to be the non-90 angle
                p=calc_p(struct.properties['alpha'],struct.properties['beta'],struct.properties['gamma'])
                if p<p_tolerance:
                    continue
                abc=total_volume/p
                sequence=int(random.uniform(0,6))
                for i in range (2):
			struct.properties[lengths[sequence_interp[sequence][i]]]=abc**(1/(3-i+0.0))*random.uniform(lattice_variance[0],lattice_variance[1])
                        abc/=struct.properties[lengths[sequence_interp[sequence][i]]]
                if abc<(total_volume/p)**(1/(3+0.0))*lattice_variance[0]:
                    continue
		struct.properties[lengths[sequence_interp[sequence][2]]]=abc
                structure_handling.set_lattice_vectors(struct)
                if not dcheck or structure_handling.cell_diagonal_rule(struct)==True:
                    break
        if btype==3:
            while True:
		struct.properties['alpha']=struct.properties['beta']=struct.properties['gamma']=90
		abc=total_volume
                sequence=int(random.uniform(0,6))
                for i in range (2):
                        struct.properties[lengths[sequence_interp[sequence][i]]]=abc**(1/(3-i+0.0))*random.uniform(lattice_variance[0],lattice_variance[1])
                        abc/=struct.properties[lengths[sequence_interp[sequence][i]]]
                if abc<(total_volume)**(1/(3+0.0))*lattice_variance[0]:
                    continue
                struct.properties[lengths[sequence_interp[sequence][2]]]=abc
                structure_handling.set_lattice_vectors(struct)
		break
        if btype==4: #Tetragonal
		struct.properties['alpha']=struct.properties['beta']=struct.properties['gamma']=90
		abc=total_volume
		picked=int(random.uniform(0,3))
		struct.properties[lengths[picked]]=abc**(1/(3-0.0))*random.uniform(lattice_variance[0],lattice_variance[1])
		for j in range (3):
			if j!=picked:
				struct.properties[lengths[j]]=(abc/struct.properties[lengths[picked]])**0.5
		structure_handling.set_lattice_vectors(struct)
	if btype==5: #Cubic
		struct.properties['alpha']=struct.properties['beta']=struct.properties['gamma']=90
		struct.properties[lengths[0]]=struct.properties[lengths[1]]=struct.properties[lengths[2]]=total_volume**(1/(3+0.0))
		structure_handling.set_lattice_vectors(struct)
"""
def calc_p (alpha,beta,gamma):
	'''
	Calculates the adjustment to the cell volume posed by the angles
	volume=p*a*b*c
	'''
	ca=numpy.cos(numpy.deg2rad(alpha))
	cb=numpy.cos(numpy.deg2rad(beta))
	cc=numpy.cos(numpy.deg2rad(gamma))
	return (1-ca*ca-cb*cb-cc*cc+2*ca*cb*cc)**0.5

