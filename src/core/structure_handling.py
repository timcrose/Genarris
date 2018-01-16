# -*- coding: utf-8 -*-
"""
Created on Fri May 29 16:09:38 2015

@author: Patrick Kilecdi

Basic modification of crystal structures
Funtions that read in a struct() and necessary arguments, and return a struct()

***Should not interact with instruct.Instruct whatsoever***

"""
import numpy
import numpy as np
#from core import structure #!!!!!!
import structure

import copy
#import parameters
molar_mass={"H":1,"C":12,"N":14,"O":16,"S":32, b'H':1, b'C':12, b'N':14, b'O':16,b'S':32,b'Cl':35.5,b'Br':80,b'F':19}
lat_interp={0:'lattice_vector_a',1:'lattice_vector_b',2:'lattice_vector_c'}
def cell_translation(struct,trans_vec,create_duplicate=True):
	'''
	Translate the entire structure by trans_vec
	'''
	if create_duplicate:
		struct=copy.deepcopy(struct)
	for i in range (len(struct.geometry)):
		for j in range (3):
			struct.geometry[i][j]+=trans_vec[j]
	return struct

def cell_reflection_z(struct,create_duplicate=True):
	'''
	Flips the cell's z axis 
	Should be sufficient to explore the enantiomers
	'''
	if create_duplicate:
		struct=copy.deepcopy(struct)
	for i in range(len(struct.geometry)):
		struct.geometry[i][2]=-struct.geometry[i][2]
	try: #A structure may not have lattice vectors
		struct.properties["lattice_vector_a"][2]=-struct.properties["lattice_vector_a"][2]
	except:
		pass
	try:
		struct.properties["lattice_vector_b"][2]=-struct.properties["lattice_vector_b"][2]
	except:
		pass
	try:
		struct.properties["lattice_vector_c"][2]=-struct.properties["lattice_vector_c"][2]
	except:
		pass
	return struct
def cell_reflection_x(struct,create_duplicate=True):
	'''
	Flips the cell's x axis 
	Should be sufficient to explore the enantiomers
	'''
	if create_duplicate:
		struct=copy.deepcopy(struct)
	for i in range(len(struct.geometry)):
		struct.geometry[i][0]=-struct.geometry[i][0]
	struct.properties["lattice_vector_a"][0]=-struct.properties["lattice_vector_a"][0]
	struct.properties["lattice_vector_b"][0]=-struct.properties["lattice_vector_b"][0]
	struct.properties["lattice_vector_c"][0]=-struct.properties["lattice_vector_c"][0]
	return struct

def cell_rotation(struct,vec=None,theta_deg=None,theta_rad=None,phi_deg=None,phi_rad=None,origin=[0,0,0],deg=None,rad=None,create_duplicate=True):
	if create_duplicate:
		struct=copy.deepcopy(struct)    
	if (deg==None) and (rad==None):
		return False
	if (vec==None) and (((theta_deg==None) and (theta_rad==None)) or ((phi_deg==None) and (phi_rad==None))):
		return False
	if rad==None:
		rad=numpy.deg2rad(deg)
	if (theta_rad==None) and (theta_deg!=None):
		theta_rad=numpy.deg2rad(theta_deg)
	if (phi_rad==None) and (phi_deg!=None):
		phi_rad=numpy.deg2rad(phi_deg)
	if vec==None:
		vec=[numpy.sin(phi_rad)*numpy.cos(theta_rad),numpy.sin(phi_rad)*numpy.sin(theta_rad),numpy.cos(phi_rad)]
	else:
		l=(vec[0]**2+vec[1]**2+vec[2]**2)**0.5
		for j in range (3):
			vec[j]/=l
	c=numpy.cos(rad); s=numpy.sin(rad)
	x,y,z=vec
	mat=[[x*x*(1-c)+c,x*y*(1-c)-z*s,x*z*(1-c)+y*s],
		 [x*y*(1-c)+z*s,y*y*(1-c)+c,y*z*(1-c)-x*s],
		 [x*z*(1-c)-y*s,y*z*(1-c)+x*s,z*z*(1-c)+c]]
	if origin!=[0,0,0]:
		cell_translation(struct,[-j for j in origin],create_duplicate=False)
	for i in range (len(struct.geometry)):
		oldpos=[0,0,0]
		for j in range (3):
			oldpos[j]=struct.geometry[i][j]
		newpos=numpy.dot(mat,oldpos)
		for j in range(3):
			struct.geometry[i][j]=newpos[j]
	try:
		struct.properties["lattice_vector_a"]=numpy.dot(mat,struct.properties["lattice_vector_a"])
	except:
		pass
	try:
		struct.properties["lattice_vector_b"]=numpy.dot(mat,struct.properties["lattice_vector_b"])
	except:
		pass
	try:
		struct.properties["lattice_vector_c"]=numpy.dot(mat,struct.properties["lattice_vector_c"])
	except:
		pass
	if origin!=[0,0,0]:
		cell_translation(struct,origin,create_duplicate=False)
	return struct

def cell_transform_mat(struct,mat,origin=[0,0,0],create_duplicate=True):
	'''
	Transform a structure through a matrix form. 
	Allows the input of an origin
	'''
	if create_duplicate:
		struct=copy.deepcopy(struct)
	if origin!=[0,0,0]:
		cell_translation(struct,[-j for j in origin],create_duplicate=False)
	for i in range (len(struct.geometry)):
		oldpos=[0,0,0]
		for j in range (3):
			oldpos[j]=struct.geometry[i][j]
		newpos=numpy.dot(mat,oldpos)
		for j in range(3):
			 struct.geometry[i][j]=newpos[j]
	if "lattice_vector_a" in struct.properties:
		struct.properties["lattice_vector_a"]=list(numpy.dot(mat,struct.properties["lattice_vector_a"]))
		struct.properties["lattice_vector_b"]=list(numpy.dot(mat,struct.properties["lattice_vector_b"]))
		struct.properties["lattice_vector_c"]=list(numpy.dot(mat,struct.properties["lattice_vector_c"]))
	if origin!=[0,0,0]:
		cell_translation(struct,origin,create_duplicate=False)
	return struct

def cell_extension(struct,extension=[[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]],create_duplicate=True):
	'''
	Extends the structure by the specified extensions
	'''

	if create_duplicate:
		struct=copy.deepcopy(struct)
	napc = len(struct.geometry)
	for extend in extension:
		#Calculate the trans_vector as a dot product of the extension and traspose of lattice_vector matrix
		trans_vector = numpy.dot(extend,[struct.properties["lattice_vector_a"],struct.properties["lattice_vector_b"],struct.properties["lattice_vector_c"]])

		#Create a new copy of each of the atom in the original geometry
		for i in range (napc):
			atom = copy.deepcopy(struct.geometry[i])
			struct.build_geo_by_atom(atom['x']+trans_vector[0],atom['y']+trans_vector[1],atom['z']+trans_vector[2],atom['element'],atom['spin'],atom['charge'],atom['fixed'])
	return struct

def lattice_lower_triangular(struct):
	'''
	Returns a list of lattice vectors that corresponds to the a, b, c, alpha, beta, gamma as specified by struct
	! In row vector form !
	'''
	lattice=[[0 for i in range (3)] for j in range (3)]
	a=struct.properties["a"]; b=struct.properties["b"]; c=struct.properties["c"]
	alpha=numpy.deg2rad(struct.properties["alpha"])
	beta=numpy.deg2rad(struct.properties["beta"])
	gamma=numpy.deg2rad(struct.properties["gamma"])
	lattice[0][0] = a
	lattice[0][1] = 0; lattice[0][2] = 0
	lattice[1][0] = numpy.cos(gamma)*b 
	lattice[1][1] = numpy.sin(gamma)*b 
	lattice[1][2] = 0 
	lattice[2][0] = numpy.cos(beta)*c 
	lattice[2][1] = (b*c*numpy.cos(alpha) - lattice[1][0]*lattice[2][0])/lattice[1][1] 
	lattice[2][2] = (c**2 - lattice[2][0]**2 - lattice[2][1]**2)**0.5 
	return lattice

def lattice_parameters(struct):
	'''
	Returns a dictionary of the lattice parameters a, b, c, alpha, beta, gamma
	'''
	parameters={}
	parameters["a"] = numpy.linalg.norm(struct.properties["lattice_vector_a"])
	parameters["b"] = numpy.linalg.norm(struct.properties["lattice_vector_b"])
	parameters["c"] = numpy.linalg.norm(struct.properties["lattice_vector_c"])
	parameters["alpha"] = numpy.rad2deg(numpy.arccos(numpy.dot(struct.properties["lattice_vector_b"],struct.properties["lattice_vector_c"])/parameters["b"]/parameters["c"]))
	parameters["beta"] = numpy.rad2deg(numpy.arccos(numpy.dot(struct.properties["lattice_vector_a"],struct.properties["lattice_vector_c"])/parameters["a"]/parameters["c"]))
	parameters["gamma"] = numpy.rad2deg(numpy.arccos(numpy.dot(struct.properties["lattice_vector_a"],struct.properties["lattice_vector_b"])/parameters["a"]/parameters["b"]))
	return parameters

def cell_lower_triangular(struct,create_duplicate=True):
	'''
	Sets the cell back to lower triangular form
	Returns a boolean to indicate whether or not the reset was required
	'''
	if create_duplicate:
		struct=copy.deepcopy(struct)
	if abs(struct.properties["lattice_vector_a"][1])<0.001 and abs(struct.properties["lattice_vector_a"][2])<0.001 and abs(struct.properties["lattice_vector_b"][2])<0.001:
		return struct

	struct.properties.update(lattice_parameters(struct)) #Add in case not calculated
	new_lattice=lattice_lower_triangular(struct)
	old_lattice=struct.get_lattice_vectors()
	rots = numpy.dot(numpy.transpose(new_lattice),numpy.linalg.inv(numpy.transpose(old_lattice)))
	struct=cell_transform_mat(struct,rots,create_duplicate=False)

	return struct

def mole_translation(struct,mn,napm,vec=None,frac=None,create_duplicate=True):
	'''
	Translates molecule #mn by a vector
	Vector expressed as in absolute coordinate (vec) or fractional coordinate (frac)
	'''
	if (vec==None) and (frac==None):
		raise RuntimeError("Please input at least one type of translation vector into structure_handling.mole_translation")
	if create_duplicate:
		struct=copy.deepcopy(struct)
	if vec==None:
		vec=[0,0,0]
		for i in range (3):
			for j in range (3):
				vec[j]+=struct.properties[lat_interp[i]][j]*frac[i]
	for i in range (mn*napm,mn*napm+napm):
		for j in range (3):
			struct.geometry[i][j]+=vec[j]
	return struct

def molecule_COM_move (molec,target=[0,0,0],create_duplicate=True):
	'''
	Moves the molecule's COM to the target position
	'''
	if create_duplicate:
		molec = copy.deepcopy(molec)
	cm = cm_calculation(molec,range(0,len(molec.geometry)))
#	print cm
	trans_vec = [target[x]-cm[x] for x in range (3)]
	cell_translation(molec,trans_vec,create_duplicate=False)
#	cm = cm_calculation(molec,range(0,len(molec.geometry)))
#	print cm
	return molec

def cm_calculation (struct,atom_list):
	'''
	Reads in a list of atom
	Find the center of mass
	'''
	cm=[0,0,0]; tm=0;
	for i in range(len(atom_list)):
#		print(struct.geometry[atom_list[i]]["element"])
		tm+=molar_mass[struct.geometry[atom_list[i]][3]]
		for j in range (3):
			cm[j]+=molar_mass[struct.geometry[atom_list[i]][3]]*struct.geometry[atom_list[i]][j]
	for j in range (3):
		cm[j]/=tm
	return cm
	
def move_molecule_in (struct,nmpc,create_duplicate=True):
	'''
	Translate the molecules by the cell vector such that their center of mass lies within the cell
	'''
	if create_duplicate:
		struct=copy.deepcopy(struct)
	napm = int(len(struct.geometry)/nmpc)
        lattice=[struct.properties["lattice_vector_a"],struct.properties["lattice_vector_b"],struct.properties["lattice_vector_c"]]
	lattice=numpy.transpose(lattice)
	latinv=numpy.linalg.inv(lattice) 
	for i in range (nmpc):
		cm=cm_calculation(struct,range(i*napm,i*napm+napm)) #Calculates the center of mass of the molecule
		frac=numpy.dot(latinv,cm) #Calculates the fractional coordinate of the center of mass
		for j in range (0,3): 
		#Go through the fractional coordinate one at a time to see if any translation is needed
			lat=lat_interp[j]
			vec=struct.properties[lat]
			if (frac[j]<-0.0001):
				kk=int(-frac[j]+1)
				for k in range(i*napm,i*napm+napm):
					for l in range (3):
						struct.geometry[k][l]+=kk*vec[l]
			elif (frac[j]>0.99999):
				kk=int(frac[j]+0.00001)
				for k in range(i*napm,i*napm+napm):
					for l in range (3):
						struct.geometry[k][l]-=kk*vec[l]
	return struct

def angle(l1,l2):
	'''
	Returns the angle as in degree between the two vectors
	'''
	return (numpy.rad2deg(numpy.arccos(numpy.dot(l1,l2)/(numpy.linalg.norm(l1)*numpy.linalg.norm(l2)))))


def COM_distance_check(struct,nmpc=None,lowerbound=None):
	'''
	Calculates and returns the minimal distance between the center of mass of the molecules in the unit cell
	if lowerbound is specified, then returns a True system if the minimal distance is greater than the lowerbound
	'''
	if nmpc==None:
		nmpc=struct.properties["nmpc"]
	napm=int(len(struct.geometry)/nmpc)
	
	cmlist=[cm_calculation(struct,range(napm*i,napm*i+napm)) for i in range (nmpc)] #Calculates the center of mass of all molecules
	tr=[[0,0,0],[0,0,1],[0,0,-1],[0,1,0],[0,-1,0],[0,1,1],[0,1,-1],[0,-1,1],[0,-1,-1],[1,0,0],[-1,0,0],[1,0,1],[1,0,-1],[-1,0,1],[-1,0,-1],[1,1,0],[1,-1,0],[-1,1,0],[-1,-1,0],[1,1,1],[1,1,-1],[1,-1,1],[1,-1,-1],[-1,1,1],[-1,1,-1],[-1,-1,1],[-1,-1,-1]] 
	#Each molecule will have to be compared with the 27 equivalent of the other one in order to conduct the distance check
	min_dist = min([numpy.linalg.norm(struct.properties["lattice_vector_a"]),numpy.linalg.norm(struct.properties["lattice_vector_b"]),numpy.linalg.norm(struct.properties["lattice_vector_c"])])
	if min_dist<lowerbound:
		return False


	for m1 in range (nmpc-1):
		for m2 in range (m1+1,nmpc):
			for tr_choice in range (27):
				new_cm=[cmlist[m2][j]+struct.properties["lattice_vector_a"][j]*tr[tr_choice][0]+struct.properties["lattice_vector_b"][j]*tr[tr_choice][1]+struct.properties["lattice_vector_c"][j]*tr[tr_choice][2] for j in range (3)] #Move the molecule by the fractional vector specified by tr[tr_choice]
				diff=[cmlist[m1][j]-new_cm[j] for j in range (3)]
				if min_dist==None or numpy.linalg.norm(diff)<min_dist:
					min_dist=numpy.linalg.norm(diff)
				if lowerbound!=None and min_dist<lowerbound:
					return False
	if lowerbound==None:
		return min_dist
	elif min_dist>=lowerbound:
		return True
	else:
		return False

def atom_distance_check_1(struct,nmpc=None,lowerbound=None):
	'''
	Calculates and returns the minimal distance between atoms in the unit cell
	***In this version 1, does not consider the distance between atoms from the same molecules***
	if lowerbound is specified, then returns a True system if the minimal distance is greater than the lowerbound
	'''
	total_atoms=len(struct.geometry)
	napm=total_atoms/nmpc
	tr=[[0,0,0],[0,0,1],[0,0,-1],[0,1,0],[0,-1,0],[0,1,1],[0,1,-1],[0,-1,1],[0,-1,-1],[1,0,0],[-1,0,0],[1,0,1],[1,0,-1],[-1,0,1],[-1,0,-1],[1,1,0],[1,-1,0],[-1,1,0],[-1,-1,0],[1,1,1],[1,1,-1],[1,-1,1],[1,-1,-1],[-1,1,1],[-1,1,-1],[-1,-1,1],[-1,-1,-1]]
	min_dist=None

	for a1 in range (total_atoms-napm):
		for tr_choice in tr:
			new_apos = [struct.geometry[a1][j]+struct.properties["lattice_vector_a"][j]*tr_choice[0]+struct.properties["lattice_vector_b"][j]*tr_choice[1]+struct.properties["lattice_vector_c"][j]*tr_choice[2] for j in range (3)]
			if tr_choice == [0,0,0]:
				start = (a1/napm+1)*napm
			else:
				start = (a1/napm)*napm #If the molecule is translated, needs to compare with itself
			for a2 in range (start,total_atoms): #Atoms should not be compared to those from the same molecule
				diff=[struct.geometry[a2][j]-new_apos[j] for j in range (3)]
				if min_dist==None or numpy.linalg.norm(diff)<min_dist:
					min_dist=numpy.linalg.norm(diff)
				if lowerbound!=None and min_dist<lowerbound:
					return False


	if lowerbound==None:
		return min_dist
	elif min_dist>=lowerbound:
		return True
	else:
		return False

	

def specific_radius_check(struct,nmpc=None,proportion=None,custom_radius={},custom_radius_pairs={}):
	'''
	Check the closeness between atoms from different molecules
	The two atoms are not allowed to be within proportion*(r1+r2),
	where r1 and r2 are the radius of the atoms
	'''
	from utilities import radius 
	ra = copy.deepcopy(radius.radius)

	total_atoms=len(struct.geometry)
	napm=total_atoms/nmpc
	min_dist=None

	#First make sure molecules won't run into its own replicas

	min_prop = 100

	tr = [[x,y,z] for x in range (0,2) 
		      for y in range (0,2) 
		      for z in range (0,2) if x!=0 or y!=0 or z!=0]

	lv_mat = [struct.properties["lattice_vector_a"],
		  struct.properties["lattice_vector_b"],
		  struct.properties["lattice_vector_c"]]
	lv_mat = np.transpose(lv_mat)


	for e in custom_radius:
		ra[e] = custom_radius[e]

	for l in custom_radius_pairs:
	#Add them so that the first atom loop will include them
		if l[0] not in ra:
			ra[l[0]] = 0.1
		if l[1] not in ra:
			ra[l[1]] = 0.1

	ra_p = custom_radius_pairs
			

	for a1 in [x for x in range (total_atoms) if struct.geometry[x][3] in ra]:
		for tr_choice in tr:
			old_apos = [struct.geometry[a1][j] for j in range(3)]
			new_apos = numpy.add(old_apos,
					     numpy.dot(lv_mat,tr_choice))

			for a2 in [x for x in range (a1-a1%napm,a1-a1%napm+napm)\
				   if struct.geometry[x][3] in ra
				   or (struct.geometry[x][3],struct.geometry[a1][3]) in ra_p
				   or (struct.geometry[a1][3],struct.geometry[x][3]) in ra_p]:

				diff = [struct.geometry[a2][j]-new_apos[j] for j in range (3)]
				dist = numpy.linalg.norm(diff)
				sum_of_rad = ra[struct.geometry[a1][3]] + ra[struct.geometry[a2][3]]
				if (struct.geometry[a1][3],struct.geometry[a2][3]) in ra_p:
					sum_of_rad = ra_p[(struct.geometry[a1][3],struct.geometry[a2][3])]
				elif (struct.geometry[a2][3],struct.geometry[a1][3]) in ra_p:
					sum_of_rad = ra_p[(struct.geometry[a2][3],struct.geometry[a1][3])]


				if proportion!=None and dist<sum_of_rad*proportion:
					return False
				else:
					min_prop = min(min_prop,dist/sum_of_rad)
					


	#Then compare each atom from one another
	tr = [[x,y,z] for x in range(-1,2)
		      for y in range(-1,2)
		      for z in range(-1,2)]

	struct = move_molecule_in(struct,total_atoms) #Move in to prevent too far away from cell
	for a1 in [x for x in range (total_atoms-napm) if struct.geometry[x][3] in ra]:
		for tr_choice in tr:
			old_apos = [struct.geometry[a1][j] for j in range(3)]
			new_apos = numpy.add(old_apos,
					     numpy.dot(lv_mat,tr_choice))
		
			for a2 in [x for x in range((a1/napm+1)*napm,total_atoms)\
				   if struct.geometry[x][3] in ra
				   or (struct.geometry[x][3],struct.geometry[a1][3]) in ra_p
				   or (struct.geometry[a1][3],struct.geometry[x][3]) in ra_p]:
			#Atoms should not be compared to those from the same molecule
				diff = [struct.geometry[a2][j]-new_apos[j] for j in range (3)]
				dist = numpy.linalg.norm(diff)

				sum_of_rad = ra[struct.geometry[a1][3]] + ra[struct.geometry[a2][3]]
				if (struct.geometry[a1][3],struct.geometry[a2][3]) in ra_p:
					sum_of_rad = ra_p[(struct.geometry[a1][3],struct.geometry[a2][3])]
				elif (struct.geometry[a2][3],struct.geometry[a1][3]) in ra_p:
					sum_of_rad = ra_p[(struct.geometry[a2][3],struct.geometry[a1][3])]

				if proportion!=None and dist < sum_of_rad*proportion:
					return False
				else:
					min_prop = min(min_prop,dist/sum_of_rad)
	return min_prop
		

def cell_diagonal_rule(struct):
	'''
	This test further ensure the orthogonality of the unit cell: not too skewed
	'''
	#The shortest face diagonal is not shorter than the longest edge of the same face
	lengths=['a','b','c']
	vectors=['lattice_vector_a','lattice_vector_b','lattice_vector_c']
	for v1 in range (2):
		for v2 in range (v1+1,3):
			sum_length=numpy.linalg.norm([struct.properties[vectors[v1]][j]+struct.properties[vectors[v2]][j] for j in range (3)])
			diff_length=numpy.linalg.norm([struct.properties[vectors[v1]][j]-struct.properties[vectors[v2]][j] for j in range (3)])
			if sum_length<struct.properties[lengths[v1]] or sum_length<struct.properties[lengths[v2]] or diff_length<struct.properties[lengths[v1]] or diff_length<struct.properties[lengths[v2]]:
				return False
	#The shortest body diagonal is not shorter than the longest edge of the cell
	length=numpy.linalg.norm([struct.properties[vectors[0]][j]+struct.properties[vectors[1]][j]+struct.properties[vectors[2]][j] for j in range (3)]) #a+b+c
	if length<struct.properties['a'] or length<struct.properties['b'] or length<struct.properties['c']:
		return False
	length=numpy.linalg.norm([struct.properties[vectors[0]][j]+struct.properties[vectors[1]][j]-struct.properties[vectors[2]][j] for j in range (3)]) #a+b-c
	if length<struct.properties['a'] or length<struct.properties['b'] or length<struct.properties['c']:
		return False
	length=numpy.linalg.norm([struct.properties[vectors[0]][j]-struct.properties[vectors[1]][j]+struct.properties[vectors[2]][j] for j in range (3)]) #a-b+c
	if length<struct.properties['a'] or length<struct.properties['b'] or length<struct.properties['c']:
		return False
	length=numpy.linalg.norm([-struct.properties[vectors[0]][j]+struct.properties[vectors[1]][j]+struct.properties[vectors[2]][j] for j in range (3)]) #-a+b+c
	if length<struct.properties['a'] or length<struct.properties['b'] or length<struct.properties['c']:
		return False
	return True

def cell_niggli_reduction(struct,napm,create_duplicate=True):

    '''
    Cell modification using Niggli reduction
    '''
    if create_duplicate:
        struct = copy.deepcopy(struct)
    lats = struct.get_lattice_vectors()
    from spglib import niggli_reduce
    reduced_lats =  niggli_reduce(lats)
    if reduced_lats is None:
        return False
    del(struct.properties["lattice_vector_a"])
    del(struct.properties["lattice_vector_b"])
    del(struct.properties["lattice_vector_c"])
    struct.set_lattice_vectors(reduced_lats)
    nmpc = len(struct.geometry)/napm
    cell_lower_triangular(struct,False)
    move_molecule_in(struct,nmpc,False)
    return struct

def cell_modification (struct,nmpc,napm,create_duplicate=True):
	'''
	Method found in the 2011 Lonie paper
	Adjust the cell so that there are no angles beyond the range of 60 and 120
	nmpc = number of molecules per cell
	napm = number of atoms per cell
	'''
	if create_duplicate:
		struct=copy.deepcopy(struct)
	if len(struct.geometry)!=nmpc*napm:
		print("cell_modification check struct.geometry list or nmpc or napm!")
		return False
	test=True
	run=0
	set_lattice_parameters(struct)	
	while test and run<10:
		test=False
		run+=1
		if (struct.properties['gamma']>120) or (struct.properties['gamma']<60):
			test=True
			if struct.properties['a']>=struct.properties['b']:
				frac=numpy.dot(struct.properties["lattice_vector_a"],struct.properties["lattice_vector_b"])/(numpy.linalg.norm(struct.properties["lattice_vector_b"])**2)
				#c=numpy.ceil(abs(frac))*numpy.sign(frac)
				c=numpy.round(frac)
				for j in range (3):
					struct.properties["lattice_vector_a"][j]-=c*struct.properties["lattice_vector_b"][j]
				struct.properties['a']=numpy.linalg.norm(struct.properties['lattice_vector_a'])
			else:
				frac=numpy.dot(struct.properties["lattice_vector_a"],struct.properties["lattice_vector_b"])/(numpy.linalg.norm(struct.properties["lattice_vector_a"])**2)
				#c=numpy.ceil(abs(frac))*numpy.sign(frac)
				c=numpy.round(frac)
				for j in range (3):
					struct.properties["lattice_vector_b"][j]-=c*struct.properties["lattice_vector_a"][j]
				struct.properties['b']=numpy.linalg.norm(struct.properties['lattice_vector_b'])
			struct.properties['alpha']=angle(struct.properties["lattice_vector_b"],struct.properties["lattice_vector_c"])
			struct.properties['beta']=angle(struct.properties["lattice_vector_a"],struct.properties["lattice_vector_c"])
			struct.properties['gamma']=angle(struct.properties["lattice_vector_a"],struct.properties["lattice_vector_b"])


		if (struct.properties['beta']>120) or (struct.properties['beta']<60):
			test=True
			if struct.properties['a']>=struct.properties['c']:
				frac=numpy.dot(struct.properties["lattice_vector_a"],struct.properties["lattice_vector_c"])/(numpy.linalg.norm(struct.properties["lattice_vector_c"])**2)
				#c=numpy.ceil(abs(frac))*numpy.sign(frac)
				c=numpy.round(frac)
				for j in range (3):
					struct.properties["lattice_vector_a"][j]-=c*struct.properties["lattice_vector_c"][j]
				struct.properties['a']=numpy.linalg.norm(struct.properties['lattice_vector_a'])
			else:
				frac=numpy.dot(struct.properties["lattice_vector_a"],struct.properties["lattice_vector_c"])/(numpy.linalg.norm(struct.properties["lattice_vector_a"])**2)
				#c=numpy.ceil(abs(frac))*numpy.sign(frac)
				c=numpy.round(frac)
				for j in range (3):
					struct.properties["lattice_vector_c"][j]-=c*struct.properties["lattice_vector_a"][j]
				struct.properties['c']=numpy.linalg.norm(struct.properties['lattice_vector_c'])
			struct.properties['alpha']=angle(struct.properties["lattice_vector_b"],struct.properties["lattice_vector_c"])
			struct.properties['beta']=angle(struct.properties["lattice_vector_a"],struct.properties["lattice_vector_c"])
			struct.properties['gamma']=angle(struct.properties["lattice_vector_a"],struct.properties["lattice_vector_b"])


		if (struct.properties['alpha']>120) or (struct.properties['alpha']<60):
			test=True
			if struct.properties['b']>=struct.properties['c']:
				frac=numpy.dot(struct.properties["lattice_vector_b"],struct.properties["lattice_vector_c"])/(numpy.linalg.norm(struct.properties["lattice_vector_c"])**2)
				c=numpy.round(frac)                
				for j in range (3):
					struct.properties["lattice_vector_b"][j]-=c*struct.properties["lattice_vector_c"][j]
				struct.properties['b']=numpy.linalg.norm(struct.properties['lattice_vector_b'])
			else:
				frac=numpy.dot(struct.properties["lattice_vector_b"],struct.properties["lattice_vector_c"])/(numpy.linalg.norm(struct.properties["lattice_vector_b"])**2)
				c=numpy.ceil(abs(frac))*numpy.sign(frac)
				for j in range (3):
					struct.properties["lattice_vector_c"][j]-=c*struct.properties["lattice_vector_b"][j]
				struct.properties['c']=numpy.linalg.norm(struct.properties['lattice_vector_c'])
			struct.properties['alpha']=angle(struct.properties["lattice_vector_b"],struct.properties["lattice_vector_c"])
			struct.properties['beta']=angle(struct.properties["lattice_vector_a"],struct.properties["lattice_vector_c"])
			struct.properties['gamma']=angle(struct.properties["lattice_vector_a"],struct.properties["lattice_vector_b"])


	struct=move_molecule_in(struct,nmpc,False)
	return struct
	
def print_aims(struct,fname):
	os=''        
	st="lattice_vector"
	for j in range(0,3):
		st=st+"{0:18.10f}".format(struct.properties["lattice_vector_a"][j])
	os=os+st+"\n"
	st="lattice_vector"
	for j in range(0,3):
		st=st+"{0:18.10f}".format(struct.properties["lattice_vector_b"][j])
	os=os+st+"\n"
	st="lattice_vector"
	for j in range(0,3):
		st=st+"{0:18.10f}".format(struct.properties["lattice_vector_c"][j])
	os=os+st+"\n"
	for i in range (len(struct.geometry)):
		st='atom'
		for j in range(0,3):
			st=st+"%18.10f" % (struct.geometry[i][j])
		os=os+st+' '+struct.geometry[i][3]+"\n"
	os += '\n'
	f=open(fname,"w")
	f.write(os)
	f.close()

def struct_properties_updates(s1,s2):
	'''
	Updates the missing properties in s1 by s2
	'''
	for key in s2.properties:
		if (not key in s1.properties) or (s1.properties[key]==None):
			s1.properties[key]=s2.properties[key]


def set_lattice_vectors(struct):
	'''
	Calculates the coordinates of the lattice vectors according to their lengths and the angles
	'''
	v=["lattice_vector_a","lattice_vector_b","lattice_vector_c"]
	a,b,c=struct.properties["a"],struct.properties["b"],struct.properties["c"]
	alpha,beta,gamma=numpy.deg2rad(struct.properties["alpha"]),numpy.deg2rad(struct.properties["beta"]),numpy.deg2rad(struct.properties["gamma"])
	for item in v:
		struct.properties[item]=[0,0,0]
	struct.properties[v[0]][0]=a
	struct.properties[v[0]][1]=struct.properties[v[0]][2]=0
	struct.properties[v[1]][0]=numpy.cos(gamma)*b
	struct.properties[v[1]][1]=numpy.sin(gamma)*b
	struct.properties[v[1]][2]=0
	struct.properties[v[2]][0]=numpy.cos(beta)*c
	struct.properties[v[2]][1]=b*c*numpy.cos(alpha)-struct.properties[v[1]][0]*struct.properties[v[2]][0]/struct.properties[v[1]][1]
	struct.properties[v[2]][2]=(c**2-struct.properties[v[2]][0]**2-struct.properties[v[2]][1]**2)**0.5

def set_lattice_parameters(struct):
	'''
	Fill in the missing lattice parameters
	a, b, c, alpha, beta, gamma
	'''
	v=["lattice_vector_a","lattice_vector_b","lattice_vector_c"]
	if "a" not in struct.properties:
		struct.properties["a"] = numpy.linalg.norm(struct.properties[v[0]])
	if "b" not in struct.properties:
		struct.properties["b"] = numpy.linalg.norm(struct.properties[v[1]])
	if "c" not in struct.properties:
		struct.properties["c"] = numpy.linalg.norm(struct.properties[v[2]])
	if "alpha" not in struct.properties:
		struct.properties["alpha"] = angle_between(struct.properties[v[1]],struct.properties[v[2]])
	if "beta" not in struct.properties:
		struct.properties["beta"] = angle_between(struct.properties[v[0]],struct.properties[v[2]])
	if "gamma" not in struct.properties:
		struct.properties["gamma"] = angle_between(struct.properties[v[0]],struct.properties[v[1]])

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

def main():
	print("You are calling structure_handling directly.")
	print("specify the tasks in main()")
	struct = structure.Structure()
#	f=open("/lustre/project/nmarom/pool_management_0.2/tmp_0.75/253052e51f/geometry.in","r")
	f=open("/home/xli20/NEW_BTM_TEST/target_13/relaxation/geometry.in.next_step","r")
#	struct.build_geo_whole_atom_format(f.read())
	struct.build_geo_whole_atom_format(f.read())
	f.close()
	molecule_COM_move(struct,create_duplicate=False)
	print struct.get_geometry_atom_format()
	
	

if __name__=='__main__':
	main()
