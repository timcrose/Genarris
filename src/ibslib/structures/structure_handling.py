
import numpy
import numpy as np

from ase.data.vdw import vdw_radii
from ase.data import atomic_numbers

from ibslib import Structure


import copy
#import parameters

molar_mass={"H":1,"C":12,"N":14,"O":16,"S":32, 'H':1, 'C':12, 'N':14, 'O':16,'S':32,'Cl':35.5,'Br':80,'F':19}
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
	trans_vec = [target[x]-cm[x] for x in range (3)]
	cell_translation(molec,trans_vec,create_duplicate=False)
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
    for i in range(int(nmpc)):
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

def full_nearest_image(struct, nmpc=2, lowerbound=None, 
                         enable_timings=False):
    '''
    Purpose:
        Computes full nearest image convention for all molecules in the cell
          excluding the molecules own image to itself. 
    lowerbound: (min_dist, min_sr)
    '''
    if enable_timings:
        start_time = time.time()
    total_atoms=len(struct.geometry)
    napm=total_atoms/nmpc
    
    radius_list = []
    for i in range(0,total_atoms):
        radius_list.append(radius.radius[struct.geometry[i][3]])
#    lattice_vector_a = np.array([x if abs(x)>0.000001 else 0 
#                                 for x in struct.properties['lattice_vector_a']])
#    lattice_vector_b = np.array([x if abs(x)>0.000001 else 0 
#                                 for x in struct.properties['lattice_vector_b']])
#    lattice_vector_c = np.array([x if abs(x)>0.000001 else 0 
#                                 for x in struct.properties['lattice_vector_c']])
    lattice_vector_a = np.array(struct.properties['lattice_vector_a'])
    lattice_vector_b = np.array(struct.properties['lattice_vector_b'])
    lattice_vector_c = np.array(struct.properties['lattice_vector_c'])
    lattice_matrix = np.array([lattice_vector_a, 
                               lattice_vector_b,
                               lattice_vector_c]).T
    inverse_lattice_matrix = np.linalg.inv(lattice_matrix)
        
    evaluations = 0
    dist_list = []
    sr_list = []
    # Don't need to compute for the last molecule
    for a1 in range(0, total_atoms-napm):
        a1_radius = radius_list[a1 % napm]
        a1_position = np.array([struct.geometry[a1][0],
                                struct.geometry[a1][1],
                                struct.geometry[a1][2]])
        a1_molecule = int(a1 / napm)
        next_molecule_index = (a1_molecule+1)*napm
        # Compare a1 to all further molecules in upper triange
        for a2 in range(next_molecule_index, total_atoms):
            a2_radius = radius_list[a2 % napm]
            a2_position = np.array([struct.geometry[a2][0],
                                    struct.geometry[a2][1],
                                    struct.geometry[a2][2]])
            distance_vector = a2_position - a1_position
            distance_crystal_basis = np.dot(inverse_lattice_matrix,
                                            distance_vector)

            for i,entry in enumerate(distance_crystal_basis):
                if entry > 0.50:
                    a2_position -= lattice_matrix[:,i]
                elif entry < -0.50:
                    a2_position += lattice_matrix[:,i]
                                
            # Calculate the final distance
            final_dist = np.linalg.norm(a2_position - a1_position)
            dist_list.append(final_dist)
            
            sum_radius = a1_radius+a2_radius
            final_sr = final_dist / sum_radius
            sr_list.append(final_sr)
            
            if lowerbound:
                if lowerbound[0]:
                    if final_dist < lowerbound[0]:
                        return False
                        #return(np.min(dist_list),np.min(sr_list))
                if lowerbound[1]:
                    if final_sr < lowerbound[1]:
                        return False
                        #return(np.min(dist_list),np.min(sr_list))
                
            evaluations += 1
    
    if enable_timings:
        end_time = time.time()
        print('Full nearest image: {} for {}'
              .format(end_time - start_time, evaluations))

    return True

def calc_atom_dist_manny(struct, nmpc=None, lowerbound=None, 
                         enable_timings=False, sr=0.85):
    '''
    Purpose:
        Constructs full supercell and compares intermolecular distances 
          in the supercell using vertorized calculations.
    '''
    if enable_timings:
        start_time = time.time()
        
    total_atoms=len(struct.geometry)
    napm=int(total_atoms/nmpc)
    if total_atoms % napm != 0:
        raise Exception("Number of molecules per cell is not compatable "+
                        "with the number of atoms.")
    
    if lowerbound:
        atom_dist = lowerbound[0]
        sr = lowerbound[1]

    # Array for translations to construct supercell
    tr=[[0,0,0],[0,0,1],[0,0,-1],
        [0,1,0],[0,-1,0],[0,1,1],
        [0,1,-1],[0,-1,1],[0,-1,-1],
        [1,0,0],[-1,0,0],[1,0,1],
        [1,0,-1],[-1,0,1],[-1,0,-1],
        [1,1,0],[1,-1,0],[-1,1,0],
        [-1,-1,0],[1,1,1],[1,1,-1],
        [1,-1,1],[1,-1,-1],[-1,1,1],
        [-1,1,-1],[-1,-1,1],[-1,-1,-1]]
    
    # Prepring array of the atoms in a single molecule
    atom_list = []
    radius_list = []
    for i in range(0,total_atoms):
        radius_list.append(vdw_radii[atomic_numbers[struct.geometry[i][3]]])
#        atom_list.append(struct.geometry[i][3])
        
    
    # Preparing numpy array of the current geometry
    current_geometry = np.array([struct.geometry[0][0],
                                 struct.geometry[0][1],
                                 struct.geometry[0][2]]).reshape(1,3)
    for atom in range(1,len(struct.geometry)):
        coords = np.array([struct.geometry[atom][0], 
                           struct.geometry[atom][1],
                           struct.geometry[atom][2]]).reshape(1,3)
        current_geometry = np.append(current_geometry, coords, axis=0)
#    print(current_geometry)
    # Checking all distances in current geometry
    rows = current_geometry.shape[0]
    x_matrix = np.repeat(current_geometry[:,0].reshape(rows,1),rows, axis=1)
    y_matrix = np.repeat(current_geometry[:,1].reshape(rows,1),rows, axis=1)
    z_matrix = np.repeat(current_geometry[:,2].reshape(rows,1),rows, axis=1)
    x_dist = np.square(x_matrix - x_matrix.T)
    y_dist = np.square(y_matrix - y_matrix.T)
    z_dist = np.square(z_matrix - z_matrix.T)
    # I just changes np.power(,0.5) to np.sqrt and save 30% of time!
    current_dist = np.sqrt(x_dist + y_dist + z_dist)

    # Extracting sr from current geometry
    columns = x_dist.shape[1]
    current_intermolecular_dist = []
    current_atom_dist = []
    current_sr_proportion = []
    for row in range(0,napm*(nmpc-1)):
        molecule = int(row/napm)+1
        atom_row = radius_list[row % napm]
        for column in range(napm*molecule, columns):
            total_dist = current_dist[row,column]
            current_atom_dist.append(total_dist)
            # Calculate sr sum list
            atom_column = radius_list[column % napm]
            sum_radius = atom_row + atom_column
            current_sr_proportion.append(total_dist / sum_radius)
    
    current_min_atom_dist = np.min(current_atom_dist)
    current_min_sr_proportion = np.min(current_sr_proportion)
    
    if enable_timings:
        current_geo_time = time.time()
        print('Time to evaluate current geometry: {}'
              .format(current_geo_time - start_time))

    if current_min_sr_proportion < sr:
        print('Exited on unit cell geometry check with atom dist '
              '{} and sr {}'.format(current_min_atom_dist, current_min_sr_proportion))
        return (current_min_atom_dist, current_min_sr_proportion)

    print('DID NOT exited on unit cell geometry check')
    
    # Preparing lattice vectors for translations to geometry
    lattice_vector_a = np.array([x if abs(x)>0.000001 else 0 
                                 for x in struct.properties['lattice_vector_a']])
    lattice_vector_b = np.array([x if abs(x)>0.000001 else 0 
                                 for x in struct.properties['lattice_vector_b']])
    lattice_vector_c = np.array([x if abs(x)>0.000001 else 0 
                                 for x in struct.properties['lattice_vector_c']])
    lattice_matrix = np.array([lattice_vector_a, 
                               lattice_vector_b,
                               lattice_vector_c]).T
    #print(lattice_matrix)
    # Preparing super cell
    tr.remove([0,0,0])
    super_cell = current_geometry
    for translation in tr:
        final_trans = np.sum(np.dot(lattice_matrix, np.diag(translation)), 
                             axis=1)
        super_cell = np.append(super_cell, current_geometry+final_trans, 
                               axis=0)
    num_supercell_atoms = super_cell.shape[0]
       
    # Until here takes 0.002 seconds, which is <2% of total calculation
    if enable_timings:
        super_cell_time = time.time()
        print('Supercell construction took {}'
              .format(super_cell_time - start_time))    
    
    # This part takes 0.02 seconds which is 20% of entire calculation
    rows = super_cell.shape[0]
    x_matrix = np.repeat(super_cell[:,0].reshape(rows,1),rows, axis=1)
    y_matrix = np.repeat(super_cell[:,1].reshape(rows,1),rows, axis=1)
    z_matrix = np.repeat(super_cell[:,2].reshape(rows,1),rows, axis=1)
    
    if enable_timings:
        mid_time = time.time()
        print('Setup for calculation has completed in {}'
              .format(mid_time - super_cell_time))
    
    x_dist = np.square(x_matrix - x_matrix.T)
    y_dist = np.square(y_matrix - y_matrix.T)
    z_dist = np.square(z_matrix - z_matrix.T)
    # I just changes np.power(,0.5) to np.sqrt and save 30% of time!
    all_dist = np.sqrt(x_dist + y_dist + z_dist)
    
    if enable_timings:
        mid_2_time = time.time()
        print('All major calculations took {}'
              .format(mid_2_time - mid_time))
    
    
    columns = x_dist.shape[1]
    intermolecular_dist = []
    sr_proportion = []
    for row in range(0,num_supercell_atoms-napm):
        molecule = int(row/napm)+1
        atom_row = radius_list[row % napm]
        for column in range(napm*molecule, columns):
            final_dist = all_dist[row,column]
            intermolecular_dist.append(final_dist)
            # Calculate sr sum list
            atom_column = radius_list[column % napm]
            sum_radius = atom_row + atom_column
            final_sr = final_dist / sum_radius
            sr_proportion.append(final_sr)
            
            if lowerbound:
                if lowerbound[0]:
                    if final_dist < lowerbound[0]:
                        return(np.min(intermolecular_dist),
                               np.min(sr_proportion))
                if lowerbound[1]:
                    if final_sr < lowerbound[1]:
                        return(np.min(intermolecular_dist),
                                np.min(sr_proportion))
            
    if enable_timings:
        end_time = time.time()
        print('Extracting results took {}'
              .format(end_time - mid_2_time))
    
        print('Total time took {}'.format(end_time - start_time))
    
    return (min(intermolecular_dist), min(sr_proportion))


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

if __name__ == '__main__':
    print("You are calling structure_handling directly.")
    print("specify the tasks in main()")