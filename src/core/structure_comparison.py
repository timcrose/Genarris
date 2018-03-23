"""
If any part of this module is used for a publication please cite:

X. Li, F. Curtis, T. Rose, C. Schober, A. Vazquez-Mayagoitia, K. Reuter,
H. Oberhofer, and N. Marom "Genarris: Random Generation of Molecular Crystal 
Structures and Fast Screening with a Harris Approximation, ",
J. Chem. Phys., DOI: 10.1063/1.5014038; arXiv 1803.02145 (2018)
"""
"""
Created on Fri May 29 13:46:41 2015

@author: Patrick Kilecdi
"""
from core import structure, structure_handling
import numpy
import copy

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


match_molecule_tolerance=0.0001
match_molecule_length_requirement=0.0001
#The above 2 should be kept the same to each other for simple molecules
match_molecule_cross_tolerance=0.001

#ui=user_input.get_config()
#is_duplicate_tolerance=ui.get_eval("comparison_settings","is_duplicate_tolerance")
#energy_tolerance=ui.get_eval("comparison_settings","energy_tolerance")
#verbose=ui.get_eval("run_settings","verbose")
#nmpc=ui.get_eval("unit_cell_settings","num_molecules")

def main(struct,structure_coll,replica):
	'''
	For a structure comparison module, the method main() must be implemented.
	The module will take a new structure and compare it to the given structure collection.
	It will return True if the structure passes a test for uniqueness and for energy
	This module directly compares the geometry of two unit cells and attempts to map one onto another
	If the residual is smaller than is_duplicate_tolerance, than a False will be returned
	'''
	list=structure_coll.get_structures()
	if len(list)==0:
		raise RuntimeError("Suprecollection length is 0 in comparison module")
	if verbose:
		message="Comparison begins: comparing with %i structures in the commmon pool" % (len(list))
		output.local_message(message,replica)	
	minimum_residual=None
	napm=int(len(struct.geometry)/nmpc)
	for i in list:
		current_struct=list[i]
		try:
			if abs(current_struct.properties["energy"]-struct.properties["energy"])>energy_tolerance:
				continue
		except:
			pass
		current_residual=is_duplicate(struct,current_struct,nmpc,napm)
		if minimum_residual==None or minimum_residual>current_residual:
			minimum_residual=current_residual
		if minimum_residual<is_duplicate_tolerance:
			if verbose:
				message="Comparison test finds the new structure to be duplicate. Resi=%f" % (minimum_residual)
				output.local_message(message,replica)
			return False #Meaning the new structures is not acceptable
	if verbose:
		try:
			message="Comparison test finds the new structure to be acceptible! Minimum resi=%f" %(minimum_residual)
		except:
			message="Comparison test finds the new structure to be acceptible! Energy unique from the other structures!"
		output.local_message(message,replica)
	return True
				


def residual_single_molecule (s1,s2,napm,mn1,mn2):
    '''
    Calculates the rough least residual of molecule mn1 of s1 and molecule mn2 of s2
    Everytime picks the closest atom to add to the residual
    Rough because measures are taken to counteract potential reordering of atoms due to TINKER
    '''
    result=0
    picked=[False for j in range (napm)]

    for a1 in range (napm):
        least_length = -1

        for a2 in range (napm):
            if (s1.geometry[mn1*napm+a1][3]==s2.geometry[mn2*napm+a2][3]) and not picked[a2]:
	    #Every time loops through all the unmapped atoms in s2
	    #In search of closest atom of the correct atom type
	    #This is to handle potential misordering of atoms of the same type
                diff=numpy.linalg.norm([s1.geometry[mn1*napm+a1][j]-s2.geometry[mn2*napm+a2][j] for j in range (3)])
                if least_length==-1 or diff<least_length:
                    least_length=diff
                    chosen=a2

        if least_length==-1:
            print("Error in structure_comparison.residual_single_molecule! Atoms not matching")
            return 10000
        result+=least_length**2
        picked[chosen]=True
    return result        

def residual_whole_cell (s1,s2,nmpc,napm,known_pairs=[],create_duplicate=True):
    '''
    Caluclate the whole cell residual of two given unit cells
    Will attempt the map s2's molecule onto s1 using only lattice vector translation
    Optionally reads in a set of known pairs
    Returns a list of: [rough least residual, [all the molecule pairings]]
    '''
    if len(s1.geometry)!=len(s2.geometry):
        raise "s1 and s2 do not have the same amount of atoms"
        return
    if create_duplicate:
        s2=copy.deepcopy(s2)
    result=0
    s2_picked=[]; s1_picked=[]
    for i in range (len(known_pairs)):
        result+=residual_single_molecule(s1,s2,napm,known_pairs[i][0],known_pairs[i][1])
        s2_picked.append(known_pairs[i][1])
        s1_picked.append(known_pairs[i][0])
    lattice=[s2.properties["lattice_vector_a"],s2.properties["lattice_vector_b"],s2.properties["lattice_vector_c"]]
    latinv=numpy.linalg.inv(numpy.transpose(lattice))    
    for m2 in range (nmpc):
        if not (m2 in s2_picked):
            best_resi=None
            for m1 in range (nmpc):
                if not (m1 in s1_picked):
                    cm1=structure_handling.cm_calculation(s1,range(m1*napm,m1*napm+napm))
                    cm2=structure_handling.cm_calculation(s2,range(m2*napm,m2*napm+napm))
                    diff=[cm1[j]-cm2[j] for j in range (3)]
                    frac=numpy.dot(latinv,diff)
                    structure_handling.mole_translation(s2,m2,napm,frac=[round(frac[j]) for j in range (3)],create_duplicate=False)
                    current_resi=residual_single_molecule(s1,s2,napm,m1,m2)
                    if best_resi==None or best_resi>current_resi:
                        best_resi=current_resi
                        chosen=m1
            known_pairs.append([m1,chosen])
            s2_picked.append(m2)
            s1_picked.append(chosen)
            result+=best_resi
    return [result,known_pairs]

def match_molecule(s1,s2,napm,mn,cm1=None,create_duplicate=True):
    '''
    Using rotation, translation, or mirror reflection, 
    tries to match the mn'th molecule of s2 onto s1
    returns the updated structure of s2
    Rotation will not be successfully found if the atoms are not matched in positions
    (ordered in the same way for the two molecules),
    because the vectors to define the rotational axis will not be effectively found
    '''
    if cm1==None:
        cm1=structure_handling.cm_calculation(s1,range(napm))
    if create_duplicate:
        s2=copy.deepcopy(s2)
    cm2=structure_handling.cm_calculation(s2,range(mn*napm,mn*napm+napm))
    trans=[cm1[j]-cm2[j] for j in range (3)]
    structure_handling.cell_translation(s2,trans,False) 
    #Fixes the cm of two molecules to be at the same place
    resi=residual_single_molecule(s1,s2,napm,0,mn)
    if resi<match_molecule_tolerance:
#        print("just a translation!")
        return s2
    lll=0
    while lll<2:
        lll+=1
        vec1=vec2=None
        i=0
        rotation_axis=[0,0,0]
        chosen=None
        while (vec2==None) and (i<napm): 
            #find two vectors of sufficient length that can help find the rotational axis
            diff=[s1.geometry[i][j]-s2.geometry[napm*mn+i][j] for j in range (3)]
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
            if vec2==None:
                vec2=[0,0,0]
            rotation_axis=numpy.cross(vec1,vec2)
        rl=numpy.linalg.norm(rotation_axis)
        if rl>match_molecule_cross_tolerance:
            for j in range (3):
                rotation_axis[j]/=rl
            if chosen==None:
                chosen=0
                while (chosen<napm) and (numpy.linalg.norm(numpy.cross(rotation_axis,[s1.geometry[chosen][j]-cm1[j] for j in range (3)]))<match_molecule_cross_tolerance):
                    chosen+=1
            v1=[s1.geometry[chosen][j]-cm1[j] for j in range (3)]
            v2=[s2.geometry[napm*mn+chosen][j]-cm1[j] for j in range (3)] #center of mass is already alligned
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
            structure_handling.cell_rotation(s2,vec=rotation_axis,origin=cm1,rad=-rad,create_duplicate=False)
            resi=residual_single_molecule(s1,s2,napm,0,mn)
            if resi<match_molecule_tolerance:
#                print("a rotation is invovled. Final resi="+str(resi))
                return s2
#        print("A mirror reflection is involved")
        structure_handling.cell_reflection_z(s2,False)
        trans=[0,0,2*cm1[2]]
        structure_handling.cell_translation(s2,trans,False)
        resi=residual_single_molecule(s1,s2,napm,0,mn)
        if resi<match_molecule_tolerance:
#            print("no further rotation is needed after mirror reflection along the xy plane, resi=", resi)
            return s2        
            

def is_duplicate (s1,s2,nmpc,napm,create_duplicate=True,extended_check=False):
    '''
    Compares two structures when the nmpc and napm are given
    Returns the least total_residual of atom positions
    If extended_check=True, then duplicate structure is automatically created
    The final residual will be divided by 8
    '''
    ex_mul=1
    if extended_check:
        s1=copy.deepcopy(s1)
        s2=copy.deepcopy(s2)
        structure_handling.cell_extension(s1,create_duplicate=False)
        structure_handling.cell_extension(s2,create_duplicate=False)
        ex_mul=8
    elif create_duplicate:
        s2=copy.deepcopy(s2)
    total_residual=None
    cm1=structure_handling.cm_calculation(s1,range(napm))
    for i in range (nmpc):
        match_molecule(s1,s2,napm,i,cm1,create_duplicate=False)
        result=residual_whole_cell(s1,s2,nmpc*ex_mul,napm,[[0,i]],create_duplicate=False)
        if total_residual==None or result[0]<total_residual:
            total_residual=result[0]
        if total_residual<is_duplicate_tolerance*ex_mul:
            break
    return total_residual/ex_mul
            
if __name__=='__main__':
    pass
