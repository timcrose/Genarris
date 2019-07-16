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
from core import user_input, structure, structure_handling
import random
import numpy



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

def lattice_vector_generation(struct,btype):
        '''
        Reads in the Bravais type of the lattice
        And create a set of Bravais lattice forming a cell with the constrained volume
        Calls diacheck to perform the diagonal rule test on the new cell
        btype= 1-->Triclinic, 2-->monoclinic, 3-->orthorhombic,
        4-->tetragonal, 5-->Cubic
	'''
        if btype==-1:
            btype=int(random.uniform(1,6))
	sequence_interp={0:[0,1,2],1:[0,2,1],2:[1,0,2],3:[1,2,0],4:[2,0,1],5:[2,1,0]} #In order for equivalence between the three lattice vectors, a number is randomized, which is then interpreted by this value into actual 

        struct.properties['bravais_system']=btype
	ui=user_input.get_config()
	standard_volume=ui.get_eval("unit_cell_settings","volume")
	range=ui.get_eval("unit_cell_settings","volume_ratio_range")
	p_tolerance=ui.get_eval("unit_cell_settings","p_tolerance")
	lattice_variance=ui.get_eval("unit_cell_settings","lattice_variance")
	total_volume=standard_volume*random.uniform(range[0],range[1])
	angles=["alpha","beta","gamma"]
	lengths=["a","b","c"]
	vectors=["lattice_vector_a","lattice_vector_b","lattice_vector_c"]
	#The above three values allows streamlined interacting with struct.properties
        if btype==1:
            while True:
		for j in range (3):
			struct.properties[angles[j]]=random.uniform(60,120)
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
		if structure_handling.cell_diagonal_rule(struct)==True:
			break
        if btype==2:
            while True:
                picked=int(random.uniform(0,3))
                struct.properties['alpha']=struct.properties['beta']=struct.properties['gamma']=90
		struct.properties[angles[picked]]=random.uniform(60,120)
                p=calc_p(struct.properties['alpha'],struct.properties['beta'],struct.properties['gamma'])
                if p<p_tolerance:
                    continue
                abc=total_volume/p
                sequence=int(random.uniform(0,6))
                for i in range (2):
			struct.properties[lengths[sequence_interp[sequence][i]]]=abc**(1/(3-i+0.0))*random.uniform(lattice_variance[0],lattice_variance[1])
                        abc/=struct.properties[lengths[sequence_interp[sequence]]]
                if abc<(total_volume/p)**(1/(3+0.0))*lattice_variance[0]:
                    continue
		struct.properties[lengths[sequence_interp[sequence][2]]]=abc
                structure_handling.set_lattice_vectors(struct)
                if structure_handling.cell_diagonal_rule(struct)==True:
                    break
        if btype==3:
            while True:
		struct.properties['alpha']=struct.properties['beta']=struct.properties['gamma']=90
           
		abc=total_volume/p
                sequence=int(random.uniform(0,6))
                for i in range (2):
                        struct.properties[lengths[sequence_interp[sequence][i]]]=abc**(1/(3-i+0.0))*random.uniform(lattice_variance[0],lattice_variance[1])
                        abc/=struct.properties[lengths[sequence_interp[sequence]]]
                if abc<(total_volume/p)**(1/(3+0.0))*lattice_variance[0]:
                    continue
                struct.properties[lengths[sequence_interp[sequence][2]]]=abc
                structure_handling.set_lattice_vectors(struct)
          break
        if btype==4: #Tetragonal
		struct.properties['alpha']=struct.properties['beta']=struct.properties['gamma']=90
		abc=total_volume
		picked=int(random.uniform(0,3))
		struct.properties[lengths[picked]]=abc**(1/(3-0.0))*random.uniform(parameters.lattice_variance[0],parameters.lattice_variance[1])
		for j in range (3):
			if j!=picked:
				struct.properties[lengths[j]]=(abc/struct.properties[lengths[picked]])**0.5
		structure_handling.set_lattice_vectors(struct)
	if btype==5: #Cubic
		struct.properties['alpha']=struct.properties['beta']=struct.properties['gamma']=90
		struct.properties[lengths[0]]=struct.properties[lengths[1]]=struct.properties[lengths[2]]=total_volume**(1/(3+0.0))
		structure_handling.set_lattice_vectors(struct)

def calc_p (alpha,beta,gamma):
	'''
	Calculates the adjustment to the cell volume posed by the angles
	volume=p*a*b*c
	'''
	ca=numpy.cos(numpy.deg2rad(alpha))
	cb=numpy.cos(numpy.deg2rad(beta))
	cc=numpy.cos(numpy.deg2rad(gamma))
	return (1-ca*ca-cb*cb-cc*cc+2*ca*cb*cc)**0.5

