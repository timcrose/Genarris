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
from core import structure, structure_handling
import numpy as np
import os


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

def extract_energy_aims (file_name):
	'''
	Extracts energy from the output file of aims calculation
	'''
	aims_out = open(file_name)
	while True:
		line = aims_out.readline()
		if line=="": 
			return False  # energy not converged
		if '  | Total energy corrected        :' in line:
			tokens = line.split()
			energy = float(tokens[5])  # converts from SI string to float
			return energy
	

def extract_geometry_aims (file_name):
	'''
	Extracts and returns a temporary geometry
	'''
	def angle(v1,v2):
		numdot = np.dot(v1,v2)
		anglerad = np.arccos(numdot/(self.leng(v1)*self.leng(v2)))
		angledeg = anglerad*180/np.pi
		return (angledeg)

	aims_out = open(os.path.join(self.working_dir, 'aims.out'))
	while True:
		line = aims_out.readline()
		if line=="": 
			return False  # not converged
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
	tmp_structure=structure.Structure()
	tmp_structure.build_geo_whole_atom_format(atom_string)

	lats= extract_lats_aims()
	latA = [float(lats[0][0]), float(lats[0][1]), float(lats[0][2])]
	latB = [float(lats[1][0]), float(lats[1][1]), float(lats[1][2])]
	latC = [float(lats[2][0]), float(lats[2][1]), float(lats[2][2])]
	tmp_structure.set_property('lattice_vector_a', latA)
	tmp_structure.set_property('lattice_vector_b', latB)
	tmp_structure.set_property('lattice_vector_c', latC)

	latA = np.asarray(latA)
	latB = np.asarray(latB)
	latC = np.asarray(latC)

	temp_vol = np.dot(np.cross(latA, latB), latC)
	alpha = angle(latB, latC)
	beta = angle(latA, latC)
	gamma = angle(latA, latB)
	a = np.linalg.norm(latA)
	b = np.linalg.norm(latB)
	c = np.linalg.norm(latC)
	tmp_structure.set_property('cell_vol', temp_vol)
	tmp_structure.set_property('alpha',alpha)
	tmp_structure.set_property('beta', beta)
	tmp_structure.set_property('gamma', gamma)
	tmp_structure.set_property('a',a)
	tmp_structure.set_property('b', b)
	tmp_structure.set_property('c', c)

	return tmp_structure

def extract_lats_aims (file_name):
	'''
	Extracts the lattice vectors from an aims calculation output
	'''
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
	
	


