"""
If any part of this module is used for a publication please cite:

X. Li, F. Curtis, T. Rose, C. Schober, A. Vazquez-Mayagoitia, K. Reuter,
H. Oberhofer, and N. Marom "Genarris: Random Generation of Molecular Crystal 
Structures and Fast Screening with a Harris Approximation, ",
J. Chem. Phys., DOI: 10.1063/1.5014038; arXiv 1803.02145 (2018)
"""
'''
This module include functions that does analysis on single structure

created by Patrick Kilecdi on 06/17/2016
'''


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

def get_properties(inst,struct):
	'''
	Retrieves and returns the desired properties of the structures
	'''
	sname = "pool_single_structure_analysis"
	properties = inst.get(sname,"properties_to_get")
	try:
		properties = eval(properties)
	except:
		pass
	if isinstance(properties,str):
		properties = [properties]
	
	result = []
	for prop in properties:
		if prop in struct.properties:
			result.append(struct.properties[prop])
		else:
			result.append(None)
	
	return struct, tuple(result)

def get_property_ratio(inst,struct):
	'''
	Calculates the ratio of properties
	'''
	sname = "pool_single_structure_analysis"
	list_of_pairs = inst.get_eval(sname,"pairs_to_get_ratio")
	if not isinstance(list_of_pairs[0],list):
		list_of_pairs = [list_of_pairs]
	
	result = []
	for pair in list_of_pairs:
		result.append(struct.properties[pair[0]]/(struct.properties[pair[1]]+0.0))
		if len(pair)==3:
			struct.properties[pair[2]] = result[-1]
	
	return struct, tuple(result)
	
def get_property_power(inst,struct):
	'''
	Calculates the ratio of properties
	'''
	sname = "pool_single_structure_analysis"
	properties = inst.get_eval(sname,"properties_to_get_power")
	if not isinstance(properties[0],list):
		properties = [properties]
	
	result = []
	for prop in properties:
		result.append(struct.properties[prop[0]]**prop[1])
		if len(prop)==3:
			struct.properties[prop[2]] = result[-1]
	
	return struct, tuple(result)

def get_property_from_list(inst,struct):
	'''
	Calculates the ratio of properties
	'''
	sname = "pool_single_structure_analysis"
	properties = inst.get_eval(sname,"properties_to_get_from_list")
	if not isinstance(properties[0],list):
		properties = [properties]
	
	result = []
	for prop in properties:
		if prop[0] in struct.properties:
			result.append(struct.properties[prop[0]][prop[1]])
		else:
			result.append(None)
		if len(prop)==3:
			struct.properties[prop[2]] = result[-1]
	
	return struct, tuple(result)

