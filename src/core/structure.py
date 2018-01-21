'''
Created on Aug 2, 2013

@author: newhouse
'''


import ast
from collections import defaultdict
import json
import math
import numpy
import numpy as np
import os

class Structure(object):
    '''
    An optimized structure with relevant information
    information includes: geometry, energy, stoichiometry, distance array, 
    '''
    
    def __init__(self):
        '''
        Creates a structure from a given geometry and it's associated properties
        Properties are stored in a dictionary
        '''
        # initialize settings and empty geometry
        self.struct_id = None
        self.input_ref = None        
        self.properties = {}
        self.geometry = numpy.zeros(0, dtype=[('x', 'float32'), ('y', 'float32'), ('z', 'float32'), \
                                      ('element', numpy.str_,16), ('spin', 'float32'), ('charge', 'float32'), ('fixed', 'bool')])
            
    # setters
        
    def build_geo_by_atom(self, x, y, z, element, spin=None, charge=None, fixed=False):
        # increase the size of the array
        size = self.geometry.size
        self.geometry.resize(size + 1)
        # assign valuesx[0
        self.geometry[size]['x'] = x
        self.geometry[size]['y'] = y
        self.geometry[size]['z'] = z
        self.geometry[size]['element'] = element
        self.geometry[size]['spin'] = spin 
        # test for non-assigned spin with math.isnan(a[i]['spin'])
        self.geometry[size]['charge'] = charge 
        # test for non-assigned charge with math.isnan(a[i]['charge'])
        self.geometry[size]['fixed'] = fixed 

    def build_geo_by_whole_atom(self, atom):
        # increase the size of the array
        size = self.geometry.size
        self.geometry.resize(size + 1)
        # assign values
        self.geometry[size] = atom

	def set_lattice_vectors(self, vectors):
		if vectors is None or vectors is False: return False
		for vector in vectors: 
			self.add_lattice_vector(vector)
	def add_lattice_vector(self, vector):
		lattice_vector_name = 'lattice_vector_a'
		if 'lattice_vector_a' in self.properties: lattice_vector_name = 'lattice_vector_b'
		if 'lattice_vector_b' in self.properties: lattice_vector_name = 'lattice_vector_c'
		if 'lattice_vector_c' in self.properties: raise Exception  # lattice vectors are full, 
		self.set_property(lattice_vector_name, list(vector))
		
	def build_geo_whole(self, geometry): self.geometry = geometry
	def build_geo_from_atom_file(self, filepath): self.build_geo_whole_atom_format(read_data(filepath))
	def build_geo_from_json_file(self, filepath): self.loads(read_data(filepath))

	def unpack_geometry(self, text): self.geometry = convert_array(text)    
	def build_geo_whole_atom_format(self, atom_string):
		'''
		Takes an "aims format" string 
		Returns the standard atom array. Returns Fasle if bad string format given.
		If initial spin moment is specified, it is stored in geometry array
		'''
		def add_previous_atom(atom):
			try: spin = atom.get('spin')
			except: spin = None
			try: charge = atom.get('charge')
			except: charge = None
			try: fixed = atom.get('fixed')
			except: fixed = False
			self.build_geo_by_atom(atom['x'], atom['y'], atom['z'],
										 atom['element'], spin, charge, fixed)
		lines_iter = iter(atom_string.split('\n'))
		atom = {}
                has_lattice = False 
		while True:
			try: line = lines_iter.next().split()  # read each line
			except: 
				try: line = lines_iter.__next__().split()
				except: 
                                    add_previous_atom(atom); 
                                    if has_lattice: self.get_unit_cell_volume()
                                    return self.geometry
			if len(line) == 0: continue
			if '#' in line[0]: continue  # skip commented lines
			if line[0] == 'lattice_vector': 
				self.add_lattice_vector([float(line[1]), float(line[2]), float(line[3])])
                                has_lattice=True
				continue
			if line[0] == 'atom':
				if not len(atom) == 0: add_previous_atom(atom)
				atom = {}
				atom['x'] = float(line[1])
				atom['y'] = float(line[2])
				atom['z'] = float(line[3])
				atom['element'] = str(line[4])
			# only affects previous atom
			if 'initial_spin' in line[0]: atom['spin'] = float(line[1])
			if 'initial_charge' in line[0]: atom['charge'] = float(line[1]) 
			if any('constrain_relaxation' in s for s in line) and any('true' in s for s in line): 
				atom['fixed'] = True

                cell_vol = self.get_unit_cell_volume()
                self.properties["cell_vol"] = cell_vol
                print "This is cell_vol: " , cell_vol
        
    def set_input_ref(self, input_ref): self.input_ref = input_ref    
    def set_property(self, key, value):
        try: self.properties[key] = ast.literal_eval(value)
        except: self.properties[key] = value
    
    # getters
    def get_geometry(self): return self.geometry
    def pack_geometry(self): return adapt_array(self.geometry)
    def get_n_atoms(self): return self.geometry.size
    def get_input_ref(self): return  self.input_ref
    def get_struct_id(self): return self.struct_id
    def get_stoic(self): return  calc_stoic(self.geometry)
    def get_stoic_str(self): return self.get_stoic().get_string()
    def get_path(self): return self.get_stoic_str() + '/' + str(self.get_input_ref()) + '/' +str(self.get_struct_id()) 
    def get_property(self, key):
        try: return self.properties.get(key)
        except:
            try: self.reload_structure()  # may not have properly read property
            except Exception as e: 
                try:
                    print(str(e))
                except:
                    pass #For version problem
                return None
    def get_lattice_vectors(self):
        if 'lattice_vector_a' not in self.properties: return False
        return_list = []
        return_list.append(self.get_property('lattice_vector_a'))
        return_list.append(self.get_property('lattice_vector_b'))
        return_list.append(self.get_property('lattice_vector_c'))
        return return_list
    
    def get_geometry_atom_format(self): 
        '''
        Takes a numpy.ndarry with standard "geometry" format.
        Returns a string with structure in standard aims format.
        If atom's spin is spedcified, it's value is located on the line below the atom's coordinates.
        similarly with charge and relaxation constraint.
        '''
        lattice_vectors = self.get_lattice_vectors()
        atom_string = ''
        if lattice_vectors is not False:
            for vector in lattice_vectors:
                atom_string += 'lattice_vector ' + ' '.join(map(str, vector)) + '\n'        
        for item in self.geometry:
            atom_string += 'atom ' + "%.5f" % item['x'] + ' ' + "%.5f" % item['y'] + ' ' + "%.5f" % item['z'] + ' ' + str(item['element']) + '\n'
            if not math.isnan(item['spin']): atom_string += 'initial_moment ' + "%.5f" % item['spin'] + '\n'
            if not math.isnan(item['charge']): atom_string += 'initial_charge ' + "%.5f" % item['charge'] + '\n'
            if item['fixed'] == True: atom_string += 'constrain_relaxation    .true.\n'
        return atom_string
        
    def get_pymatgen_structure(self):
        '''
        Inputs: self
        Outputs: A pymatgen core structure object with basic geometric properties
        '''
        frac_data = self.get_frac_data()
        coords = frac_data[0] # frac coordinates
        atoms = frac_data[1] # site labels
        lattice = LatticeP.from_parameters(a=frac_data[2],
                           b=frac_data[3],
                           c=frac_data[4], 
                           alpha=frac_data[5],
                           beta=frac_data[6],
                           gamma=frac_data[7])
        structp = StructureP(lattice, atoms, coords)
        return structp  


    def get_frac_data(self):
        '''
        Inputs: A np.ndarry structure with standard "geometry" format
        Outputs:  Fractional coordinate data in the form of positions (list), atom_types (list), lattice vector a magnitude, lattice vector b magnitude, lattice vector c magnitude, alpha beta, gamma.
        '''
        geo = self.geometry
        A = self.get_property('lattice_vector_a')
        B = self.get_property('lattice_vector_b')
        C = self.get_property('lattice_vector_c')
        alpha, beta, gamma = self.get_lattice_angles()
        a, b, c = self. get_lattice_magnitudes()
        atoms = [i for i in range(len(geo))]
        lattice_vector = np.transpose([A,B,C])
        latinv = np.linalg.inv(lattice_vector)
        coords = []
        for i in range(len(geo)):
            atoms[i] = geo[i][3]
            coords.append(np.dot(latinv,[geo[i][0],geo[i][1],geo[i][2]]))
        return coords, atoms, a, b, c, alpha, beta, gamma

    def get_lattice_angles(self):
        A = self.get_property('lattice_vector_a')
        B = self.get_property('lattice_vector_b')
        C = self.get_property('lattice_vector_c')
        a = np.linalg.norm(A)
        b = np.linalg.norm(B)
        c = np.linalg.norm(C)
        alpha = self.angle(B,C)
        beta = self.angle(A,C)
        gamma = self.angle(A,B)
        return alpha, beta, gamma

    def get_lattice_magnitudes(self):
        A = self.get_property('lattice_vector_a')
        B = self.get_property('lattice_vector_b')
        C = self.get_property('lattice_vector_c')
        a = np.linalg.norm(A)
        b = np.linalg.norm(B)
        c = np.linalg.norm(C)
        self.set_property("a", a)
        self.set_property("b", b)
        self.set_property("c", c)
        return a, b, c

    def get_unit_cell_volume(self):
        if "cell_vol" in self.properties:
            return self.properties["cell_vol"]
        if "unit_cell_volume" in self.properties:
            return self.properties["unit_cell_volume"]
        A = self.get_property('lattice_vector_a')
        B = self.get_property('lattice_vector_b')
        C = self.get_property('lattice_vector_c')
        self.properties["unit_cell_volume"] = np.linalg.det([A,B,C])
        self.properties["cell_vol"] = np.linalg.det([A,B,C])
        return self.properties["unit_cell_volume"]
        

    def get_atom_distance(self,a1,a2):
        return np.linalg.norm([self.geometry[a1][k]-self.geometry[a2][k] for k in range(3)])

    def angle(self, v1, v2):
        numdot = np.dot(v1,v2)
        anglerad = np.arccos(numdot/(np.linalg.norm(v1)*np.linalg.norm(v2)))
        angledeg = anglerad*180/np.pi
        return angledeg
    
    # json data handling packing
    def dumps(self):
        if "lattice_vector_a" in self.properties:
            self.properties["lattice_vector_a"] == \
                    list(self.properties["lattice_vector_a"])               
        if "lattice_vector_b" in self.properties:
            self.properties["lattice_vector_b"] == \
                    list(self.properties["lattice_vector_b"])
        if "lattice_vector_c" in self.properties:
            self.properties["lattice_vector_c"] == \
                    list(self.properties["lattice_vector_c"])
 
        data_dictionary = {}
        data_dictionary['properties'] = self.properties
        data_dictionary['struct_id'] = self.struct_id
        data_dictionary['input_ref'] = self.input_ref
        data_dictionary['geometry'] = self.geometry.tolist()
        return json.dumps(data_dictionary, indent=4)
        
    def loads(self, json_string):
        data_dictionary = json.loads(json_string)
        self.properties = data_dictionary['properties']
        try:
            self.struct_id = data_dictionary['struct_id']
        except: pass
        try:
            self.input_ref = data_dictionary['input_ref']
        except: pass # if input reference from initial pool then skip this part
        self.build_geo_whole(convert_array(data_dictionary['geometry']))
        
    def file_loads(self,fname):
        f=open(fname,"r")
        st=f.read()
        data_dictionary=json.loads(st)
        self.properties=data_dictionary['properties']
        self.build_geo_whole(convert_array(data_dictionary['geometry']))       

    def loads_try_both(self,fname):
        f=open(fname,"r")
        st = f.read()
        try:
            self.loads(st)
        except:
            try:
                self.build_geo_whole_atom_format(st)
            except:
                return False
        return True 
        
    def update_local_data(self):
        """
        updates local data with that from the shared filesystem
        """
        struct_path = os.path.join(structure_dir, self.get_stoic_str(), str(self.input_ref))
        self.loads(read_data(os.path.join(struct_path, str(self.struct_id)), 'data.json'))
        
    def update_shared_data(self):
        """
        updates local data with that from the shared filesystem
        """
        struct_path = os.path.join(structure_dir, self.get_stoic_str(), str(self.input_ref))
        write_data(os.path.join(struct_path, str(self.struct_id)), 'data.json', self.dumps())                         
        
class StoicDict(defaultdict):    
    
    def __hash__(self):
        return str(self).__hash__()    
    
    def get_string(self):
        keys = self.keys()
        keys.sort()
        stoic_string = ''
        for item in keys:
            stoic_string += str(item) + ':' + str(self[item]) + '_'
        stoic_string = stoic_string[:-1]  # remove last underscore
        return stoic_string
    
def calc_stoic(geo):
    '''
    returns a dictionary representing the stoichiometries
    '''
    stoic = StoicDict(int)
    for item in geo:
        stoic[item['element']] += 1
    return stoic

def get_geo_from_file(file_name):
    ''' 
    given the path to a geometry-style file, returns the geometry in proper format
    '''
    tmp_struct = Structure()
    atom_file = open(file_name, 'r')
    geo = tmp_struct.build_geo_whole_atom_format(atom_file.read())
    atom_file.close()
    return geo

def adapt_array(arr):
    return json.dumps(arr.tolist())

def convert_array(list_of_list):
    """
    takes the array stored in json format and return it to a numpy array with proper dtype
    """
    geometry = numpy.zeros(len(list_of_list), dtype=[('x', 'float32'), ('y', 'float32'), ('z', 'float32'), \
                                   ('element', numpy.str_,16), ('spin', 'float32'), ('charge', 'float32'), ('fixed', 'bool')])
    for i in range(len(list_of_list)): 
        geometry[i]['x'] = list_of_list[i][0]
        geometry[i]['y'] = list_of_list[i][1]
        geometry[i]['z'] = list_of_list[i][2]
        geometry[i]['element'] = str(list_of_list[i][3])
        try:
            geometry[i]['spin'] = list_of_list[i][4]
        except: geometry[i]['spin'] = None
        try:
            geometry[i]['charge'] = list_of_list[i][5]
        except: geometry[i]['charge'] = None
        try:
            geometry[i]['fixed'] = list_of_list[i][6]
        except: geometry[i]['fixed'] = None
    return geometry

def read_data(filepath, filename=None):
    if filename is not None: full_filepath = os.path.join(filepath, filename)
    else: full_filepath = filepath
    d_file = open(full_filepath, 'r')
    contents_string = d_file.read()
    d_file.close()
    return contents_string
