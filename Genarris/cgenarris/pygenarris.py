# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.10
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.



from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_pygenarris', [dirname(__file__)])
        except ImportError:
            import _pygenarris
            return _pygenarris
        if fp is not None:
            try:
                _mod = imp.load_module('_pygenarris', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _pygenarris = swig_import_helper()
    del swig_import_helper
else:
    import _pygenarris
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0



def generate_molecular_crystals(*args):
  return _pygenarris.generate_molecular_crystals(*args)
generate_molecular_crystals = _pygenarris.generate_molecular_crystals

def find_allowed_positions_using_molecular_symmetry(*args):
  return _pygenarris.find_allowed_positions_using_molecular_symmetry(*args)
find_allowed_positions_using_molecular_symmetry = _pygenarris.find_allowed_positions_using_molecular_symmetry

def allocate_xtal(*args):
  return _pygenarris.allocate_xtal(*args)
allocate_xtal = _pygenarris.allocate_xtal
class COMPATIBLE_SPG(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, COMPATIBLE_SPG, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, COMPATIBLE_SPG, name)
    __repr__ = _swig_repr
    __swig_setmethods__["spg"] = _pygenarris.COMPATIBLE_SPG_spg_set
    __swig_getmethods__["spg"] = _pygenarris.COMPATIBLE_SPG_spg_get
    if _newclass:spg = _swig_property(_pygenarris.COMPATIBLE_SPG_spg_get, _pygenarris.COMPATIBLE_SPG_spg_set)
    __swig_setmethods__["num_allowed_pos"] = _pygenarris.COMPATIBLE_SPG_num_allowed_pos_set
    __swig_getmethods__["num_allowed_pos"] = _pygenarris.COMPATIBLE_SPG_num_allowed_pos_get
    if _newclass:num_allowed_pos = _swig_property(_pygenarris.COMPATIBLE_SPG_num_allowed_pos_get, _pygenarris.COMPATIBLE_SPG_num_allowed_pos_set)
    __swig_setmethods__["allowed_pos"] = _pygenarris.COMPATIBLE_SPG_allowed_pos_set
    __swig_getmethods__["allowed_pos"] = _pygenarris.COMPATIBLE_SPG_allowed_pos_get
    if _newclass:allowed_pos = _swig_property(_pygenarris.COMPATIBLE_SPG_allowed_pos_get, _pygenarris.COMPATIBLE_SPG_allowed_pos_set)
    __swig_setmethods__["pos_overlap_list"] = _pygenarris.COMPATIBLE_SPG_pos_overlap_list_set
    __swig_getmethods__["pos_overlap_list"] = _pygenarris.COMPATIBLE_SPG_pos_overlap_list_get
    if _newclass:pos_overlap_list = _swig_property(_pygenarris.COMPATIBLE_SPG_pos_overlap_list_get, _pygenarris.COMPATIBLE_SPG_pos_overlap_list_set)
    def __init__(self): 
        this = _pygenarris.new_COMPATIBLE_SPG()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _pygenarris.delete_COMPATIBLE_SPG
    __del__ = lambda self : None;
COMPATIBLE_SPG_swigregister = _pygenarris.COMPATIBLE_SPG_swigregister
COMPATIBLE_SPG_swigregister(COMPATIBLE_SPG)

class crystal(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, crystal, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, crystal, name)
    __repr__ = _swig_repr
    __swig_setmethods__["lattice_vectors"] = _pygenarris.crystal_lattice_vectors_set
    __swig_getmethods__["lattice_vectors"] = _pygenarris.crystal_lattice_vectors_get
    if _newclass:lattice_vectors = _swig_property(_pygenarris.crystal_lattice_vectors_get, _pygenarris.crystal_lattice_vectors_set)
    __swig_setmethods__["Xcord"] = _pygenarris.crystal_Xcord_set
    __swig_getmethods__["Xcord"] = _pygenarris.crystal_Xcord_get
    if _newclass:Xcord = _swig_property(_pygenarris.crystal_Xcord_get, _pygenarris.crystal_Xcord_set)
    __swig_setmethods__["Ycord"] = _pygenarris.crystal_Ycord_set
    __swig_getmethods__["Ycord"] = _pygenarris.crystal_Ycord_get
    if _newclass:Ycord = _swig_property(_pygenarris.crystal_Ycord_get, _pygenarris.crystal_Ycord_set)
    __swig_setmethods__["Zcord"] = _pygenarris.crystal_Zcord_set
    __swig_getmethods__["Zcord"] = _pygenarris.crystal_Zcord_get
    if _newclass:Zcord = _swig_property(_pygenarris.crystal_Zcord_get, _pygenarris.crystal_Zcord_set)
    __swig_setmethods__["atoms"] = _pygenarris.crystal_atoms_set
    __swig_getmethods__["atoms"] = _pygenarris.crystal_atoms_get
    if _newclass:atoms = _swig_property(_pygenarris.crystal_atoms_get, _pygenarris.crystal_atoms_set)
    __swig_setmethods__["spg"] = _pygenarris.crystal_spg_set
    __swig_getmethods__["spg"] = _pygenarris.crystal_spg_get
    if _newclass:spg = _swig_property(_pygenarris.crystal_spg_get, _pygenarris.crystal_spg_set)
    __swig_setmethods__["wyckoff_position"] = _pygenarris.crystal_wyckoff_position_set
    __swig_getmethods__["wyckoff_position"] = _pygenarris.crystal_wyckoff_position_get
    if _newclass:wyckoff_position = _swig_property(_pygenarris.crystal_wyckoff_position_get, _pygenarris.crystal_wyckoff_position_set)
    __swig_setmethods__["num_atoms_in_molecule"] = _pygenarris.crystal_num_atoms_in_molecule_set
    __swig_getmethods__["num_atoms_in_molecule"] = _pygenarris.crystal_num_atoms_in_molecule_get
    if _newclass:num_atoms_in_molecule = _swig_property(_pygenarris.crystal_num_atoms_in_molecule_get, _pygenarris.crystal_num_atoms_in_molecule_set)
    __swig_setmethods__["Z"] = _pygenarris.crystal_Z_set
    __swig_getmethods__["Z"] = _pygenarris.crystal_Z_get
    if _newclass:Z = _swig_property(_pygenarris.crystal_Z_get, _pygenarris.crystal_Z_set)
    __swig_setmethods__["Zp"] = _pygenarris.crystal_Zp_set
    __swig_getmethods__["Zp"] = _pygenarris.crystal_Zp_get
    if _newclass:Zp = _swig_property(_pygenarris.crystal_Zp_get, _pygenarris.crystal_Zp_set)
    def __init__(self): 
        this = _pygenarris.new_crystal()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _pygenarris.delete_crystal
    __del__ = lambda self : None;
crystal_swigregister = _pygenarris.crystal_swigregister
crystal_swigregister(crystal)


def generate_crystal(*args):
  return _pygenarris.generate_crystal(*args)
generate_crystal = _pygenarris.generate_crystal

def create_crystal_from_array(*args):
  return _pygenarris.create_crystal_from_array(*args)
create_crystal_from_array = _pygenarris.create_crystal_from_array

def print_crystal(*args):
  return _pygenarris.print_crystal(*args)
print_crystal = _pygenarris.print_crystal

def free_xtal(*args):
  return _pygenarris.free_xtal(*args)
free_xtal = _pygenarris.free_xtal

def c_check_structure(*args):
  return _pygenarris.c_check_structure(*args)
c_check_structure = _pygenarris.c_check_structure

def check_structure_with_vdw_matrix(*args):
  return _pygenarris.check_structure_with_vdw_matrix(*args)
check_structure_with_vdw_matrix = _pygenarris.check_structure_with_vdw_matrix

def generate_molecular_crystals_with_vdw_cutoff_matrix(*args):
  return _pygenarris.generate_molecular_crystals_with_vdw_cutoff_matrix(*args)
generate_molecular_crystals_with_vdw_cutoff_matrix = _pygenarris.generate_molecular_crystals_with_vdw_cutoff_matrix

def num_compatible_spacegroups(*args):
  return _pygenarris.num_compatible_spacegroups(*args)
num_compatible_spacegroups = _pygenarris.num_compatible_spacegroups
# This file is compatible with both classic and new-style classes.

