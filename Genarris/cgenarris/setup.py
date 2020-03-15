
"""
setup.py file 
"""

from distutils.core import setup, Extension


package = 'numpy'
try:
    __import__(package)
except ImportError:
    print("Please install numpy python package")
    exit()

import numpy

sources_spglib = ['arithmetic.c',
           'cell.c',
           'delaunay.c',
           'determination.c',
           'hall_symbol.c',
           'kgrid.c',
           'kpoint.c',
           'mathfunc.c',
           'niggli.c',
           'overlap.c',
           'pointgroup.c',
           'primitive.c',
           'refinement.c',
           'sitesym_database.c',
           'site_symmetry.c',
           'spacegroup.c',
           'spin.c',
           'spg_database.c',
           'spglib.c',
           'symmetry.c']

source_dir = "spglib_src"
include_dirs = [source_dir, ]
for i, s in enumerate(sources_spglib):
    sources_spglib[i] = "%s/%s" % (source_dir, s)


pygenarris = Extension('_pygenarris',include_dirs= ['./', numpy.get_include()], sources=['pygenarris.i', 'pygenarris.c', 
'combinatorics.c', 'molecule_placement.c', 'algebra.c', 'molecule_utils.c',
'spg_generation.c', 'lattice_generator.c', 'crystal_utils.c', 'check_structure.c', 'read_input.c', 'randomgen.c']+sources_spglib,
extra_compile_args=["-fopenmp", "-std=gnu99"], extra_link_args=['-lgomp'])

setup (name = 'pygenarris',
       version = '0.1',
       author      = "xxx",
       description = """yyy""",
       ext_modules = [pygenarris],
       py_modules = ["pygenarris"],
       )
