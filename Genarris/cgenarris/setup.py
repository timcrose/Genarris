
"""
setup.py file 
"""

from setuptools import setup, Extension
from distutils import sysconfig
import os

os.environ["CC"] = "mpicc" 
ldshared = sysconfig.get_config_var('LDSHARED')
#remove gcc or icc and paste mpicc for linker
ldshared = "mpicc " + ldshared.partition(' ')[2]
os.environ["LDSHARED"] = ldshared
print(ldshared)


package = 'numpy'
try:
    __import__(package)
except ImportError:
    print("Please install numpy python package")
    exit()
import numpy

package = 'mpi4py'
try:
    __import__(package)
except ImportError:
    print("Please install mpi4py python package")
    exit()
import mpi4py


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

pygenarris_mpi = Extension('_pygenarris_mpi',
                           include_dirs= ['./', numpy.get_include(), mpi4py.get_include()],
                           sources=['pygenarris_mpi.i', 'pygenarris_mpi.c', 'combinatorics.c', 'molecule_placement.c',
                           'algebra.c', 'molecule_utils.c','spg_generation.c', 'lattice_generator.c', 'crystal_utils.c',
                           'check_structure.c', 'read_input.c', 'randomgen.c']+sources_spglib,
                           extra_compile_args=["-std=gnu99"])

setup (name = 'pygenarris_mpi',
       version = '0.1',
       author      = "Rithwik Tom",
       description = """email:rtom@andrew.cmu.edu""",
       ext_modules = [pygenarris_mpi],
       py_modules = ["pygenarris_mpi"],
       )
