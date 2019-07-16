import os, shutil
from setuptools import find_packages
from setuptools import setup
os.system('pip install spherical-functions')
setup(name = 'Genarris_2.0',
      description = 'Generation of Molecular Crystal Structures',
      author = 'T. Rose, T. Rithwik',
      url = 'http://noamarom.com',
      packages = find_packages(),
      install_requires = ['ase', 'pymatgen'],
      version = '2.0'
)
if os.path.exists('aimsutils'):
    shutil.rmtree('aimsutils')
os.mkdir('aimsutils')
os.system('git clone https://gitlab.lrz.de/theochem/aimsutils.git aimsutils')
shutil.copyfile('aimsutils_replacement_files/wigner_fast.py', 'aimsutils/aimsutils/rotate/wigner_fast.py')
shutil.copyfile('aimsutils_replacement_files/parser.py', 'aimsutils/aimsutils/parser.py')
shutil.copyfile('aimsutils_replacement_files/restarts.py', 'aimsutils/aimsutils/restarts.py')
shutil.copyfile('aimsutils_replacement_files/rotate2.py', 'aimsutils/aimsutils/rotate/rotate2.py')
shutil.copyfile('aimsutils_replacement_files/ylms.py', 'aimsutils/aimsutils/ylms.py')
os.chdir('aimsutils')
os.system('python setup.py install')
os.chdir('..')

os.system('git clone https://github.com/ritwit/cgenarris.git')
os.chdir('cgenarris')
os.system('python setup.py install')
