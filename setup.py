import os, shutil
from setuptools import find_packages
from setuptools import setup
os.system('pip install spherical-functions')
setup(name = 'Genarris_2.0',
      description = 'Generation of Molecular Crystal Structures',
      author = 'T. Rose, T. Rithwik',
      url = 'http://noamarom.com',
      packages = find_packages(),
      install_requires = ['pymatgen'],
      version = '2.0'
)
cwd = os.path.abspath(os.getcwd())
os.chdir(os.path.join('Genarris', 'ibslib'))
os.system('python setup.py install')
os.chdir(cwd)
if os.path.exists('ase'):
   shutil.rmtree('ase')
os.system('pip uninstall -y ase')
os.system('git clone https://gitlab.com/ase/ase.git')
os.chdir('ase')
os.system('python setup.py install')
os.chdir(cwd)
if os.path.exists('aimsutils'):
    shutil.rmtree('aimsutils')
os.mkdir('aimsutils')
os.system('git clone https://gitlab.lrz.de/theochem/aimsutils.git aimsutils')
shutil.copyfile(os.path.join('aimsutils_replacement_files', 'wigner_fast.py'), os.path.join('aimsutils', 'aimsutils', 'rotate', 'wigner_fast.py'))
shutil.copyfile(os.path.join('aimsutils_replacement_files', 'parser.py'), os.path.join('aimsutils', 'aimsutils', 'parser.py'))
shutil.copyfile(os.path.join('aimsutils_replacement_files', 'restarts.py'), os.path.join('aimsutils', 'aimsutils', 'restarts.py'))
shutil.copyfile(os.path.join('aimsutils_replacement_files', 'rotate2.py'), os.path.join('aimsutils/aimsutils/rotate/rotate2.py'))
shutil.copyfile(os.path.join('aimsutils_replacement_files', 'ylms.py'), os.path.join('aimsutils', 'aimsutils', 'ylms.py'))
os.chdir('aimsutils')
os.system('python setup.py install')
os.chdir(cwd)

if os.path.exists('cgenarris'):
    shutil.rmtree('cgenarris')
os.system('git clone https://github.com/ritwit/cgenarris.git')
os.chdir('cgenarris')
os.system('python setup.py build_ext --inplace')
