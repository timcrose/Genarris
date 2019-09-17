import os, shutil, platform
from setuptools import find_packages
from setuptools import setup
python_version = float(platform.python_version()[:3])
if python_version < 3.5:
    raise Exception('Need python 3.5 or above')
if os.path.exists('build'):
    shutil.rmtree('build')
if os.path.exists('Genarris_2.0.egg-info'):
    shutil.rmtree('Genarris_2.0.egg-info')
os.system('pip install spherical-functions')
setup(name = 'Genarris_2.0',
      description = 'Generation of Molecular Crystal Structures',
      author = 'T. Rose, T. Rithwik',
      url = 'http://noamarom.com',
      packages = find_packages(),
      install_requires = ['pymatgen', 'spglib'],
      version = '2.0'
)
cwd = os.path.abspath(os.getcwd())
package = 'mpi4py'
try:
    __import__(package)
except ImportError:
    os.system('pip install https://bitbucket.org/mpi4py/mpi4py/get/master.tar.gz')
ibslib_path = os.path.join('Genarris', 'ibslib')
if os.path.exists(ibslib_path):
    shutil.rmtree(ibslib_path)
os.chdir(os.path.dirname(ibslib_path))
os.system('git clone https://github.com/manny405/ibslib.git')
os.chdir(os.path.basename(ibslib_path))
os.system('python setup.py install')
os.chdir(cwd)
str_py_version = str(int(python_version * 10))
torch_whl = 'https://download.pytorch.org/whl/cpu/torch-1.1.0-cp' + str_py_version + '-cp' + str_py_version + 'm-linux_x86_64.whl'
os.system('pip install ' + torch_whl)
torch_whl = 'https://download.pytorch.org/whl/cpu/torchvision-0.3.0-cp' + str_py_version + '-cp' + str_py_version + 'm-linux_x86_64.whl'
os.system('pip install ' + torch_whl)
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
os.chdir('Genarris')
if os.path.exists('cgenarris'):
    shutil.rmtree('cgenarris')
os.system('git clone https://github.com/ritwit/cgenarris.git')
os.chdir('cgenarris')
os.system('python setup.py install')
os.chdir(cwd)
