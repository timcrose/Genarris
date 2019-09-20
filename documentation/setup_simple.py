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
setup(name = 'Genarris_2.0',
      description = 'Generation of Molecular Crystal Structures',
      author = 'T. Rose, T. Rithwik',
      url = 'http://noamarom.com',
      packages = find_packages(),
      install_requires = ['pymatgen', 'spglib'],
      version = '2.0'
)
