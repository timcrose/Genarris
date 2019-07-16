# -*- coding: utf-8 -*-

from setuptools import find_packages
from distutils.core import setup

setup(
      name='ibslib',
      version='0.0',
      packages=['ibslib',
                'ibslib/analysis',
                'ibslib/crossover',
                'ibslib/descriptor',
                'ibslib/io',
                'ibslib/molecules',
                'ibslib/motif',
                'ibslib/structures',
                'ibslib/sgroup',
                'ibslib/database',
                'ibslib/calculators'
                ],
      #find_packages(exclude=[]),
      requires=['sklearn', 'pandas','pymatgen', 'scipy', 'pymongo']
      )
