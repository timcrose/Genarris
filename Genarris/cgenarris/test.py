#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 16:36:06 2019

@author: ritwit
"""

from pygenarris_mpi import *
import numpy as np
from mpi4py import MPI

cutoff = np.load('cutoff_gly4mpc.npy');
cutoff=np.array(cutoff, dtype='float32');
print(cutoff);

seedstate = 39944;
filename='geometry.out';
num_structures = 10;
Z = 4;
volume_mean = 380;
volume_std = 18;
tol = 0.0001;
max_attempts = 2000000000;
comm = MPI.COMM_WORLD;

mpi_generate_molecular_crystals_with_vdw_cutoff_matrix(cutoff, num_structures, Z, volume_mean, volume_std, tol, max_attempts, comm);

