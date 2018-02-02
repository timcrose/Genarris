#!/bin/bash
#COBALT -n 1
#COBALT -t 60
#COBALT -q debug-chace-quad
#COBALT -A HybridPV_tesp

python ../../../../genarris/src/genarris_master.py volume_estimate_generation.conf
