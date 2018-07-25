"""
If any part of this module is used for a publication please cite:

X. Li, F. Curtis, T. Rose, C. Schober, A. Vazquez-Mayagoitia, K. Reuter,
H. Oberhofer, and N. Marom "Genarris: Random Generation of Molecular Crystal 
Structures and Fast Screening with a Harris Approximation, ",
J. Chem. Phys., DOI: 10.1063/1.5014038; arXiv 1803.02145 (2018)
"""
import os

__author__ = "Xiayue Li, Timothy Rose, Christoph Schober, and Farren Curtis"
__copyright__ = "Copyright 2018, Carnegie Mellon University and "+\
                "Fritz-Haber-Institut der Max-Planck-Gessellschaft"
__credits__ = ["Xiayue Li", "Luca Ghiringhelli", "Farren Curtis", "Tim Rose",
               "Christoph Schober", "Alvaro Vazquez-Mayagoita",
               "Karsten Reuter", "Harald Oberhofer", "Noa Marom"]
__license__ = "BSD-3"
__version__ = "1.0"
__maintainer__ = "Timothy Rose"
__email__ = "trose@andrew.cmu.edu"
__url__ = "http://www.noamarom.com"

structure_suffix = ".rcd"
structure_dir = "./rcd_evaluated"
rcd_output_file = "./rcd_difference_matrix.out"

files = [x for x in os.listdir(structure_dir)
         if x[-len(structure_suffix):]==structure_suffix]

file_sorted = files.sort()
for x in files:
    print(x)
    f = open(os.path.join(structure_dir,x))
    l = f.read().split("\n")
    f.close()
    l.pop()
    diff = [float(value.split()[2]) for value in l]
    f = open(rcd_output_file,"a")
    f.write(" ".join(map(str,diff))+"\n")

    
    

