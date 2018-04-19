"""
If any part of this module is used for a publication please cite:

X. Li, F. Curtis, T. Rose, C. Schober, A. Vazquez-Mayagoitia, K. Reuter,
H. Oberhofer, and N. Marom "Genarris: Random Generation of Molecular Crystal 
Structures and Fast Screening with a Harris Approximation, ",
J. Chem. Phys., DOI: 10.1063/1.5014038; arXiv 1803.02145 (2018)
"""
import sys, glob, os



__author__ = "Xiayue Li, Timothy Rose, Christoph Schober, and Farren Curtis"
__copyright__ = "Copyright 2018, Carnegie Mellon University and "+\
                "Fritz-Haber-Institut der Max-Planck-Gessellschaft"
__credits__ = ["Xiayue Li", "Luca Ghiringhelli", "Farren Curtis", "Tim Rose",
               "Christoph Schober", "Alvaro Vazquez-Mayagoita",
               "Karsten Reuter", "Harald Oberhofer", "Noa Marom"]
__license__ = "BSD-3"
__version__ = "180324"
__maintainer__ = "Timothy Rose"
__email__ = "trose@andrew.cmu.edu"
__url__ = "http://www.noamarom.com"

def get_root_file_name(path_name):
    #get just the file_name without extension
    return path_name[path_name.rfind('/') + 1: path_name.rfind('.')]

def get_xyzs(xyz_names):

    xyzs = []
    for name in xyz_names:
        with open(name) as f:
            lines = f.readlines()
        coord_lines = [line.split() for line in lines if len(line.split()) == 4]
        xyzs.append(coord_lines)

    return xyzs

def mkdir_if_dne(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def write_geo(xyzs, xyz_names, out_dir):
    #write geometry.in file format

    for i, xyz in enumerate(xyzs):
        n_atoms = len(xyzs)
        xyz_root_name = get_root_file_name(xyz_names[i])
        with open(os.path.join(out_dir, xyz_root_name + '.in'), 'w') as f:
            for atom in xyz:
                f.write('atom ' + str(atom[1]) + ' ' + str(atom[2]) + ' ' + str(atom[3]) + ' ' + str(atom[0]) + '\n')

def main():
    #folder containing all xyz.in files
    xyz_dir = sys.argv[1]

    #get names of all xyz files
    xyz_names = glob.glob(xyz_dir + '/*.xyz')

    #get atomic coordinates and atomic symbol of each atom into a list of lists
    xyzs = get_xyzs(xyz_names)

    out_dir = sys.argv[2]

    mkdir_if_dne(out_dir)

    write_geo(xyzs, xyz_names, out_dir)

if __name__ == '__main__':
    main()
