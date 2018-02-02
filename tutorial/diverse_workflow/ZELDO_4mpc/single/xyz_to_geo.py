import sys, glob, os

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
