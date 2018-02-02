import os
structure_suffix = ".json"
structure_dir = "ZELDO_4_500_he_rcd"
rcd_output_file = "./rcd_difference_matrix_diverse_500_4_nmpc.info"

files = [x for x in os.listdir(structure_dir)
         if x[-len(structure_suffix):]==structure_suffix
         and x[-4:]!=".rcd"]

file_sorted = files.sort()
for x in files:
    f = open(os.path.join(structure_dir,x+".rcd"))
    l = f.read().split("\n")
    f.close()
    l.pop()
    diff = [float(x.split()[2]) for x in l]
    f = open(rcd_output_file,"a")
    f.write(" ".join(map(str,diff))+"\n")

    
    

