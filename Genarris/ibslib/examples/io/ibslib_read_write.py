

from ibslib.io import read,write

# See help for detailed information
print(help(read))
print(help(write))

struct_dir = "cif"
struct_dict = read(struct_dir,file_format="")

output_dir = "geo"
file_format = "geo"
write(output_dir,struct_dict,"geo",overwrite=False)
