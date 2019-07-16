

import os
from ibslib.io import read,write


struct_dir = "rubrene_200_2"
mbd_dir = "mbd_xyz_files"

struct_dict = read(struct_dir)
for struct_id,struct in struct_dict.items():
    print(struct_id)
    file_path = os.path.join(mbd_dir,struct_id)
    write(file_path, struct, file_format="mbd")
