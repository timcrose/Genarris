

from ibslib.io.io_structures import import_geo, output_struct
from ibslib.analysis.sg_analysis import pymatgen_sg_struct

struct_path = 'relaxation/geometry.in.next_step'
struct_id = 'FUQJIK_relaxed'
output_json = 'FUQJIK_relaxed.json'
# 0.001 works well for structures with well defined space groups
# 0.85 works well when the space group of the structure is less well defined
symprec = 0.1

struct = import_geo(struct_path, struct_id)
psg = pymatgen_sg_struct(struct, symprec=symprec)
struct.set_property('space_group', psg)
output_struct(output_json, struct)
print('Pymatgen Space Group: {}'.format(psg))