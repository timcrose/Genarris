

from ibslib.io import read,write
from ibslib.molecules import UniqueMolecules


# Define structure file to use
struct_file = "TETCEN.cif"
# Can change output format to "cif"
output_format = "geo"
# folder for molecule outputs
molecule_folder = "molecules"
# folder for whole molecule reconstructed structure files (rstruct)
rstruct_folder = "rstruct"
# overwrite controls if structure files can be overwritten 
overwrite = True
# Read in to Structure object
struct = read(struct_file)

# Identify unique molecules and construct whole molecules
um = UniqueMolecules(struct)
# Write single molecule found in structure
um.write_unique(molecule_folder, file_format=output_format,
                overwrite=overwrite)
# Write structure with whole molecules 
um.write_rstruct(rstruct_folder, file_format=output_format,
                 overwrite=overwrite)
