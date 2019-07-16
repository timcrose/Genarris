from ibslib.io import read,write
from ibslib.molecules import UniqueMolecules

struct_path = "/Users/ibier/Research/Results/Molecule_Volume_Estimation/PAHs/PAHs_crystal/ABECAL/ABECAL.cif"
struct_test = read(struct_path)

um = UniqueMolecules(struct_test)
um.write_unique("Unique_Molecule", file_format="geo")
um.write_rstruct("rstruct", file_format="geo")
