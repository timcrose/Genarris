
from ibslib.io import read,write
from ibslib.io.extract import extract
from ibslib.io.aims_extractor import name_abs_path,name_from_path

# Input directory to extract from
struct_dir = "../../../rubrene/rubrene_relaxed_slabs/batch_run/200_2/"
# Output directory for structure files
output_dir = "test_rubrene_200_2"
# Type of output file
output_format = "json"
# name_func determins how a name for the structure will come from the file_path.
#   Users can define their own naming function which takes a file_path as an
#   argument and returns a name. name_abs_path and name_from_path are two 
#   example function. name_abs_path converts the absolute path into a name. 
#   name_from_path uses last directory as the name.
name_func = name_from_path
# All properties: "energy","vdw_energy","time","space_group","hirshfeld_volumes"
aims_property = ["energy", "vdw_energy", "hirshfeld_volumes"]
# Name for stored total energy value
energy_name = "energy"

# Definition of keyword arguments dictionary to be passed to extractor
extractor_kwargs = \
    {
        "aims_property": aims_property, 
        "energy_name": energy_name, 
        "log_file": None,
        "name_func": name_func
    }

# Extractor returns StructDict
struct_dict = extract(struct_dir, extractor="aims", 
                      extractor_kwargs=extractor_kwargs)

write(output_dir, struct_dict, file_format=output_format)
