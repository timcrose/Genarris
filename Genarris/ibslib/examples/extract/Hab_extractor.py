

from ibslib.io import write
from ibslib.io.hab_extractor import hab_extractor

Hab_calc_path = "C:\\Users\\manny\\Research\\Datasets\\Hab_Dimers\\" \
                "FUQJIK_4mpc_tight\\batchrun"


test = hab_extractor(Hab_calc_path, extract_property="get_dimers")
dimer_dict = test.run_extractor()

write("C:\\Users\\manny\\Research\\Datasets\\Hab_Dimers\\" 
      "FUQJIK_4mpc_tight\\dimers", dimer_dict, overwrite=True)