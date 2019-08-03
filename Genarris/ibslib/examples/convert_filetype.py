# -*- coding: utf-8 -*-


from ibslib.io import import_structures,import_geo,import_json,\
  output_struct_dict, output_struct, output_geo, check_overwrite \
  
                                    

###############################################################################
# Example: Convert JSON directory to FHI-aims geometry directory              #
###############################################################################
#struct_dir  = 'initial_pool_Arjuna_BrN'    # Input directory
#file_format = 'geometry'                   # Target file format
#output_dir  = 'initial_pool_Arjuna_geo'    # Output directory
#overwrite   = False                        # Overwrite directory if it exists
#
#struct_dict = import_structures(struct_dir)
#output_struct_dict(output_dir, struct_dict, file_format=file_format,
#                   overwrite=overwrite)

###############################################################################
# Example: Convert JSON structure file to FHI-aims geometry file              #
###############################################################################
#struct_file = 'IBUZIP_ccdc.json'           # Input geometry file
#struct_id   = 'IBUZIP'                     # Structure ID for structure object
#output_file = 'IBUZIP_ccdc.in  '           # Output JSON structure file name
#overwrite   = False                        # Overwrites file if it exists
#
#struct = import_json(struct_file)
#check_overwrite(output_file)
#output_geo(output_file, struct)

###############################################################################
# Example: Convert FHI-aims geometry file to JSON structure file              #
###############################################################################
#struct_file = 'IBUZIP_ccdc.in'             # Input geometry file
#struct_id   = 'IBUZIP'                     # Structure ID for structure object
#output_file = 'IBUZIP_ccdc.json'           # Output JSON structure file name
#overwrite   = False                        # Overwrites file if it exists
#
#struct = import_geo(struct_file, struct_id=struct_id)
#check_overwrite(output_file)
#output_struct(output_file, struct)
