import os, glob, json, shutil

'''
Purpose is to clean tmp dir of all structures that were calculated correctly
'''
root = '/home/ibier/genarris-runs/DUXYUQ/3_HA'
# Directory of structures already evaluated
evaluated_dir = os.path.join(root,'HA_structures')
# HA tmp directory that needs to be cleaned to save memory on computer
tmp_dir = os.path.join(root,'tmp_HA')
# Key for the harris energy in the json files
harris_energy_key = 'energy_harris'

def main():
    complete_list = complete_correct_list(evaluated_dir)
    print(len(complete_list))
    delete_tmp_dir(tmp_dir, complete_list)

def complete_correct_list(evaluated_dir):
    # Gets the list of structures that have a Harris Energy value
    complete_list = []
    for struct_file in os.listdir(evaluated_dir):
        struct_path = os.path.join(evaluated_dir, struct_file)
        struct_json = load_json(struct_path)
        # Check if HA has been completed correctly
        if struct_json['properties'][harris_energy_key] == False:
            continue
        else:
            struct_name = struct_file.replace('.json', '')
            complete_list.append(struct_name)
    return(complete_list)

def delete_tmp_dir(tmp_dir, complete_list):
    original_dir = os.getcwd()
    os.chdir(tmp_dir)
    for struct_name in complete_list:
        # glob returns data as list
        struct_tmp_dir = glob.glob('*'+struct_name+'*') 
        try: struct_tmp_dir = struct_tmp_dir[0]
        except: 
            print('Harris folder for {} may have already been removed'
                  .format(struct_name))
            continue
        struct_tmp_dir_abs = os.path.join(tmp_dir, struct_tmp_dir)
        print(struct_tmp_dir_abs)
        HA_check = check_HA_folder(struct_tmp_dir_abs)
        if not HA_check:
            print('Check that {} ran correctly.'
                  .format(struct_tmp_dir_abs))
            continue
        else:
            print('HA for {} ran correctly. Removing.'
                  .format(struct_tmp_dir_abs))
            shutil.rmtree(struct_tmp_dir_abs)
    os.chdir(original_dir)

def check_HA_folder(struct_tmp_dir_abs):
    HA_check = False
    harris_dir = os.path.join(struct_tmp_dir_abs, 'harris')
    aims_path = os.path.join(harris_dir, 'aims.out')        
    f = open(aims_path, 'r')
    for line in f.readlines():
        if 'density from restart information' in line:
            HA_check = True
            break
    f.close()
    return(HA_check)

def load_json(filepath):
    with open(filepath, 'r') as f:
        py_json = json.load(f)
    return(py_json)

if __name__ == '__main__':
    main()
