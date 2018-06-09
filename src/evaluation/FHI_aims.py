"""
If any part of this module is used for a publication please cite:

X. Li, F. Curtis, T. Rose, C. Schober, A. Vazquez-Mayagoitia, K. Reuter,
H. Oberhofer, and N. Marom "Genarris: Random Generation of Molecular Crystal 
Structures and Fast Screening with a Harris Approximation, ",
J. Chem. Phys., DOI: 10.1063/1.5014038; arXiv 1803.02145 (2018)
"""
'''
Created by Patrick Kilecdi on 01/31/2017

Conducts FHI-aims calculations
'''


import os, random, time, shutil, sys
from core.structure import Structure
from utilities import parallel_run, misc
import copy
#sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../external_libs'))
sys.path.append('/home/trose/mpi4py_aimsdir')
import aims_w


__author__ = "Xiayue Li, Timothy Rose, Christoph Schober, and Farren Curtis"
__copyright__ = "Copyright 2018, Carnegie Mellon University and "+\
                "Fritz-Haber-Institut der Max-Planck-Gessellschaft"
__credits__ = ["Xiayue Li", "Luca Ghiringhelli", "Farren Curtis", "Tim Rose",
               "Christoph Schober", "Alvaro Vazquez-Mayagoita",
               "Karsten Reuter", "Harald Oberhofer", "Noa Marom"]
__license__ = "BSD-3"
__version__ = "180324"
__maintainer__ = "Timothy Rose"
__email__ = "trose@andrew.cmu.edu"
__url__ = "http://www.noamarom.com"

def print0(line):
    with open('Genarris.log', 'a') as f:
        if type(line) is str:
            f.write(line + '\n')
        else:
            for el in line:
                f.write(str(el) + ' ')
            f.write('\n')

def fhi_aims_batch_run(inst, comm):
    '''
    Main module that interacts with inst file
    '''
    print0(aims_w.__file__)
    
    sname = "fhi_aims_batch_run"
    structure_folder = inst.get(sname, "structure_folder")
    structure_suffix = inst.get_with_default(sname, "structure_suffix","")
    output_folder = inst.get(sname,"output_folder")
    output_suffix = inst.get_with_default(sname, "output_suffix", None,eval=True)
    output_format = inst.get_with_default(sname,"output_format",["json"],eval=True)
    #binary_path = inst.get(sname, "binary_path")
    control_path = inst.get(sname, "control_path")
    tmp_folder = inst.get_with_default(sname,"tmp_folder", "./FHI_aims_tmp")
    clean_tmp = not inst.get_boolean(sname,"disable_clean_tmp")
    #job_command = inst.get_with_default(sname, "job_command", "mpirun")
    #additional_args = inst.get_with_default(sname, "additional_args", None, eval=True)
    #nodes_per_partition = inst.get_with_default(sname, "nodes_per_partition", None, eval=True)
    #processes_per_partition = inst.get_with_default(sname, "processes_per_partition", None, 
    #                                                eval=True)
    #processes_per_node = inst.get_with_default(sname,"processes_per_node",None, eval=True)
    number_of_partitions = inst.get_with_default(sname, "number_of_partitions", 1, eval=True)
    #job_launch_wait = inst.get_with_default(sname, "job_launch_wait", None, eval=True)
    #job_update_wait = inst.get_with_default(sname, "job_update_wait", None, eval=True)
    absolute_success = inst.get_boolean(sname,"absolute_success")
    energy_key = inst.get_with_default(sname, "energy_key", "energy")

    split_comm = comm.Split(color=comm.Get_rank() % number_of_partitions)
    
    if split_comm.Get_rank() == 0:
        print0("Beginning FHI-aims batch relaxation run")
    
    conduct_relaxation(split_comm, structure_folder, structure_suffix,
                       output_folder, output_suffix, output_format,
                       control_path, tmp_folder, clean_tmp,
                       energy_key, absolute_success)
    
    split_comm.Free()

def conduct_relaxation(split_comm, structure_folder, structure_suffix,
                       output_folder, output_suffix, output_format, 
                       control_path, tmp_folder="./FHI_aims_tmp", clean_tmp=True,
                       energy_key="energy", absolute_success=False):
    
    list_of_files_remaining = get_list_of_files_remaining(structure_folder,
                                                          structure_suffix="")
    
    while list_of_files_remaining != []:
    
        struct, struct_path = _get_structure(structure_folder,
                                             list_of_files_remaining)
        
        calc_path = _submit_job(struct, struct_path, control_path, tmp_folder,
                                split_comm)
        
        if split_comm.Get_rank() == 0:
            print0("Job ended!")
        
            job_list = [struct, struct_path, calc_path, time.time(), 0]
        
            _conclude_job(job_list, output_folder, output_format, output_suffix,
                      absolute_success=absolute_success,energy_key=energy_key,
                      clean_tmp=clean_tmp)
    
        list_of_files_remaining = get_list_of_files_remaining(structure_folder,
                                                          structure_suffix="")
    
def fhi_aims_scavenge(inst):
    '''
    Main module that interacts with inst
    '''
    sname = "fhi_aims_scavenge"
    tmp_folder = inst.get_with_default(sname,"tmp_folder","./FHI_aims_tmp")
    wait_time = inst.get_with_default(sname,"wait_time",15,eval=True)
    clean_up = not inst.get_boolean(sname,"disable_cleanup")
    update_structure = not inst.get_boolean(sname,"disable_structure_update")
    update_energy_key = inst.get_with_default(sname,"update_energy_key","")
    scavenge_tmp(tmp_folder, wait_time, clean_up, update_structure, update_energy_key)

def scavenge_tmp(tmp_folder, wait_time=15, clean_up=True, update_structure=True, update_energy_key=""):
    if not os.path.isdir(tmp_folder):
        return
    folders = [os.path.join(tmp_folder,x) 
               for x in os.listdir(tmp_folder) 
               if os.path.isdir(os.path.join(tmp_folder,x))]

    print("Folders to be inspected: " + str(folders))
    for folder in folders:
        scavenge_folder(folder, wait_time=wait_time, update_structure=update_structure, update_energy_key=update_energy_key)
        if clean_up:
            misc.safe_rmdir(folder)

def scavenge_folder(folder, wait_time=15, update_structure=True, update_energy_key=""):

    print("Scavenging folder: ", folder)
   
    
    if time.time()-misc.read_active(folder)<wait_time or \
       os.path.isfile(os.path.join(folder,"success.info")) or \
       os.path.isfile(os.path.join(folder,"scavenged.info")) or \
       not os.path.isfile(os.path.join(folder,"struct.info")):
        print("Folder still active/already scavenged/was a success/struct file missing/")
        return False

    sfile = open(os.path.join(folder,"scavenged.info"),"w")
    sfile.write("Folder scavenged at: %s \n" % time.strftime("%Y-%m-%d %H:%M:%S"))

    f = open(os.path.join(folder,"struct.info"),"r")
    struct_path = f.read()
    f.close()
    if not os.path.isfile(struct_path) or os.path.isfile(struct_path+".fin") \
                                       or os.path.isfile(struct_path+".next"):
        sfile.write("Original structure missing, finished or updated\n")
        print("Original structure missing, finished or updated\n")
        sfile.close()
        return False

    _remove_job_status_file(struct_path)

    if not os.path.isfile(os.path.join(folder,"geometry.in.next_step")):
        _remove_job_status_file(struct_path)
        sfile.write("No geometry update happened\n")
        print("No geometry update happened")
        sfile.close()
        return True
            
#    try:
    struct = Structure()
    struct.loads_try_both(struct_path)
#    except:
#        sfile.write("Original structure file corrupted\n")
#        sfile.close()
#        print "Original structure file corrupted"
#        return False

#    try:
    if update_energy_key == "":
        update_energy_key = None
    _update_structure(struct, folder, energy_key = update_energy_key, create_duplicate=False)
#    except:
#        sfile.write("Attempt to update structure from .next_step failed")
#        sfile.close()
#        print "Attempt to update structure from .next_step failed"
#        return True

    print("Scavenging succeeded!")
    if update_structure:
        f = open(struct_path+".next","w")
        f.write(struct.dumps())
        f.close()


def _remove_job_status_file(struct_path):
    try:
        os.remove(struct_path+".lock")
    except:
        pass
    try:
        os.remove(struct_path+".failed")
    except:
        pass
    try:
        os.remove(struct_path+".hung")
    except:
        pass
        
def _is_running(jobs):
    for job in jobs:
        if job!=None:
            return True
    return False
'''
def _submit_job(struct, struct_path, binary, control_path, tmp_folder, 
                job_command, partition=None, additional_args=None):
    """
    Submits an FHI-aims job
    """
    if struct.struct_id!=None:
        name = struct.struct_id + "_" + misc.get_random_index()
    else:
        name = misc.get_random_index()

    calc_path = os.path.abspath(os.path.join(tmp_folder,name))
    print("This is calc_path: " + calc_path)
    misc.safe_make_dir(calc_path)
    
    f = open(os.path.join(calc_path, "geometry.in"),"w")
    f.write(struct.get_geometry_atom_format())
    f.close()

    print("This is control_path: ", control_path)
    shutil.copyfile(control_path,os.path.join(calc_path,"control.in"))

    f = open(os.path.join(calc_path, "struct.info"),"w")
    f.write(os.path.abspath(struct_path))
    f.close()

    misc.write_active(calc_path)

    args = [job_command]
    if job_command == "mpirun":
        if partition!=None:
            args += parallel_run._mpirun_arguments(partition)
        args += ["-wdir", calc_path]
        if additional_args!=None:
            args += additional_args
        args.append(binary)
    elif job_command == "aprun":
        if partition!=None:
            args += parallel_run._aprun_arguments(partition)
        if additional_args!=None:
            args += additional_args
        args.append(binary) 
    
    
    aimsout = open(os.path.join(calc_path,"aims.out"),"w")
    aimserr = open(os.path.join(calc_path,"aims.err"),"w")
    print "Submitting job: ", " ".join(map(str,args))
    current_dir = os.getcwd()
    os.chdir(calc_path)
    p = subprocess.Popen(args, stdout=aimsout, stderr=aimserr)
#    print "Job submitted in Popen"
    time.sleep(0.2)
#    print "Time waited"
    os.chdir(current_dir)
#    print "directory changed"

    return p, calc_path, aimsout, aimserr
'''

def _submit_job(struct, struct_path, control_path, tmp_folder, split_comm):
    """
    Submits an FHI-aims job
    """
    if split_comm.Get_rank() == 0:
        if struct.struct_id!=None:
            name = struct.struct_id + "_" + misc.get_random_index()
        else:
            name = misc.get_random_index()
    
        calc_path = os.path.abspath(os.path.join(tmp_folder,name))
    
    
        print0("This is calc_path: " + calc_path)
        misc.safe_make_dir(calc_path)
        
        f = open(os.path.join(calc_path, "geometry.in"),"w")
        f.write(struct.get_geometry_atom_format())
        f.close()
    
        print0("This is control_path: " + control_path)
        shutil.copyfile(control_path,os.path.join(calc_path,"control.in"))
    
        f = open(os.path.join(calc_path, "struct.info"),"w")
        f.write(os.path.abspath(struct_path))
        f.close()
    
        misc.write_active(calc_path)
    else:
        calc_path = None
    
    calc_path = split_comm.bcast(calc_path, root=0)
    
    split_comm.barrier()

    current_dir = os.getcwd()
    
    '''while not os.path.exists(calc_path):
        time.sleep(5)
        print(calc_path + ' DNE')
    '''
    os.chdir(calc_path)
    
    split_comm.barrier()
    
    commf = split_comm.py2f()
    
    split_comm.barrier()
    
    
    #os.system('ldd /home/trose/test/aims_py/aims_w.cpython-36m-x86_64-linux-gnu.so')
    #os.environ['LD_LIBRARY_PATH'] = '/opt/intel/compilers_and_libraries_2018.2.199/linux/mpi/intel64/lib:/opt/intel/compilers_and_libraries_2018.2.199/linux/mpi/mic/lib:/opt/intel/compilers_and_libraries/linux/mpi/intel64/lib:/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64:/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64_lin:/opt/intel/compilers_and_libraries/linux/mpi/intel64/lib:/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64::/opt/intel/compilers_and_libraries_2018.2.199/linux/mpi/intel64/lib:/opt/intel/compilers_and_libraries_2018.2.199/linux/mpi/mic/lib'
    
    aims_w.aims_w(commf)
    
    split_comm.barrier()
#    print "Job submitted in Popen"
    time.sleep(0.2)
#    print "Time waited"
    os.chdir(current_dir)
#    print "directory changed"
    return calc_path
     
def _conclude_job(job, output_folder, output_format, output_suffix, 
                  absolute_success=False, energy_key="energy", clean_tmp=True):
    '''
    Post-process a job
    '''
    
    struct = job[0]
    struct_path = job[1]
    calc_path = job[3]
    
    if not _is_successful(os.path.join(calc_path,"aims.out"), absolute_success):
        print0(time.strftime("%Y-%m-%d %H:%M:%S") + " : Job failed at " + calc_path)
        f = open(struct_path+".failed","a")
        f.write("Calculation failed at: " + time.strftime("%Y-%m-%d %H:%M:%S")+"\n")
        f.write(calc_path)
        f.close()
        return False
    else:
        f = open(struct_path+".fin","w")
        f.write("Successfully completed at: " + time.strftime("%Y-%m-%d %H:%M:%S")+"\n")
        f.write(calc_path)
        f.close()
        f = open(os.path.join(calc_path,"success.info"),"w")
        f.write("Successfully completed at: " + time.strftime("%Y-%m-%d %H:%M:%S")+"\n")
        f.write(calc_path)
        f.close()
        

    _update_structure(struct, calc_path, energy_key, False)
    try:
        os.remove(struct_path+".lock")
    except:
        pass

    misc.dump_structure(struct, output_folder, output_format, output_suffix)

    if clean_tmp:
        shutil.rmtree(calc_path)

def _handle_hung_job(job):
    '''
    Handles a job that is hung
    '''
    p = job[2]; calc_path = job[3]
    try:
        p.send_signal(2)
    except:
        pass
    
    p.wait()

    outfile = job[6]; outfile.close()
    errfile = job[7]; errfile.close()

    f = open(job[1]+".hung","w")
    f.write("Job hung at: " + time.strftime("%Y-%m-%d %H:%M:%S")+"\n")
    f.write(calc_path)
    f.close()
    
   
def _handle_failed_launch(job):
    '''
    Handles a job that fails to launch
    '''
    try:
        os.remove(job[1]+".lock")
    except:
        pass

    p = job[2]; calc_path = job[3]
    try:
        p.send_signal(2)
    except:
        pass

    p.wait()
 
    print(calc_path+": Job failed to launch at " + time.strftime("%Y-%m-%d %H:%M:%S")+"\n")
    
    outfile = job[6]; outfile.close()
    errfile = job[7]; errfile.close()
    
def get_list_of_files_remaining(folder, structure_suffix=""):
    list_of_files = [name for name in os.listdir(folder) \
                     if os.path.isfile(os.path.join(folder,name)) and 
                     (name[len(name)-len(structure_suffix):]==structure_suffix or
                     name[len(name)-5:]==".next") and
                     name[len(name)-5:]!=".lock" and name[len(name)-4:]!=".fin" and
                     not os.path.isfile(os.path.join(folder,name+'.lock')) and 
                     not os.path.isfile(os.path.join(folder,name+".fin")) and
                     not os.path.isfile(os.path.join(folder,name+".next"))]
    
    return list_of_files

def _get_structure(folder, list_of_files):
    '''
    Finds a structure currently under the tmp directory 
    Places a lock on it for use
    '''

    found = False
    while len(list_of_files) > 0:
        chosen = int(random.uniform(0, len(list_of_files)))
        name = list_of_files[chosen]

        if not os.path.isfile(os.path.join(folder, name+".lock")):

            f = open(os.path.join(folder, name+".lock"),"w")
            f.write(str(time.time()))
            f.close()

            f = open(os.path.join(folder,name),"r")
            struct = Structure()

            fread = f.read()
            try:
                struct.loads(fread)
            except:
                try:
                    struct.build_geo_whole_atom_format(fread)
                except:
                    print("Failed to load structure: " + 
                          os.path.join(folder, name))

            if struct.struct_id == None:
                l = name[:]
                while l[len(l)-5:]==".next":
                    l = l[:len(l)-5]
                struct.struct_id = l

            return struct, os.path.join(folder,name)

        else:
	    
            list_of_files=list_of_files[:chosen]+list_of_files[chosen+1:]

    return False, False

def _is_successful(aims_path,absolute_success=False):
    '''
    checks if relaxation/optimization was successful
    '''
    try:
        aims_out = open(aims_path,"r")
    except:
        return False #File missing

    counter = 0
    while True:
        line = aims_out.readline()
        if (not absolute_success and "Leaving FHI-aims" in line)\
            or "Have a nice day" in line:
                return True
        elif line == '':
            counter += 1
        else:
            counter = 0
        if counter > 10: #Allowing 10 empty lines in a row before determining eof
            break
    return False


def _update_structure(struct, calc_path, energy_key="energy", create_duplicate=True):
    if create_duplicate:
        struct = copy.deepcopy(struct)

    if os.path.isfile(os.path.join(calc_path,"geometry.in.next_step")):
        new_struct = Structure()
        f = open(os.path.join(calc_path,"geometry.in.next_step"), "r")
        new_struct.build_geo_whole_atom_format(f.read())
        f.close()
        struct.geometry = copy.deepcopy(new_struct.geometry)
        struct.properties["unit_cell_volume"] = new_struct.get_property("unit_cell_volume")
        struct.properties["cell_vol"] = new_struct.get_property("cell_vol")
        struct.properties["a"] = new_struct.get_property("a")
        struct.properties["b"] = new_struct.get_property("b")
        struct.properties["c"] = new_struct.get_property("c")
        struct.properties["lattice_vector_a"] = new_struct.get_property("lattice_vector_a")
#        print struct.properties["lattice_vector_a"]
        struct.properties["lattice_vector_b"] = new_struct.get_property("lattice_vector_b")
        struct.properties["lattice_vector_c"] = new_struct.get_property("lattice_vector_c")
   
    if energy_key!=None: 
        struct.properties[energy_key] = _extract_energy(os.path.join(calc_path,"aims.out"))
    return struct

def _extract_energy(path):
    '''
    reads the resulting energy from the FHI-aims log
    Returns: float if successful, False if unsuccessful
    '''
    aims_out = open(path)
    while True:
        line = aims_out.readline()
        if not line: return False  # energy not converged
        if '  | Total energy of the DFT / Hartree-Fock s.c.f. calculation      :' in line:
            tokens = line.split()
            energy = float(tokens[11])  # converts from SI string to float
            return energy


