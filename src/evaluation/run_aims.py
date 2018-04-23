"""
If any part of this module is used for a publication please cite:

X. Li, F. Curtis, T. Rose, C. Schober, A. Vazquez-Mayagoitia, K. Reuter,
H. Oberhofer, and N. Marom "Genarris: Random Generation of Molecular Crystal 
Structures and Fast Screening with a Harris Approximation, ",
J. Chem. Phys., DOI: 10.1063/1.5014038; arXiv 1803.02145 (2018)
"""
'''
Created by Patrick Kilecdi 

This module handles the job calling of FHI-aims

'''

import shutil, os, time
import sys
from core import structure
from external_libs.filelock import FileLock
import ast, subprocess
import utilities.write_log as wl



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

def aims_single_run (inst):
	'''
	Calls and run one fhi-aims instance
	Returns True if the **execution** is successful, False if not
	Will not check the output file to see if the job is produces the desired result
	'''
	info_level = inst.get_info_level()
	working_dir = inst.get("aims_single_run","working_dir")
	bin = inst.get("FHI-aims","aims_bin_path")
    #aims_lib_path = inst.get("aims_single_run","aims_lib_path")
    execute_command = inst.get("FHI-aims","execute_command")
    execute_style = inst.get_with_default("FHI-aims","execute_style","safe_subprocess")
    multiple_launches = int(inst.get("FHI-aims","multiple_launches"))

    if execute_style == "safe_subprocess" or execute_style == "safe_system":
        update_poll_interval = int(inst.get("aims_single_run","update_poll_interval"))
        update_poll_times = int(inst.get("aims_single_run","update_poll_times"))

    if inst.has_option("aims_single_run","structure_path") and inst.get("aims_single_run","structure_path")!="": #Meaning a structure needs to be copied into inst
        structure_path = inst.get("aims_single_run","structure_path")
        structure_format = inst.get("aims_single_run","structure_format")
        if structure_format == "geometry":
            if structure_path!=os.path.join(working_dir,"geometry.in"):
                shutil.copyfile(structure_path,os.path.join(working_dir,"geometry.in"))
        elif molecule_format == "json":
            struct = structure.Structure()
            f=open(structure_path,"r")
            struct.load(f.read())
            f.close()
            f=open(os.path.join(output_dir,"geometry.in"))
            f.write(struct.get_geometry_atom_format())
            f.close()

    if inst.has_option("aims_single_run","control_path"): #Move the control.in to working_dir
        control_path=inst.get("aims_single_run","control_path")
        if control_path!=os.path.join(working_dir,"control.in"):
            shutil.copyfile(control_path,os.path.join(working_dir,"control.in"))

    if inst.has_option("aims_single_run","output_name"): #Allowing alternative output file name
        aimsout = os.path.join(working_dir,inst.get("aims_single_run","output_name"))
    else:
        aimsout = os.path.join(working_dir,"aims.out")

    if inst.has_option("aims_single_run","write_active_enabled"):
        if inst.get("aims_single_run","write_active_enabled")!="TRUE":
            raise ValueError("write_active_enabled set to an unknown value; for security, if present, it should be TRUE")
        write_active_enabled = True
        if inst.has_option("aims_single_run","write_active_path"):
            write_active_path = inst.get("aims_single_run","write_active_path")
        else:
            write_active_path = os.path.join(working_dir,"active.info")
    else:
        write_active_enabled = False
        write_active_path = ""

    if inst.has_option("aims_single_run","write_log_path"):
        write_log_enabled = True
        write_log_path = inst.get("aims_single_run","write_log_path")
        if inst.has_option("aims_single_run","write_log_name"):
            write_log_name = inst.get("aims_single_run","write_log_name")
        else:
            write_log_name = working_dir
    else:
        write_log_enabled = False

    if inst.has_option("FHI-aims","call_interval"): #Each binary call will be seperated by such value
        call_interval = float(inst.get("FHI-aims","call_interval"))
        if inst.has_option("FHI-aims","execute_info_path"):
            execute_info_path = inst.get("FHI-aims","execute_info_path")
        else:
            execute_info_path = os.path.join(inst.get("Genarris_master","working_dir"),"execute.info")
    else:
        call_interval = 0

    if inst.has_option("FHI-aims","execute_command_source"): #In case user want to launch another version of the execute command
        arglist = ["sh",inst.get("FHI-aims","execute_command_source")]
    else:
        arglist = [execute_command]

    if execute_command == "mpirun":
        arglist += ["-wdir",working_dir]
        if inst.has_option("aims_single_run","mpirun_processes"):
            arglist += ["-n",inst.get("aims_single_run","mpirun_processes")]
        if inst.has_option("aims_single_run","mpirun_hosts"):
            hostlist = inst.get_eval("aims_single_run","mpirun_hosts")
            hoststring = hostlist[0]
            for i in range (1,len(hostlist)):
                hoststring += ","+hostlist[i]
            arglist += ["--host",hoststring]
        arglist.append(bin)
        if info_level>=3:
            wl.write_master_log(inst,"This is arglist in aims_single_run "+str(arglist))

    elif execute_command == "runjob":
        block = inst.get("aims_single_run","runjob_block")
        nodes = int(inst.get("aims_single_run","runjob_nodes"))
        modes = int(inst.get("FHI-aims","runjob_modes"))
        thres = int(inst.get("FHI-aims","runjob_thres"))

    arglist += ["--np",str(modes*nodes),"-p",str(modes),"--envs","OMP_NUM_THREADS="+str(thres),"--verbose","INFO","--block",block,"--cwd",working_dir,"--exe",bin]

        if inst.has_option("aims_single_run","runjob_corner"):
            corner = inst.get("aims_single_run","runjob_corner")
            shape = inst.get("aims_single_run","runjob_shape")
            arglist += ["--corner",corner,"--shape",shape]
    else:
        raise ValueError("Unknown type of execute_command")

    if inst.has_option("aims_single_run","additional_arguments"):
        addi = inst.get_eval("aims_single_run","additional_arguments")
        arglist += [str(x) for x in addi]

    inst.grant_permission(working_dir)
    if execute_style == "fast_subprocess":
        #Uses subprocess.call, which will block until the command is complete
        outfile = open(aimsout,"w")
        errfile = open(inst.get_master_err_path(),"a")
        if call_interval!=0:
            get_execute_clearance(execute_info_path,write_active_path,call_interval)

        result = subprocess.call(arglist,stdout=outfile,stderr=errfile)
        inst.grant_permission(working_dir)
        if result==0:
            return True
        else:
            return False
        outfile.close()
        errfile.close()

    elif execute_style == "fast_system":
        command = arglist_to_command(arglist) + ">> " +aimsout
        if call_interval!=0:
            get_execute_clearance(execute_info_path,write_active_path,call_interval)

        os.system(command)
        inst.grant_permission(working_dir)
        return True #Unable to determine if the job is successfull with fast_system

    elif execute_style == "safe_subprocess":
        #Uses subprocess.Popen to launch to job and tracks the job's progress
        for i in range (multiple_launches):
            outfile = open(aimsout,"w")
            errfile = open(inst.get_master_err_path(),"a")
            if call_interval!=0:
                get_execute_clearance(execute_info_path,write_active_path,call_interval)
            if info_level>=3:
                wl.write_master_log(inst,"aims_single_run about to submit mpirun command in Popen")		

            p = subprocess.Popen(arglist,stdout=outfile,stderr=errfile)
            if write_log_enabled:
                write_log(write_log_path,"%s: Attempt %i to launch aims calculation" % (write_log_name,i))
                inst.grant_permission(write_log_path)
            time.sleep(1)

            if inst.has_option("aims_single_run","launch_time_out"):
                time_limit = int(inst.get("aims_single_run","launch_time_out"))
            else:
                time_limit = 60

            for j in range (time_limit): #Allow 60 seconds for aims to start outputting
                if (p.poll()!=None) or (os.stat(aimsout).st_size>512):
                    break
                if write_active_enabled:
                    write_active(write_active_path)
                    inst.grant_permission(write_active_path)
                inst.grant_permission(working_dir)
                time.sleep(1)

            if (os.stat(aimsout).st_size>512):
                if write_log_enabled:
                    write_log(write_log_path,write_log_name+": Aims calculation begins to output")
                    inst.grant_permission(write_log_path)
                break
            outfile.close()
            errfile.close()
            if write_log_enabled:
                write_log(write_log_path,write_log_name+": Aims calculation failed to launch")
                inst.grant_permission(write_log_path)
            try:
                p.send_signal(2)
            except:
                if write_log_enabled:
                    write_log(write_log_path,write_log_name+": Aims calculation unable to be killed")
                    inst.grant_permission(write_log_path)
                raise RuntimeError("Unable to kill subprocess for aims calculation")

            if write_active_enabled:
                active_sleep(60,write_active_path)
                inst.grant_permission(write_active_path)
            else:
                time.sleep(60)
            if i == multiple_launches-1:
                if write_log_enabled:
                    write_log(write_log_path,write_log_name+": WARNING, repeated failure to launch Aims calculation")
                    inst.grant_permission(write_log_path)

        #Now start monitoring the output to see if it is frequently updated
        counter=0; last=os.stat(aimsout).st_size
        while counter<update_poll_times and p.poll()==None: #The output file needs to update at least once in every 5 minutes
            if write_active_enabled:
                write_active(write_active_path)
                inst.grant_permission(write_active_path)
            inst.grant_permission(working_dir)
            time.sleep(update_poll_interval)
            if os.stat(aimsout).st_size>last:
                last=os.stat(aimsout).st_size
                counter=0
            else:
                counter+=1

        if counter==update_poll_times:
            if write_log_enabled:
                write_log(write_log_path,write_log_name+": Aims job hung")
                inst.grant_permission(write_log_path)
            try:
                p.send_signal(2)
            except:
                if write_log_enabled:
                    write_log(write_log_path,write_log_name+": Aims job cannot be killed")
                    inst.grant_permission(write_log_path)
                raise RuntimeError("Unable to kill aims job")
                if write_active_enabled:
                    active_sleep(60,write_active_path)
                    inst.grant_permission(write_active_path)
                else:
                    time.sleep(60)
            outfile.close()
            errfile.close()
            inst.grant_permission(working_dir)
            return False

        outfile.close()
        errfile.close()
        if write_log_enabled:
            write_log(write_log_path,write_log_name+": Aims job exited with status "+str(p.poll()))
            inst.grant_permission(write_log_path)

    else:
        raise ValueError("Unsupported mode of FHI-aims execute style "
				+ execute_style)
    if os.path.isfile(aimsout): 
        inst.grant_permission(aimsout)
    if os.path.isdir(working_dir):
        inst.grant_permission(working_dir)
	
def arglist_to_command(arglist):
    '''
    Given a list of arguments, concatenate them into a single string with one space between two arguments
    '''
    result=""
    for arg in arglist:
        result+=arg+" "
    return result	

def set_permission(working_dir):
    os.system("chmod -R g=u "+working_dir)

def write_active(file_path):
    f = open(file_path,"w")
    f.write(str(time.time()))
    f.close()
    os.system("chmod g=u "+file_path)

def active_sleep(total_time,folder,write_active_interval=10):
    '''
    Sleep for the total time but continues to write_active in the folder
    '''
    times = int(total_time/write_active_interval)
    for i in range (times):
        write_active(folder)
        time.sleep(write_active_interval)
    write_active(folder)
    try:
        time.sleep(total_time-write_active_interval*times)
        write_active(folder)
    except:
        pass

def write_log(log_path,message,time_stamp=True):
    dirname = os.path.dirname(log_path)
    filename = log_path[len(dirname)+1:]
    if time_stamp:
        message = time.strftime("%Y-%m-%d %H:%M:%S")+" "+message
    with FileLock(filename,dirname,10):
        f = open(log_path,"a")
        f.write(message+"\n")
        f.close()
        os.system("chmod g=u "+log_path)

def get_execute_clearance(execute_info_path,write_active_path="",buffer=3,max_wait=1000):
    '''
    Reads the execute_info_path as the time stamp of the last time the aims binary is called.
    Gets gets clearance for executing commands such as runjob and mpirun
    Will write_active to write_active_path if specified
    '''
    dirname = os.path.dirname(execute_info_path)
    filename = execute_info_path[len(dirname)+1:]
    for i in range (max_wait):
        with FileLock(filename,dirname,600):
            if not os.path.exists(execute_info_path):
                data_file=open(execute_info_path,"w")
                data_file.write(str(time.time()))
                data_file.close()
                os.system("chmod g=u "+execute_info_path)
                return True

            data_file=open(execute_info_path,"r")
            last_time=float(data_file.read())
            data_file.close()
            current_time=time.time()
            if (current_time-last_time>buffer) or (i==max_wait):
                data_file=open(execute_info_path,"w")
                data_file.write(str(time.time()))
                data_file.close()
                if i==max_wait and current_time-last_time<=buffer:
                    return False
                return True
        if write_active_path!="":
            write_active(write_active_path)
        time.sleep(buffer)