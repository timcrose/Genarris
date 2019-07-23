"""
If any part of this module is used for a publication please cite:

X. Li, F. Curtis, T. Rose, C. Schober, A. Vazquez-Mayagoitia, K. Reuter,
H. Oberhofer, and N. Marom "Genarris: Random Generation of Molecular Crystal 
Structures and Fast Screening with a Harris Approximation, ",
J. Chem. Phys., DOI: 10.1063/1.5014038; arXiv 1803.02145 (2018)
"""
'''
Handles all kinds of log writing
'''

import os,time
from Genarris.external_libs.filelock import FileLock


__author__ = "Xiayue Li, Timothy Rose, Christoph Schober, and Farren Curtis"
__copyright__ = "Copyright 2018, Carnegie Mellon University and "+\
                "Fritz-Haber-Institut der Max-Planck-Gessellschaft"
__credits__ = ["Xiayue Li", "Luca Ghiringhelli", "Farren Curtis", "Tim Rose",
               "Christoph Schober", "Alvaro Vazquez-Mayagoita",
               "Karsten Reuter", "Harald Oberhofer", "Noa Marom"]
__license__ = "BSD-3"
__version__ = "1.0"
__maintainer__ = "Timothy Rose"
__email__ = "trose@andrew.cmu.edu"
__url__ = "http://www.noamarom.com"

def set_global_output(stdout, stderr):
    global STDOUT
    STDOUT = stdout
    global STDERR
    STDERR = stderr

def write_log(log_path,message,time_stamp=True):
	dirname = os.path.dirname(log_path)
	filename = log_path[len(dirname)+1:]
	if time_stamp:
		message = time.strftime("%Y-%m-%d %H:%M:%S")+" "+message
	with FileLock(filename,dirname,600):
		f = open(log_path,"a")
		f.write(message+"\n")
		f.close()
		#os.system("chmod g=u "+log_path)

def print_time_log(message):
    STDOUT.write(time.strftime("%Y-%m-%d %H:%M:%S") + " " + message + "\n")

def print_time_warning(message):
    message = time.strftime("%Y-%m-%d %H:%M:%S") + \
            " WARNING: " + message + "\n"
    STDOUT.write(message)
    STDERR.write(message)

def print_time_error(message):
    message = time.strftime("%Y-%m-%d %H:%M:%S") + \
            " ERROR: " + message + "\n"
    STDOUT.write(message)
    STDERR.write(message)


def write_master_log(inst,message,time_stamp=True,additional_item=None):
	'''
	Writes log to Genarris_master.master_log_path
	Accepts additional item to be listed in message
	Should be a list [[section_1,option_1],[section_2,option_2],...]
	'''
	if additional_item != None:
		for section, option in additional_item:
			if inst.has_option(section,option):
				message += " "+section+"."+option+"="+inst.get(section,option)
	write_log(inst.get("Genarris_master","master_log_path"),message,time_stamp)
	inst.grant_permission(inst.get("Genarris_master","master_log_path"))

def write_master_err(inst,message,time_stamp=True,additional_item=None,write_master_log=True):
	'''
	Writes log to Genarris_master.master_log_path
	'''
	if additional_item != None:
		for section, option in additional_item:
			if inst.has_option(section,option):
				message += " "+section+"."+option+"="+inst.get(section,option)
	
	write_log(inst.get("Genarris_master","master_err_path"),message,time_stamp)
	inst.grant_permission(inst.get("Genarris_master","master_err_path"))

	if write_master_log:
		write_log(inst.get("Genarris_master","master_log_path"),"ERROR: "+message,time_stamp)
		inst.grant_permission(inst.get("Genarris_master","master_log_path"))

def write_master_err_and_raise(inst,message,time_stamp=True,additional_item=None,write_master_log=True):
	'''
	Writes error log and raises ValueError
	'''
	write_master_err(inst,message,time_stamp,additional_item,write_master_log)
	inst.grant_permission(inst.get("Genarris_master","master_err_path"))
	raise ValueError(message)
