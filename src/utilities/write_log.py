'''
Handles all kinds of log writing
'''

import os,time
from external_libs.filelock import FileLock
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
