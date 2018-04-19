"""
If any part of this module is used for a publication please cite:

X. Li, F. Curtis, T. Rose, C. Schober, A. Vazquez-Mayagoitia, K. Reuter,
H. Oberhofer, and N. Marom "Genarris: Random Generation of Molecular Crystal 
Structures and Fast Screening with a Harris Approximation, ",
J. Chem. Phys., DOI: 10.1063/1.5014038; arXiv 1803.02145 (2018)
"""
'''
This set of utilities functions contain the necessary functions that carry out parallelization of the code
'''
import os,socket,time
from utilities import misc
import subprocess, math



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

def get_nodelist(system):
	'''
	Given the environment of the current Python code,
	Returns the list of available nodes in this job
	Currently only supports Cypress's SLURM system.
	'''
	def interp(s):
		if "-" in s:
			result=[]
			for i in range (int(s[:s.index("-")]),int(s[s.index("-")+1:])+1):
				result.append(str(i).zfill(3))
		else:
			result = [s]
		return result		
		
	if system == "SLURM":
		try:
			result = os.environ["SLURM_NODELIST"]
		except KeyError:
			try:
				result = os.environ["SLURM_JOB_NODELIST"]
			except:
				return [socket.gethostname()]
	else:
		return [socket.gethostname()]

	if not "[" in result:
		return [result]
	else:
		name = result[:result.index("[")]
		inside = result[result.index("[")+1:result.index("]")]
		if not "," in inside:
			tags = interp(inside)
		else:
			items = inside.split(",")
			tags =[]
			for item in items:
				tags += interp(item)
		return [name+tag for tag in tags]

def ssh_spawn_python_subprocess(hosts_and_insts,modules=[],additional_commands="",wait_time=3):
	'''
	Given a list of [[host_1,inst_1],[host_2,inst_2],...]
	And a list of modules to be launched after the ssh occurs
	For each inst, must have Genarris_master, working_dir and script_path options specified
	Spawn replica Genarris_master subprocesses
	And return the list of subprocesses for monitoring
	'''
	processes = []
	tmp_inst_list = []
	for (host,inst) in hosts_and_insts:
		#First will generate a temporary file 
		misc.safe_make_dir(inst.get("Genarris_master","working_dir"))
		tmp_inst_file = os.path.join(inst.get("Genarris_master","working_dir"),misc.get_random_index()+".conf")
		f = open(tmp_inst_file,"w")
		inst.write(f)
		f.close()
		tmp_inst_list.append(tmp_inst_file)

		tmp_exec_file = os.path.join(inst.get("Genarris_master","working_dir"),misc.get_random_index()+".exe")
		f = open(tmp_exec_file,"w")
		f.write("#!/bin/sh \n ssh %s << EOF \n" % host)
		if len(modules)>0:
			for module in modules:
				f.write("module load "+module+"\n")
		f.write(additional_commands+"\n")
		f.write("python %s %s  \n" % (inst.get("Genarris_master","script_path"),os.path.abspath(tmp_inst_file)))
		f.close()

		tmp_inst_list.append(tmp_exec_file)	
		
		#Now spawn a new python instance
#		p = subprocess.Popen(["ssh",host,"python "+inst.get("Genarris_master","script_path")+" "+os.path.abspath(tmp_inst_file)])
		p = subprocess.Popen(["sh",tmp_exec_file])
		processes.append(p)

	time.sleep(wait_time)
	for tmp_inst_file in tmp_inst_list:
		os.remove(tmp_inst_file)
	return processes

def get_partitions(command, hostlist=None, nodes_per_partition=None,
                   processes_per_partition=None, number_of_partitions=None):
    if command=="mpirun":
        processes = _get_all_processes(command, hostlist)
        return _partition_processes(processes, nodes_per_partition,
                                    processes_per_partition, number_of_partitions)
    elif command=="aprun":
        if processes_per_partition==None or number_of_partitions==None:
            raise ValueError("aprun demands setting processes_per_partition and number_of_partitions")
        

def _get_all_processes(command,hostlist=None):
    '''
    This function returns all the processors available for the current job
    Specify hostlist if only wishes to access certain nodes
    '''
    if command!="mpirun" and command!="srun":
        raise ValueError("Unsupported command for get_all_hosts; only supporting mpirun and srun")
    arglist = [command]
    if hostlist!=None:
        if command == "mpirun":
            arglist += ["--host",",".join(map(str,hostlist))]
        elif command == "srun":
            arglist += ["-w",",".join(map(str,hostlist))]

    print_host = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "print_host.py")
    arglist += ["python",print_host]
    print "This is arglist in _get_all_processes: " + str(arglist)
    
    p = subprocess.Popen(arglist,stdout=subprocess.PIPE)
    time.sleep(2)
    p.wait()
    out , err = p.communicate()

    try:
        out = str(out,"utf-8") #Python 3
    except:
	pass

    hosts = out.split("\n")
    hosts.pop() #Last line empty
    print "These are the hosts:" + str(hosts)
    return hosts

def _partition_processes(processes, nodes_per_partition=None, 
                         processes_per_partition=None, number_of_partitions=None):
    '''
    Given a list of processes, partition them according to given commands
    In the comment below, a partition is also called a replica
    ''' 	
    nodes = list(set(processes))
    nodes.sort()
    noh = len(processes) #Number of processes/hosts
    non = len(nodes) #Number of nodes
    hpn = int(noh/non) #Processes/hosts per node
    npp = [] #Nodes per partition
    ppp = [] #Processes per partition
    nop = 0 #Number of partitions

    if nodes_per_partition!=None:
        npp = nodes_per_partition
        if processes_per_partition!=None:
            ppp = processes_per_partition
        else:
            ppp = npp*hpn

        if number_of_partitions!=None:
            nop = number_of_partitions
        else:
            nop = non/npp
        npp = [npp]*nop
        ppp = [ppp]*nop

    elif processes_per_partition!=None:
        ppp = processes_per_partition

        if ppp >= hpn:
            npp = int(math.ceil(ppp/(hpn+0.0)))

            if number_of_partitions!=None:
                nop = number_of_partitions
            else:
                nop = max(non/npp,1)
            npp = [npp]*nop
            ppp = [ppp]*nop

        else:
            ppn = hpn/ppp
            if number_of_partitions!=None:
                nop = number_of_partitions
            else:
                nop = non*ppn
            npp = [int(0+(x%ppn)==0) for x in range(1,nop+1)]
	    ppp = [ppp]*nop

    elif number_of_partitions!=None:
        nop = number_of_partitions
        if nop > non: #Multiple replicas have to share same node
            ppn = nop / non
            add = nop % non #Certain nodes may have to handle more replicas
            npp = ([0]*ppn+[1])*add + ([0]*(ppn-1)+[1])*(non-add)
            #Only the last replica on the node gets a 1
            #Later, when assigning nodes, this is when the next node is accessed
            #First assign the number of processes for each replica on the extra-loaded nodes

            ppp = ([hpn/(ppn+1)+1]*(hpn%(ppn+1))\
				+[hpn/(ppn+1)]*(ppn+1-(hpn%(ppn+1))))*add

            ppp+= ([hpn/ppn+1]*(hpn%ppn)+[hpn/ppn]*(ppn-(hpn%ppn)))*(non-add)

        else:
	    npp = non / nop
            add = non % nop
            npp = [npp+1]*add + [npp]*(nop-add)
            ppp = [hpn*x for x in npp]
    else:
        raise KeyError("No partition variable is specified")

    print("Number of parallel replicas: "+str(nop))
    print("Nodes assigned to each replica (0 indicates that this replica is assigned a fraction of a node and is not the last replica on the node: " + " ".join(map(str,npp)))
    print("Processes assigned to each replica: "+" ".join(map(str,ppp)))
#
    print("Total available nodes: %i; assigned nodes: %i" % (non,sum(npp)))

    if sum(npp) > non:
        print("Oversubscription of node has occured; Check the compatibility of number_of_replicas and nodes_per_replica/processes_per_node; aborting...")
        raise ValueError("Oversubscription of node has occured; Check the compatibility of number_of_replicas and nodes_per_replica/processes_per_node")

    elif sum(npp) < non:
        print("Not all nodes are utilized; try setting only number_of_replicas or setting a replica_per_node that divides the total number of nodes")

    print("Total available processes: %i; assigned processes: %i" % (noh,sum(ppp)))

    if sum(ppp)>noh:
        print("Oversubscription of processes has occured; This should be avoided in general, but GAtor will proceed for now.")
    elif sum(ppp)<noh:
        print("Undersubscription of processes has occured; If wish to utilize all processes, try setting only number_of_replicas or setting a processes_per_replica that fits the processes_per_node")

    partitions = []
    n_index = 0
    for i in range(nop):
        if npp[i] <= 1:
            partitions.append([(nodes[n_index],ppp[i])])
        else:
            partitions.append([(x, hpn) for x in nodes[n_index:n_index+npp[i]]])
        
        n_index += npp[i]
    
    return partitions

def _mpirun_arguments(partition):
    '''
    Given a partition, returns the proper arguments to set to mpirun
    '''
    args = ["-host",",".join(map(str,[x[0] for x in partition]))]
    args += ["-n",str(sum([x[1] for x in partition]))]
    return args 

def _aprun_arguments(partition):
    args = ["-n",str(partition[0]*partition[1]),"-N",str(partition[1])]
    return args    

#processes = ["cypress-01"]*20+["cypress-02"]*20 + ["cypress-03"]*20
#print(_partition_processes(processes,nodes_per_partition=None,number_of_partitions=2,processes_per_partition=None))

