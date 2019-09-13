# Genarris
Molecular Crystal Structure Generation
To install
1) Setup MPI and MKL
If already installed and modules exist, load them.
e.g. (gnu can be used instead of intel)
module load intel
module load impi
If already installed and modules do not exist, include the mpi and mkl directories in your environment variables.
e.g. if $intel is /opt/ohpc/pub/intel/intel18/compilers_and_libraries_2018.3.222/linux then
export PATH=$intel/mpi/intel64/bin_ohpc:$intel/mpi/intel64/bin:$intel/bin/intel64:$PATH
export LD_LIBRARY_PATH=$intel/mpi/intel64/lib:$intel/mpi/mic/lib:$intel/compiler/lib/intel64:$intel/compiler/lib/intel64_lin:$intel/ipp/lib/intel64:$intel/mkl/lib/intel64_lin:$intel/tbb/lib/intel64/gcc4.1:/opt/ohpc/pub/intel/intel18/debugger_2018/iga/lib:/opt/ohpc/pub/intel/intel18/debugger_2018/libipt/intel64/lib:$intel/daal/lib/intel64_lin:$intel/tbb/lib/intel64_lin/gcc4.4

Additional environment variables
export LD_PRELOAD=$intel/mkl/lib/intel64_lin/libmkl_intel_lp64.so:$intel/mkl/lib/intel64_lin/libmkl_sequential.so:$intel/mkl/lib/intel64_lin/libmkl_core.so:$intel/mkl/lib/intel64_lin/libmkl_blacs_intelmpi_lp64.so:$intel/mkl/lib/intel64_lin/libmkl_scalapack_lp64.so:$intel/mpi/intel64/lib/libmpi.so.12
2) create a python 3.5+ virtual environment
e.g.
download and install anaconda (after changing directory to some desired path denoted as $anaconda)
e.g.
wget https://repo.anaconda.com/archive/Anaconda3-2019.07-Linux-x86_64.sh
chmod +x Anaconda3-2019.07-Linux-x86_64.sh
./Anaconda3-2019.07-Linux-x86_64.sh

Include anaconda's binary in PATH
e.g.
export PATH=$anaconda/anaconda3/bin:$PATH

Follow instructions at the following website to get intelpython3_core
https://software.intel.com/en-us/articles/using-intel-distribution-for-python-with-anaconda
3) direct your path variables to include your python and site-packages directories. e.g.
export PYTHONPATH=$anaconda/anaconda3/envs/intelpython3_core/lib/python3.6:$anaconda/anaconda3/envs/intelpython3_core/lib/python3.6/site-packages:$PYTHONPATH

Reset PATH to use the python virtual environment
export PATH=$intel/mpi/intel64/bin_ohpc:$intel/mpi/intel64/bin:$intel/bin/intel64:$anaconda/anaconda3/envs/intelpython3_full/bin:$anaconda/anaconda3/bin:$PATH

4) run the following in this directory 
python setup.py install
5) Compile FHI-aims as a library according to the instructions in this package in a file called installation_instructions_arjuna.log. This is if you are on the Arjuna cluster. Otherwise, you will need to adapt the instructions for your cluster.

Notes
1) Your previously installed ase version will be removed and replaced with the git version. This is why it's recommended to create a python virtual env in step 2.
2) The tutorial example has a preload_scripts.sh file which is called by the sub_genarris.sh slurm submission script. If your MPI is not Intel, you might need to call python directly in the submission script rather than calling ./preload_scripts.sh
3) If you don't specify volume parameters like volume_mean and volume_std, and do have the Estimate_Unit_Cell_Volume procedure before a structure generation procedure in the procedures list, then volume parameters will be estimated and the estimates will be used by the generation procedure and printed to the Genarris.log file.
4) Some procedures, such as the ones that call FHI-aims, use one rank as a master rank. e.g. if you say num_cores = 57 in run_fhi_aims_batch and num_replicas = 1 then 56 cores will work on an aims job at a time.
5) Some options are listed in the example ui.conf as under the inferred header. These will be gotten from other defined sections in a reasonable order if not defined in the current section. Usually the order is in the order of the workflow, though sometimes it is in the reverse order when it makes sense. e.g. AP comes before RCD when run_fhi_aims_batch is trying to infer structure_dir as the output from those sections.
