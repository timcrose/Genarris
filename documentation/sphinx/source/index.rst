.. Test documentation master file, created by
   sphinx-quickstart on Mon Sep 16 09:29:10 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Genarris 2.0 Beta Documentation
===============================

.. toctree::
   :maxdepth: 2
   :numbered:

Installation
------------
To install
1) Setup MPI and MKL
If already installed and modules exist, load them after unloading all conflicting modules. Note, in this installation tutorial we will use intel including intel's parallel studio package, but other program environments such as gnu will also work.
e.g.
module unload gnu
module unload openmpi
module load intel
module load impi

If MKL and MPI are already installed but modules do not exist, include the MPI and MKL directories in your environment variables.
e.g.
export $intel=/opt/ohpc/pub/intel/intel18/compilers_and_libraries_2018.3.222/linux  #Change to your parallel studio path
export PATH=$intel/mpi/intel64/bin_ohpc:$intel/mpi/intel64/bin:$intel/bin/intel64:$PATH
export LD_LIBRARY_PATH=$intel/mpi/intel64/lib:$intel/mpi/mic/lib:$intel/compiler/lib/intel64:$intel/compiler/lib/intel64_lin:$intel/ipp/lib/intel64:$intel/mkl/lib/intel64_lin:$intel/tbb/lib/intel64/gcc4.1:/opt/ohpc/pub/intel/intel18/debugger_2018/iga/lib:/opt/ohpc/pub/intel/intel18/debugger_2018/libipt/intel64/lib:$intel/daal/lib/intel64_lin:$intel/tbb/lib/intel64_lin/gcc4.4

Also export LD_PRELOAD to load the parallel studio MKL and Scalapack so importing FHI-aims and numpy does not cause conflict.
e.g.
export LD_PRELOAD=$intel/mkl/lib/intel64_lin/libmkl_intel_lp64.so:$intel/mkl/lib/intel64_lin/libmkl_sequential.so:$intel/mkl/lib/intel64_lin/libmkl_core.so:$intel/mkl/lib/intel64_lin/libmkl_blacs_intelmpi_lp64.so:$intel/mkl/lib/intel64_lin/libmkl_scalapack_lp64.so:$intel/mpi/intel64/lib/libmpi.so.12

2) create a python 3.5+ virtual environment
e.g.
export $anaconda=${HOME}/anaconda  #Change this to your desired anaconda install path
mkdir $anaconda
cd $anaconda

download and install anaconda
e.g.
wget https://repo.anaconda.com/archive/Anaconda3-2019.07-Linux-x86_64.sh
chmod +x Anaconda3-2019.07-Linux-x86_64.sh
./Anaconda3-2019.07-Linux-x86_64.sh

Include anaconda's binary in PATH
e.g.
export PATH=$anaconda/anaconda3/bin:$PATH

Make a python environment called e.g. genarris_env by installing intelpython3_core.
e.g.
conda config --add channels intel
conda create -n genarris_env intelpython3_core python=3

3) direct your path variables to include the new env
e.g.
export PYTHONPATH=$anaconda/anaconda3/envs/genarris_env/lib/python3.6:$anaconda/anaconda3/envs/genarris_env/lib/python3.6/site-packages:$PYTHONPATH
export PATH=$intel/mpi/intel64/bin_ohpc:$intel/mpi/intel64/bin:$intel/bin/intel64:$anaconda/anaconda3/envs/intelpython3_full/bin:$anaconda/anaconda3/bin:$PATH

4) Extract Genarris_v2.tar.gz into a desired directory and enter it
e.g.
export $genarris=${HOME}/genarris
mkdir $genarris
cp Genarris_v2.tar.gz $genarris
cd $genarris
tar -xzf Genarris_v2.tar.gz

5) Install Genarris. Note, one reason we recommend to create a python virutal env earlier is that running this installation script will remove the ase installation (if any) in the currently active python environment.
e.g.
cd $genarris/Genarris
python setup.py install

Genarris is now installed. We will first test that Genarris imports and MPI is working correctly with the following test and then the next step will be to compile FHI-aims as a python-importable library if you desire to use FHI-aims.
6) Test that Genarris imports and MPI is working correctly. Modify the submission script for your backend (here, we used slurm).
cd $genarris/documentation/mpi_and_genarris_test
sbatch mpi_and_genarris_test.sh

The desired output is that each rank reports a unique number.

7) Compile libaims into a python-importable library
# Set ulimit to avoid any possible memory problems
ulimit -s unlimited
ulimit -v unlimited

# Set OMP_NUM_THREADS to 1
export OMP_NUM_THREADS=1

# Obtain FHI-aims from https://aims-git.rz-berlin.mpg.de/aims/FHIaims If you don't have permissions, ask Volker Blum at volker.blum@duke.edu
export $aims=${HOME}/aims  #Change to your desired location for FHI-aims

# In its src directory ($aims/src), make sure the Makefile has all compilation flags (user defined settings) commented out
# Copy the make.sys file in the documentation directory of Genarris into FHI-aims' src directory. The make.sys is pasted here for reference.
cp $genarris/documentation/make.sys $aims/src
# Note, this make.sys assumes you are using intel's parallel studio and that your cluster's backend is intel. If this isn't the case, you'll need to set the flags accordingly.

# make.sys
###############
# Basic Flags #
###############
FC = mpiifort
FFLAGS = -O3 -ip -fp-model precise -fPIC
F90FLAGS = $(FFLAGS)
ARCHITECTURE = Generic
LAPACKBLAS = -L${MKLROOT}/lib/intel64 \
             -lmkl_intel_lp64 \
             -lmkl_sequential \
             -lmkl_core \
             -lmkl_blacs_intelmpi_lp64 \
             -lmkl_scalapack_lp64
F90MINFLAGS = -O0 -fp-model precise -fPIC

#########################
# Parallelization Flags #
#########################
USE_MPI = yes
MPIFC = ${FC}
SCALAPACK = ${LAPACKBLAS}

###############
# C,C++ Flags #
###############
CC = icc
CFLAGS = -O3 -ip -fp-model precise -fPIC

# EOF make.sys

# Compile FHI-aims as a shared library object:
cd $aims/src
make -j 20 libaims.scalapack.mpi
# where the "20" is however many cores you'd like to use for compilation

# Make a directory for compiling FHI-aims as a python library
e.g.
mkdir $aims/aims_as_python_lib
cd $aims/aims_as_python_lib

# Copy the Makefile and aims_w.f90 in the Genarris documentation directory to this directory. A copy of it has been pasted here for reference. Note that you will need to change the libaims version (currently shown as 190522). Again, you'll need to change the f90exec and/or fcompiler flags if your backend is not intel. aims_w.f90 is a wrapper script to interface with FHI-aims.
e.g.
cp $genarris/Genarris/documentation/Makefile $aims/aims_as_python_lib
cp $genarris/Genarris/documentation/aims_w.f90 $aims/aims_as_python_lib

# Makefile
LIBAIMS=${aims}/lib/libaims.190522.scalapack.mpi.so
include_dir=${anaconda}/anaconda3/envs/genarris_env/include

aims_w.so: aims_w.f90
	f2py --f90exec=mpiifort --fcompiler=intelem -m aims_w -c aims_w.f90 ${LIBAIMS} -I${include_dir}

clean:
	rm aims_w.*.so
# EOF Makefile

# compile FHI-aims as an importable python library!
make

8) Test that FHI-aims can run a job
# Modify the submission script in the $genarris/documentation/aims_test directory to run on your backend (here we used slurm).
export PYTHONPATH=$PYTHONPATH:$aims/aims_as_python_lib
cd $genarris/documentation/aims_test
sbatch aims_test.sh


Introduction
------------

Configuration File
^^^^^^^^^^^^^^^^^^

Genarris is a random crystal structure generation code that can be adapted to 
perform *ab initio* crystal structure prediction. The modularity of Genarris
is achieved through the sequential execution of procedures. The execution of 
Genarris is controlled by a `configuration`_ file. Below is a small example
of a configuration file for Genarris.::

    [Genarris_master]
    procedures = ["Pygenarris_Structure_Generation"]
    
    [pygenarris_structure_generation]
    # Path to the single molecule file to used for crystal structure generation
    molecule_path = relaxed_molecule.in
    # Number of cores (MPI ranks) to run this section with
    num_cores = 56
    # Number of OpenMP Threads
    omp_num_threads = 2
    num_structures = 5000
    Z = 4
    sr = 0.85
    tol = 0.00001
    max_attempts_per_spg_per_rank = 1000000000
    geometry_out_filename = glycine_4mpc.out
    output_format = json
    output_dir = glycine_4mpc_raw_jsons

**Sections** of the configuration file are denoted by square brakets, ``[...]``.
All parameters that are specified below a section are called **options**. The 
workflow of Genarris can be precisely controlled by the user by specifying the 
order of the desired procedures in ``[Genarris_master]``. The user must also
include the corresponding section for each procedure listed in 
``[Genarris_master]``. Each section may have many options which are required,
optional, or inferred.

This document details the options for procedures that are executed in the Genarris 2.0
*Robust* workflow. In order these are::
    
    ["Relax_Single_Molecule", 
     "Estimate_Unit_Cell_Volume",
     "Pygenarris_Structure_Generation", 
     "Run_Rdf_Calc", 
     "Affinity_Propagation_Fixed_Clusters",
     "FHI_Aims_Energy_Evaluation", 
     "Affinity_Propagation_Fixed_Clusters", 
     "Run_FHI_Aims_Batch"]
     
There are many options that can be specified and modified for each section. 
All of these options are specified in this document under the
**Configuration File Options** section of each procedure.


.. _category:

Option Category
^^^^^^^^^^^^^^^

There are three *categories* of **Configuration File Options**. These are *required*,
*optional*, and *inferred*. In the **Configuration File Options**, these categories 
are specified after the *type* of the option, such as *int*, *float*, or *bool*.

1. *Required* options have no category placed after the type in the 
   documentation. These options are required to be in the configuration 
   file for execution of Genarris. 

2. *Optional* arguments are specified after the option *type*. 
   These areguments have default settings built into the code perform 
   well in general. The user may specify these *optional* arguments 
   in the configuration file to have more control over the program 
   executing. 
     
3. *Inferred* options are specified after the option *type*. These options 
   may be present in multiple different procedures. For example, the option 
   ``aims_lib_dir`` is needed in the ``Relax_Single_Molecule``, 
   ``FHI_Aims_Energy_Evaluation``, and ``Run_FHI_Aims_Batch``. 
   But, because it is an inferred parameter, it only needs to be specified 
   once in the earliest procedure in which occurs and then it will be 
   inferred by all further procedures. Options which are inferred are thus 
   optional in all proceeding sections. 



..
    Description of the meaning of Sections, functions, Configuration file options, arguments.
    How the API ties all these together. 
    
    Most functions do not have many arguments. Control of the execution of the function is 
    typically controlled using an Instruct object which parses the configuration file. 
    Some functions may have many arguments, such as run_fhi_aims_batch. These arguments
    typically overlap with options which would typically be found in the configuration file. 
    However, these optional arguments can be provided to run it as a standalone function.
    
    Configuration file parameter inferred parameters...

..
    Add description of output file formats

.. 
    Add description of ibslib
    
.. Hypderlinks to be included in the document
    
.. _configuration: https://docs.python.org/3.4/library/configparser.html


Genarris 2.0 Procedures for Robust Workflow
-------------------------------------------

Description
^^^^^^^^^^^
This section details all arguments and configuration file
options for the procedures executed by the Robust Genarris 2.0 workflow. Each 
procedure is a class function of the of the ``Genarris`` master class.
The documentation follows a standard format for each procedure. The name
of the procedure is given first followed by a short description of the function 
the function it performs. Below the description is the the configuration file 
options subsection. This section gives the name, the data type, 
the :ref:`category`, and a description of each option which is accepted by the 
procedure. By referencing this documentation, the user can obtain precise 
control over the execution of Genarris procedures.

Genarris Procedures
^^^^^^^^^^^^^^^^^^^

.. autoclass:: Genarris.genarris_master.Genarris
    :members: Relax_Single_Molecule, 
              Estimate_Unit_Cell_Volume, 
              Pygenarris_Structure_Generation, 
              Run_Rdf_Calc,
              Affinity_Propagation_Fixed_Clusters,
              FHI_Aims_Energy_Evaluation,
              Run_FHI_Aims_Batch



.. 
    Genarris 2.0 Callable Functions
    -------------------------------
    
    .. autofunction:: Genarris.evaluation.run_fhi_aims.run_fhi_aims_batch
    
    
    

..
    Code Improvements
    -----------------
    Affinity propagation routine needs to be written in a general way to accept 
    two or more operations of clustering seamlessly. For this, I recommend to 
    allow for procedure names such as ``Affinity_Propagation_Fixed_Clusters`` for 
    only a single calculation and ``Affinity_Propagation_Fixed_Clusters_1`` and
    ``Affinity_Propagation_Fixed_Clusters_2``, and so on for more than one 
    calculation. This can be handled easily if the procedure name is parsed 
    before execution. Then, set ``sname`` and ``self.run_num`` in the ``APHandler`` 
    class to the corresponding value. This allows for more two executions of 
    AP in a simple way.
    
    Need to implement a default.conf file. Instruct will parse this file first
    and then parse the user provided configuration file, thus overwriting the 
    settings of default.conf if they are provided. The default procedures will
    be an empty list, but all default settings in their respective settings
    will be provided.
    