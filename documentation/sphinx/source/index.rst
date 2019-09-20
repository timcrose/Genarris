.. Test documentation master file, created by
   sphinx-quickstart on Mon Sep 16 09:29:10 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Genarris 2.0 Beta Documentation
===============================

.. toctree::
   :maxdepth: 2
   :numbered:

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
    