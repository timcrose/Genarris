.. Test documentation master file, created by
   sphinx-quickstart on Mon Sep 16 09:29:10 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Genarris 2.0 Beta Documentation
===============================

.. toctree::
   :maxdepth: 2
   :numbered:

How to use API Documentation
----------------------------

The execution of Genarris is controlled by a `configuration`_ file. The configuration
file specifies the execution of Genarris which is broken down into procedures, 
such as ``Pygenarris_Structure_Generation``. Each procedure has a corresponding
section in the configuration file, for our example ``pygenarris_structure_generation``.
The section contains options which control the operations performed by each 
procedure. 

This document details the options for procedures that are executed in the Genarris 2.0
*Robust* workflow. In order these are, ``Relax_Single_Molecule, Estimate_Unit_Cell_Volume, 
Pygenarris_Structure_Generation, Run_Rdf_Calc, Affinity_Propagation_Fixed_Clusters,
FHI_Aims_Energy_Evaluation, Affinity_Propagation_Fixed_Clusters, Run_FHI_Aims_Batch``.
There are many options that can be specified and modified for each section. 
All of these options are specified in this document under the
**Configuration File Options** section of each procedure.


There are three categories of **Configuration File Options**. These are *required*,
*optional*, and *inferred*. In the **Configuration File Options**, these categories 
are specified after the *type* of the option, such as *int*, *float*, or *bool*.
*Required* options have no category placed after the type. Both *optional* and 
*inferred* are specified after the type. *Optional* arguments are those that
have default settings that in general perform well. The user may specify these
*optional* arguments to have more control over the program executing. 
*Inferred* options are those that may be present in multiple different procedures. 
For example, the option ``aims_lib_dir`` is needed in the ``Relax_Single_Molecule``, 
``FHI_Aims_Energy_Evaluation``, and ``Run_FHI_Aims_Batch``. But, because it is 
an inferred parameter, it only needs to be specified once in the earliest procedure 
in which occurs and then it will be inferred by all further procedures. Options which 
are inferred are thus optional in all proceeding sections. 



..
    Description of the meaning of Sections, functions, Configuration file options, arguments.
    How the API ties all these together. 
    
    Most functions do not have many arguments. Control of the execution of the function is 
    typically controlled using an Instruct object which parses the configuration file. 
    Some functions may have many arguments, such as run_fhi_aims_batch. These arguments
    typically overlap with options which would typically be found in the configuration file. 
    However, these optional arguments can be provided to run it as a standalone function.
    
    Configuration file parameter infered parameters...


Genarris 2.0 Procedures for Robust Workflow
-------------------------------------------

.. autoclass:: Genarris.genarris_master.Genarris
    :members: Relax_Single_Molecule, 
              Estimate_Unit_Cell_Volume, 
              Pygenarris_Structure_Generation, 
              Run_Rdf_Calc,
              Affinity_Propagation_Fixed_Clusters,
              FHI_Aims_Energy_Evaluation,
              Run_FHI_Aims_Batch



Genarris 2.0 Callable Functions
-------------------------------

.. autofunction:: Genarris.evaluation.run_fhi_aims.run_fhi_aims_batch




.. Hypderlinks to be included in the document

.. _configuration: https://docs.python.org/3.4/library/configparser.html