.. Test documentation master file, created by
   sphinx-quickstart on Mon Sep 16 09:29:10 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Test's documentation!
================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Please look at the Genarris 2.0 Sections for all possible sections which have
been implemented in Genarris 2.0. 

How to use API Documentation
============================
Description of the meaning of Sections, functions, Configuration file options, arguments.
How the API ties all these together. 

Most functions do not have many arguments. Control of the execution of the function is 
typically controlled using an Instruct object which parses the configuration file. 
Some functions may have many arguments, such as run_fhi_aims_batch. These arguments
typically overlap with options which would typically be found in the configuration file. 
However, these optional arguments can be provided to run it as a standalone function.


Genarris 2.0 Sections
=====================

.. autofunction:: Genarris.evaluation.run_fhi_aims.run_fhi_aims_batch

