# Genarris
Molecular Crystal Structure Generation
To install
1) create a python 3 virtual environment
e.g. see https://software.intel.com/en-us/articles/using-intel-distribution-for-python-with-anaconda
2) direct your PYTHONPATH to include your python and site-packages directories. e.g. /home/trose/anaconda/anaconda3/envs/intelpython3_core/lib/python3.6:/home/trose/anaconda/anaconda3/envs/intelpython3_core/lib/python3.6/site-packages
3) Get permissions for getting cgenarris
4) run the following in this directory 
$ python setup.py install
or
$ python setup.py develop
5) Compile FHI-aims as a library according to the instructions in this package in a file called installation_instructions_arjuna.log. This is if you are on the Arjuna cluster. Otherwise, you will need to adapt the instructions for your cluster.

Notes
1) Your previously installed ase version will be removed. This is why it's recommended to create a python virtual env in step 1.
2) The tutorial example has a preload_scripts.sh file which is called by the sub_genarris.sh slurm submission script. If your MPI is not Intel, you might need to call python directly in the submission script rather than calling ./preload_scripts.sh
