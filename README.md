# Genarris
Molecular Crystal Structure Generation
To install
1) create a python virutal environment and direct your PYTHONPATH to include your python and site-packages directories. e.g. /home/trose/anaconda/anaconda3/envs/intelpython3_core/lib/python3.6:/home/trose/anaconda/anaconda3/envs/intelpython3_core/lib/python3.6/site-packages
2) Get permissions for getting cgenarris
3) run the following in this directory 
$ python setup.py install
or
$ python setup.py develop
4) Compile FHI-aims as a library according to the instructions in this package in a file called installation_instructions_arjuna.log. This is if you are on the Arjuna cluster. Otherwise, you will need to adapt the instructions for your cluster.
Please note that your previously installed ase version will be removed. This is why it's recommended to create a python virtual env in step 1.
