# Genarris
Molecular Crystal Structure Generation

Genarris is a python package for generating molecular crystal structures in general and special positions. It's workflow includes relaxing the single molecule, estimating the solid state volume, generating structures, computing feature vector descriptors, clustering with Affinity Propagation, evaluating the energy of cluster exemplars, clustering the exemplars, and finally, relaxing the minimum energy structures from each of those clusters. The workflow is completely automated and uses MPI to parallelize each procedure. The DFT calculator we use is FHI-aims, but this may be simply exchanged for another calculator that accepts an MPI communicator.

Please see the installation_instructions.log file in the documentation directory and the tutorial in the tutorial directory.
