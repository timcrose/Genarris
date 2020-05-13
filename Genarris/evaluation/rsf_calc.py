
from ibslib.io import read,write
from ibslib.acsf import Driver, RSF
from ibslib.libmpi import ParallelCalc

import torch, os, glob
import numpy as np
from scipy.spatial.distance import pdist, squareform
from Genarris.core.instruct import get_last_active_procedure_name

def rsf_calc(structure_path, comm=None, device=torch.device("cpu"),
                 acsf_kwargs=
                     {
                        "n_D_inter": 12,
                        "init_scheme": "shifted",
                        "cutoff": 12,
                        "eta_range": [0.05,0.5],
                        "Rs_range": [0.1,12],
                        "learn_rep": False
                     }, pdist_distance_type='euclidean', dist_mat_fpath='rsf_distance_matrix.dat',
                     output_dir='no_new_output_dir', save_envs=False, normalize_rsf_vectors=True,
                     standardize_distance_matrix=True):
    '''
    comm: MPI communicator or None

    structure_path: str
        Directory containing structures. 

    device: torch.device
        Keep at cpu for now

    acsf_kwargs: dict
        n_D_inter: int
            number of gaussians for each interaction
        init_scheme: str
            shifted works better
        cutoff: int
            angstrom distance cutoff for neighborhood of atoms
        eta_range: list of float
            ?
        Rs_range: list of numbers
            ?
        learn_rep: bool
            ?

    pdist_distance_type: str
        input parameter for the pdist function

    dist_mat_fpath: str
        path to file to write distance matrix to.

    output_dir: str
        Path of directory to write structures to (will create if it DNE). If 'no_new_output_dir' then
        input structures will be overwritten.

    save_envs: bool
        Whether to save the environment vectors

    normalize_rsf_vectors: bool
        Whether to normalize the rsf vectors over the columns of the feature 
        matrix before using them to compute the distance matrix. 
        Note that each rsf vector (row of the feature matrix) is not normalized.

    standardize_distance_matrix: bool
        Whether to standardize the distance matrix. The method is to simply divide all elements by the max
        value in the distance matrix. Because it is a distance matrix and thus all elements are positive,
        the standardized elements will be in the range [0, 1]

    Return None

    Purpose: calculation RSF vector for every atomic species pairwise combination. Update each json file 
        with its RSF vector and then write the structure pairwise distance matrix to a file for later use.

    Notes:
        1) This function is parallelized because a large number of files and/or a large number of different
            atomic species may benefit from parallelization. However, it has a serial option because in general,
            it is fast.
    '''
    if output_dir != 'no_new_output_dir':
        try:
            os.makedirs(output_dir)
        except:
            pass
    else:
        output_dir = structure_path
    
    struct_dict = read(structure_path)
    rsf = RSF(struct_dict, force=False, del_neighbor=True, 
              **acsf_kwargs)
    driver_kw = \
        {
            "file_format": "struct",
            "cutoff": acsf_kwargs["cutoff"],
            "prop_name": "",
        }
    driver = Driver(rsf, **driver_kw)
    
    
    if comm is None or comm.size == 1 or os.path.isfile(structure_path):
        # Serial version
        if comm is not None and comm.size > 1 and os.path.isfile(structure_path) and comm.rank > 0:
            return
        
        driver.output_dir = output_dir
        driver.calc(struct_dict)
        
    else:
        # Split up the work among ranks
        if not os.path.isdir(structure_path):
            raise Exception('structure_path must be a valid directory. Got:', structure_path)
        struct_files = glob.glob(os.path.join(structure_path, '*.json'))
        num_struct_files = len(struct_files)
        if comm.size > num_struct_files:

            if comm.rank >= num_struct_files:
                color = 1
            else:
                color = 0
            split_comm = comm.Split(color)
            if color == 1:
                return
        else:
            split_comm = comm
            
        pc = ParallelCalc(structure_path, output_dir, driver, comm=split_comm, 
                          verbose=False,
                          overwrite=True)
        pc.calc()
        
        split_comm.barrier()
        
        if split_comm.rank == 0:
            struct_dict = read(output_dir)
            
            ## Get dimension of features from that stored in the driver
            acsf_dim = driver.features.shape[0]
            
            ### Initialize array to store results
            all_rsf_vecs = np.zeros((len(struct_dict), acsf_dim))
            
            i = 0
            for struct_id,struct in struct_dict.items():
                all_rsf_vecs[i] = struct.properties["RSF"]
                i += 1
                
        ## Only want root rank to continue
        if split_comm.rank > 0:
            return
    
    if split_comm.rank > 0:
        return
    
    # Calculate pairwise distance matrix
    if normalize_rsf_vectors:
        norm_of_rsf_vecs_columnwise = np.linalg.norm(all_rsf_vecs, axis=0)
        all_rsf_vecs = all_rsf_vecs / norm_of_rsf_vecs_columnwise

    dist_mat = squareform(pdist(all_rsf_vecs, pdist_distance_type))
    if standardize_distance_matrix:
        dist_mat = dist_mat / np.max(dist_mat)
    
    print(dist_mat.shape)
    fp = np.memmap(dist_mat_fpath, dtype='float32', mode='write', shape=dist_mat.shape)
    fp[:] = dist_mat[:]
    

def run_rsf_calc(inst, comm):
    """
    Performs RSF calculation using ibslib.
    
    
    """
    sname = 'run_rsf_calc'
    last_section = get_last_active_procedure_name(inst, sname)

    sname_list = [sname, last_section, 
                  'fhi_aims_energy_evaluation', 
                  'harris_approximation_batch', 
                  'pygenarris_structure_generation', 
                  'structure_generation_batch']
    structure_dir = inst.get_inferred(sname, sname_list, 
        ['structure_dir'] + (len(sname_list) - 1) * ['output_dir'], type_='dir', required=True)

    n_D_inter = inst.get_with_default(sname, 'n_D_inter', 12, eval=True)
    init_scheme = inst.get_with_default(sname, 'init_scheme', 'shifted')
    cutoff = inst.get_with_default(sname, 'cutoff', 12, eval=True)
    eta_range = inst.get_with_default(sname, 'eta_range', [0.05,0.5], eval=True)
    Rs_range = inst.get_with_default(sname, 'Rs_range', [1,10], eval=True)
    pdist_distance_type = inst.get_with_default(sname, 'pdist_distance_type', 'euclidean')
    device = inst.get_with_default(sname, 'device', 'cpu')
    dist_mat_fpath = inst.get_with_default(sname, 'dist_mat_fpath', 'rsf_distance_matrix.dat')
    if not dist_mat_fpath.endswith('.dat'):
        raise Exception('Only supporting distance matrices saved as an np memmap with .dat extension')
    output_dir = inst.get_with_default(sname, 'output_dir', 'no_new_output_dir')
    save_envs = inst.get_boolean(sname, 'save_envs')
    normalize_rsf_vectors = inst.get_with_default(sname, 'normalize_rsf_vectors', 'TRUE')
    inst.set(sname, 'normalize_rsf_vectors', normalize_rsf_vectors)
    normalize_rsf_vectors = inst.get_boolean(sname, 'normalize_rsf_vectors')
    standardize_distance_matrix = inst.get_with_default(sname, 'standardize_distance_matrix', 'TRUE')
    inst.set(sname, 'standardize_distance_matrix', standardize_distance_matrix)
    standardize_distance_matrix = inst.get_boolean(sname, 'standardize_distance_matrix')

    acsf_kwargs = {
                        "n_D_inter": n_D_inter,
                        "init_scheme": init_scheme,
                        "cutoff": cutoff,
                        "eta_range": eta_range,
                        "Rs_range": Rs_range,
                     }

    rsf_calc(structure_dir, comm=comm, device=torch.device(device), 
             acsf_kwargs=acsf_kwargs, pdist_distance_type=pdist_distance_type, 
             dist_mat_fpath=dist_mat_fpath, output_dir=output_dir,
             save_envs=save_envs, normalize_rsf_vectors=normalize_rsf_vectors, 
             standardize_distance_matrix=standardize_distance_matrix)
