from ibslib.analysis.diversity import DiversityAnalysis
from ibslib.io import read,write
import torch, os, glob
import numpy as np
from scipy.spatial.distance import pdist, squareform
from Genarris.core.instruct import get_last_active_procedure_name

def rdf_calc(structure_path, comm=None, device=torch.device("cpu"),
                 acsf_kwargs=
                     {
                        "n_D_inter": 12,
                        "init_scheme": "shifted",
                        "cutoff": 12,
                        "eta_range": [0.05,0.5],
                        "Rs_range": [0.1,12],
                        "learn_rep": False
                     }, pdist_distance_type='euclidean', dist_mat_fpath='rdf_distance_matrix.dat',
                     output_dir='no_new_output_dir', save_envs=False, normalize_rdf_vectors=True,
                     standardize_distance_matrix=True):
    '''
    comm: MPI communicator or None

    structure_path: str
        Either file containing a structure or directory containing structures with
        .json file extension

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

    normalize_rdf_vectors: bool
        Whether to normalize the rdf vectors over the columns of the feature matrix before using them to compute the
        distance matrix. Note that each rdf vector (row of the feature matrix) is not normalized.

    standardize_distance_matrix: bool
        Whether to standardize the distance matrix. The method is to simply divide all elements by the max
        value in the distance matrix. Because it is a distance matrix and thus all elements are positive,
        the standardized elements will be in the range [0, 1]

    Return None

    Purpose: calculation RDF vector for every atomic species pairwise combination. Update each json file 
        with its RDF vector and then write the structure pairwise distance matrix to a file for later use.

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
    if comm is None or comm.size == 1 or os.path.isfile(structure_path):
        # Serial version
        if comm is not None and comm.size > 1 and os.path.isfile(structure_path) and comm.rank > 0:
            return
        s = read(structure_path)
        da = DiversityAnalysis(device=device, acsf_kwargs=acsf_kwargs)
        all_rdf_vecs,_ = da.calc(s)
        all_rdf_vecs = all_rdf_vecs.detach().numpy()
        if os.path.isfile(structure_path):
            s.properties['rdf'] = list(np.array(all_rdf_vecs, dtype=str))
            if output_dir != 'no_new_output_dir':
                struct_file = os.path.join(output_dir, struct_file)
            if not save_envs:
                del s.properties['R']
            write(struct_file, s, file_format='json', overwrite=True)
            return
        else:
            struct_files = glob.glob(os.path.join(structure_path, '*.json'))
            for i,struct_id in enumerate(s):
                s[struct_id].properties['rdf'] = list(np.array(all_rdf_vecs[i,:], dtype=str))
                struct_file = struct_files[i]
                if output_dir != 'no_new_output_dir':
                    struct_file = os.path.join(output_dir, os.path.basename(struct_file))
                if not save_envs:
                    del s[struct_id].properties['R']
                write(struct_file, s[struct_id], file_format='json', overwrite=True)
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
        num_structs_per_rank = int(num_struct_files / split_comm.size)
        #print('num_structs_per_rank', num_structs_per_rank, flush=True)
        #print('split_comm.rank, split_comm.size', split_comm.rank, split_comm.size, flush=True)
        #print('num_struct_files', num_struct_files, flush=True)
        if split_comm.rank == split_comm.size - 1:
            my_struct_files = struct_files[split_comm.rank * num_structs_per_rank :]
        else:
            my_struct_files = struct_files[split_comm.rank * num_structs_per_rank : (split_comm.rank + 1) * num_structs_per_rank]
        # Read in my json files
        my_struct_dict = {}
        for struct_file in my_struct_files:
            struct = read(struct_file)
            my_struct_dict[struct.struct_id] = struct
        # calcluate RDF
        #print('my_struct_dict', my_struct_dict, flush=True)
        da = DiversityAnalysis(device=device, acsf_kwargs=acsf_kwargs)
        my_rdf,_ = da.calc(my_struct_dict)
        #print('my_rdf', my_rdf, flush=True)
        #print('before', type(my_rdf), flush=True)
        my_rdf = my_rdf.detach().numpy()
        #print('numpy', type(my_rdf), flush=True)
        #print('shape', my_rdf.shape, flush=True)
        # Write RDF vector to corresponding json files
        for i,struct_id in enumerate(my_struct_dict):
            my_struct_dict[struct_id].properties['rdf'] = list(np.array(my_rdf[i,:], dtype=str))
            for root_atom in my_struct_dict[struct_id].properties['R']:
                for pair in my_struct_dict[struct_id].properties['R'][root_atom]:
                    tensor = np.array(my_struct_dict[struct_id].properties['R'][root_atom][pair].detach().numpy(), dtype=str)
                    if len(tensor.shape) == 1:
                        tensor_list = list(tensor)
                    else:
                        tensor_list = list(map(list, tensor))
                    my_struct_dict[struct_id].properties['R'][root_atom][pair] = tensor_list
            struct_file = my_struct_files[i]
            if output_dir != 'no_new_output_dir':
                struct_file = os.path.join(output_dir, os.path.basename(struct_file))
            #print('my_struct_dict[struct_id].properties', my_struct_dict[struct_id].properties, flush=True)
            if not save_envs:
                del my_struct_dict[struct_id].properties['R']
            write(struct_file, my_struct_dict[struct_id], file_format='json', overwrite=True)
        # Root rank gather all rdf vectors in order to write pairwise distance matrix
        all_rdf_vecs = split_comm.gather(my_rdf, root=0)
        if split_comm.rank > 0:
            return
        all_rdf_vecs = np.vstack(all_rdf_vecs)
    
    if comm.rank > 0:
        return
    # Calculate pairwise distance matrix
    if normalize_rdf_vectors:
        norm_of_rdf_vecs_columnwise = np.linalg.norm(all_rdf_vecs, axis=0)
        all_rdf_vecs = all_rdf_vecs / norm_of_rdf_vecs_columnwise
    # Could use torch if you want...but no python pairwise distance matrix calculator is as fast as pdist
    #all_rdf_vecs = torch.from_numpy(all_rdf_vecs)
    #dist_mat = torch.norm(all_rdf_vecs[:,None] - all_rdf_vecs, dim=-1)
    #dist_mat = dist_mat.numpy()

    # pdist is a boss! It doesn't need parallelization! Shout out to pdist! Yes, it's 12:30 AM and yes I've been working too long lol
    dist_mat = squareform(pdist(all_rdf_vecs, pdist_distance_type))
    if standardize_distance_matrix:
        dist_mat = dist_mat / np.max(dist_mat)
    #print('dist_mat', dist_mat, flush=True)
    #print('type(dist_mat)', type(dist_mat), flush=True)
    fp = np.memmap(dist_mat_fpath, dtype='float32', mode='write', shape=dist_mat.shape)
    fp[:] = dist_mat[:]
    

def run_rdf_calc(inst, comm):
    sname = 'run_rdf_calc'
    last_section = get_last_active_procedure_name(inst, sname)

    sname_list = [sname, last_section, 'fhi_aims_energy_evaluation', 'harris_approximation_batch', 'pygenarris_structure_generation', 'structure_generation_batch']
    structure_dir = inst.get_inferred(sname, sname_list, ['structure_dir'] + (len(sname_list) - 1) * ['output_dir'], type_='dir', required=True)

    n_D_inter = inst.get_with_default(sname, 'n_D_inter', 12, eval=True)
    init_scheme = inst.get_with_default(sname, 'init_scheme', 'shifted')
    cutoff = inst.get_with_default(sname, 'cutoff', 12, eval=True)
    eta_range = inst.get_with_default(sname, 'eta_range', [0.05,0.5], eval=True)
    Rs_range = inst.get_with_default(sname, 'Rs_range', [0.1,12], eval=True)
    learn_rep = inst.get_boolean(sname, 'learn_rep')
    pdist_distance_type = inst.get_with_default(sname, 'pdist_distance_type', 'euclidean')
    device = inst.get_with_default(sname, 'device', 'cpu')
    dist_mat_fpath = inst.get_with_default(sname, 'dist_mat_fpath', 'rdf_distance_matrix.dat')
    if not dist_mat_fpath.endswith('.dat'):
        raise Exception('Only supporting distance matrices saved as an np memmap with .dat extension')
    output_dir = inst.get_with_default(sname, 'output_dir', 'no_new_output_dir')
    save_envs = inst.get_boolean(sname, 'save_envs')
    normalize_rdf_vectors = inst.get_with_default(sname, 'normalize_rdf_vectors', 'TRUE')
    inst.set(sname, 'normalize_rdf_vectors', normalize_rdf_vectors)
    normalize_rdf_vectors = inst.get_boolean(sname, 'normalize_rdf_vectors')
    standardize_distance_matrix = inst.get_with_default(sname, 'standardize_distance_matrix', 'TRUE')
    inst.set(sname, 'standardize_distance_matrix', standardize_distance_matrix)
    standardize_distance_matrix = inst.get_boolean(sname, 'standardize_distance_matrix')

    acsf_kwargs = {
                        "n_D_inter": n_D_inter,
                        "init_scheme": init_scheme,
                        "cutoff": cutoff,
                        "eta_range": eta_range,
                        "Rs_range": Rs_range,
                        "learn_rep": learn_rep
                     }

    rdf_calc(structure_dir, comm=comm, device=torch.device(device), acsf_kwargs=acsf_kwargs, pdist_distance_type=pdist_distance_type, dist_mat_fpath=dist_mat_fpath, output_dir=output_dir,
                    save_envs=save_envs, normalize_rdf_vectors=normalize_rdf_vectors, standardize_distance_matrix=standardize_distance_matrix)