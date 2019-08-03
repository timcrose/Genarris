from ibslib.analysis.diversity import DiversityAnalysis
from ibslib.io import read,write
import torch, os, glob
import numpy as np
from scipy.spatial.distance import pdist, squareform
from Genarris.core.instruct import get_last_active_procedure_name

def rdf_calc(comm=None, structure_path, device=torch.device("cpu"),
                 acsf_kwargs =
                     {
                        "n_D_inter": 12,
                        "init_scheme": "shifted",
                        "cutoff": 12,
                        "eta_range": [0.05,0.5],
                        "Rs_range": [0.1,12],
                        "learn_rep": False
                     }, pdist_distance_type='cosine', dist_mat_fpath='rdf_distance_matrix.npy',
                     output_dir='no_new_output_dir')
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
        path to file to write distance matrix to. It should end in .npy or .npz

    output_dir: str
        Path of directory to write structures to (will create if it DNE). If 'no_new_output_dir' then
        input structures will be overwritten.

    Return None

    Purpose: calculation RDF vector for every atomic species pairwise combination. Update each json file 
        with its RDF vector and then write the structure pairwise distance matrix to a file for later use.

    Notes:
        1) This function is parallelized because a large number of files and/or a large number of different
            atomic species may benefit from parallelization. However, it has a serial option because in general,
            it is fast.
    '''
    if comm is None or comm.size == 1 or os.path.isfile(structure_path):
        # Serial version
        if comm is not None and comm.size > 1 and os.path.isfile(structure_path) and comm.rank > 0:
            return
        s = read(structure_path)
        da = DiversityAnalysis(device=device, acsf_kwargs=acsf_kwargs)
        all_rdf_vecs,_ = da.calc(s)
    else:
        # Split up the work among ranks
        if not os.path.isdir(structure_path):
            raise Exception('structure_path must be a valid directory. Got:', structure_path)
        struct_files = glob.glob(os.path.join(structure_path, '*.json'))
        num_struct_files = len(struct_files)
        num_structs_per_rank = int(num_struct_files / comm.size)
        if comm.rank == comm.size - 1:
            my_struct_files = struct_files[comm.rank * num_structs_per_rank :]
        else:
            my_struct_files = struct_files[comm.rank * num_structs_per_rank : (comm.rank + 1) * num_structs_per_rank]
        # Read in my json files
        my_struct_dict = {}
        for struct_file in my_struct_files:
            struct = read(struct_file)
            my_struct_dict[struct.struct_id] = struct
        # calcluate RDF
        da = DiversityAnalysis(device=device, acsf_kwargs=acsf_kwargs)
        my_rdf,_ = da.calc(my_struct_dict)
        # Write RDF vector to corresponding json files
        for i,struct_id in enumerate(my_struct_dict):
            my_struct_dict[struct_id]['rdf'] = my_rdf[i,:]
            struct_file = my_struct_files[i]
            if output_dir != 'no_new_output_dir':
                struct_file = os.path.join(output_dir, os.path.basename(struct_file))
            write(struct_file, my_struct_dict[struct_id], file_format='json', overwrite=True)
        # Root rank gather all rdf vectors in order to write pairwise distance matrix
        all_rdf_vecs = comm.gather(my_rdf, root=0)
        if comm.rank > 0:
            return
        all_rdf_vecs = np.vstack(all_rdf_vecs)
    # Calculate pairwise distance matrix
    dist_mat = squareform(pdist(all_rdf_vecs, pdist_distance_type))
    np.save('rdf_distance_matrix.npy', dist_mat)


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
    pdist_distance_type = inst.get_with_default(sname, 'pdist_distance_type', 'cosine')
    device = inst.get_with_default(sname, 'device', 'cpu')
    dist_mat_fpath = inst.get_with_default(sname, 'dist_mat_fpath', 'rdf_distance_matrix.npy')
    output_dir = inst.get_with_default(sname, 'output_dir', 'no_new_output_dir')

    acsf_kwargs = {
                        "n_D_inter": n_D_inter,
                        "init_scheme": init_scheme,
                        "cutoff": cutoff,
                        "eta_range": eta_range,
                        "Rs_range": Rs_range,
                        "learn_rep": learn_rep
                     }

    rdf_calc(comm, structure_dir, device=torch.device(device), acsf_kwargs=acsf_kwargs, pdist_distance_type=pdist_distance_type, dist_mat_fpath=dist_mat_fpath, output_dir=output_dir)