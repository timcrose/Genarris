# -*- coding: utf-8 -*-

import os

import numpy as np
import torch

from ibslib.io import read,write
from ibslib.analysis.diversity import DiversityAnalysis

# Structure directory to test
struct_dir = "/Users/ibier/Research/Genarris/2019_Publication/Test_Systems/BENZEN/benzene_2mpc_raw_jsons"
# Directory to output pairs of structures
out_dir = "test_structures"

try: s
except:
    s = read(struct_dir)

da = DiversityAnalysis()
descriptor,prop = da.calc(s)
    
    
pairwise_dist = torch.norm(descriptor - descriptor[:,None], dim=-1)
pairwise_dist = pairwise_dist.detach().numpy()

upper_idx = np.triu_indices(pairwise_dist.shape[0],k=1)
upper_result = pairwise_dist[upper_idx]


zero_idx = np.where(upper_result == 0)[0]
zero_idx_pairs = []
for idx in zero_idx:
    zero_idx_pairs.append((upper_idx[0][idx],
                          upper_idx[1][idx]))

for i,pair in enumerate(zero_idx_pairs):
    folder = os.path.join(out_dir, "{}".format(i))
    struct_1 = s[keys[pair[0]]]
    struct_2 = s[keys[pair[1]]]
    temp_dict = {"{}".format(keys[pair[0]]): struct_1,
                 "{}".format(keys[pair[1]]): struct_2}
    write(folder,temp_dict,"geo")