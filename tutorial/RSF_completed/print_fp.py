import numpy as np

fpath = 'affinity_matrix1.dat'
fp = np.memmap(fpath, dtype='float32', mode='r')
fp.resize((int(np.sqrt(len(fp))), int(np.sqrt(len(fp)))))
print(fp)
