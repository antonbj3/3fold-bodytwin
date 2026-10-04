from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from util import PAYLOAD, read

def attach(t, arrays):
    t = dict(t)
    if t['status'] != 'READY':
        return t
    for key in ['xy', 'faces', 'uv', 'prior', 'ceiling', 'weights']:
        t[key] = arrays[key]
    n = len(t['xy'])
    t['A'] = csr_matrix((arrays['a_data'], arrays['a_indices'], arrays['a_indptr']), shape=(int(arrays['a_rows']), n))
    t['obstacle_b'] = arrays['b'] - t['requirements']['clearance_mm']
    t['preparation_z'] = np.full(n, t['preparation_height_mm'])
    return t

def load_local(t, public):
    if t['status'] != 'READY':
        return dict(t)
    with np.load(Path(public) / t['geometry_file'], allow_pickle=False) as a:
        return attach(t, dict(a))
