from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from .common import read, sha

def load_task(path):
    path = Path(path)
    t = read(path)
    if t.get('status') != 'READY':
        return t
    name = t['geometry_file']
    if Path(name).name != name:
        raise ValueError('task path traversal')
    p = path.parent / 'scenes' / name
    if sha(p) != t['geometry_sha256']:
        raise ValueError('public geometry hash')
    with np.load(p, allow_pickle=False) as a:
        for key in ['xy', 'faces', 'prior', 'ceiling', 'weights', 'uv', 'neighbor_points', 'antagonist']:
            t[key] = a[key]
        t['A'] = csr_matrix((a['a_data'], a['a_indices'], a['a_indptr']), shape=(int(a['a_rows']), len(t['xy'])))
        t['obstacle_b'] = a['b'] - t['requirements']['clearance_mm']
    t['preparation_z'] = np.full(len(t['xy']), t['preparation_height_mm'])
    return t
