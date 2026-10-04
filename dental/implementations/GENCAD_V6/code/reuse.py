"""Only pure functions from parent lanes. No upstream writes/main calls."""
from common import *
saved = sys.modules.get('common')
path = list(sys.path)
try:
    vc = module('v5_common_readonly', V5 / 'gencad_bench_v5/common.py')
    sys.modules['common'] = vc
    sys.path.insert(0, str(V5 / 'gencad_bench_v5'))
    VF = module('v5_functional_readonly', V5 / 'gencad_bench_v5/functional.py')
finally:
    if saved is not None:
        sys.modules['common'] = saved
    else:
        sys.modules.pop('common', None)
    sys.path[:] = path

def grid(tri, h=0.25):
    v = tri.reshape(-1, 3)
    axes = [np.arange(v[:, k].min(), v[:, k].max() + h / 2, h) for k in range(2)]
    (xx, yy) = np.meshgrid(*axes, indexing='ij')
    i = np.arange(xx.size).reshape(xx.shape)
    faces = np.concatenate([np.stack([i[:-1, :-1], i[1:, :-1], i[:-1, 1:]], -1).reshape(-1, 3), np.stack([i[1:, :-1], i[1:, 1:], i[:-1, 1:]], -1).reshape(-1, 3)])
    return (np.c_[xx.ravel(), yy.ravel()], faces)
