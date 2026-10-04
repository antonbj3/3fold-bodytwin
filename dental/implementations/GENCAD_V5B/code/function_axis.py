"""Read-only reuse of v5 exact triangle contact maps on the adapted exterior."""
from common import *
import importlib, sys
_saved = sys.modules['common']
_path = list(sys.path)
try:
    vc = load_module('v5_common_readonly', V5 / 'gencad_bench_v5/common.py')
    sys.modules['common'] = vc
    sys.path.insert(0, str(V5 / 'gencad_bench_v5'))
    VF = load_module('v5_functional_readonly', V5 / 'gencad_bench_v5/functional.py')
finally:
    sys.modules['common'] = _saved
    sys.path[:] = _path
CACHE = {}

def query(row, rec):
    key = rec['key']
    z = load_npz(row['mesh_path'])
    ext = z['vertices'][z['faces'][z['face_roles'] == 0]]
    if key not in CACHE:
        a = load_npz(V4 / 'payload/whole_private' / key / 'reference.npz')
        p = load_npz(rec['prep_path'])
        tri = a['source_triangles']
        tri = tri[tri.mean(1)[:, 2] >= float(p['margin_z'])]
        (xy, faces) = VF.grid(tri)
        ceiling = VF.height(a['antagonist_triangles'], xy, True)
        CACHE[key] = (a, xy, faces, ceiling)
    (a, xy, faces, ceiling) = CACHE[key]
    roof = VF.height(ext, xy)
    gap = ceiling - roof
    fn = {str(s): VF.functional_map(xy, faces, gap + s) for s in [-0.05, 0, 0.05]}
    em = VF.make_mesh(ext)
    prox = {}
    for side in ['mesial', 'distal']:
        pts = VF.sample(a[side + '_triangles'], 256)
        if len(pts):
            d = VF.distance(em, pts)
            prox[side] = dict(unsigned_min_gap_mm=float(d.min()), unsigned_p05_gap_mm=float(np.quantile(d, 0.05)), resolution='PER_SURFACE_REGION', overlap_status='UNKNOWN_UNSIGNED_DISTANCE')
        else:
            prox[side] = dict(status='UNKNOWN_NO_NEIGHBOUR')
    return dict(function=fn, proximal=prox, source_v5_function=rec.get('function'), source_v5_proximal=rec.get('proximal'), scope='Recomputed after adaptation; static geometric contact, not force or clinical adequacy', resolution='PER_SURFACE_REGION')
