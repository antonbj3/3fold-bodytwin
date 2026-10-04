import os
os.environ.update(OMP_NUM_THREADS='4', OPENBLAS_NUM_THREADS='4', MKL_NUM_THREADS='4')
import json, hashlib, time
from pathlib import Path
import numpy as np
import trimesh
ROOT = Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def closest(mesh, points):
    rr = [trimesh.proximity.closest_point(mesh, p) for p in np.array_split(points, max(1, int(np.ceil(len(points) / 64))))]
    return tuple((np.concatenate([r[k] for r in rr]) for k in range(3)))

def run():
    st = time.perf_counter()
    rows = []
    for e in json.loads((ROOT / 'inputs/R4_FROZEN_EXPORTS.json').read_text())['exports']:
        assert sha(e['mesh_path']) == e['mesh_sha256']
        a = dict(np.load(e['mesh_path'], allow_pickle=False))
        (v, f, r) = (a['vertices'], a['faces'], a['roles'])
        outer = v[f[r == 0]]
        inner = trimesh.Trimesh(v, f[r == 1], process=False)
        points = outer.mean(1)
        (q, d, j) = closest(inner, points)
        area = np.linalg.norm(np.cross(outer[:, 1] - outer[:, 0], outer[:, 2] - outer[:, 0]), axis=1) / 2
        radius = np.linalg.norm(outer - points[:, None, :], axis=2).max(1)
        dest = ROOT / 'raw' / f"{e['family']}_wall.npz"
        np.savez_compressed(dest, points=points, nearest_inner=q, distance=d, nearest_face=j, area=area, radius=radius, outer=outer, vertices=v, faces=f, roles=r)
        k = int(d.argmin())
        row = {'family': e['family'], 'key': e['key'], 'points': len(points), 'min_centroid_distance_mm': float(d[k]), 'quantiles_mm': np.quantile(d, [0, 0.01, 0.05, 0.5, 0.95, 1]).tolist(), 'min_point_mm': points[k].tolist(), 'nearest_inner_mm': q[k].tolist(), 'min_lipschitz_lower_mm': float((d - radius).min()), 'array_path': str(dest), 'array_sha256': sha(dest), 'mesh_sha256': e['mesh_sha256'], 'resolution': 'PER_POINT', 'not_qualification': True}
        rows.append(row)
        print(json.dumps(row), flush=True)
    result = {'rows': rows, 'seconds': time.perf_counter() - st, 'source_uncertainty': 'UNKNOWN', 'roundoff_enclosure': 'MISSING', 'claim_type': 'information_link'}
    (ROOT / 'raw/INITIAL_GEOMETRY.json').write_text(json.dumps(result, indent=2) + '\n')
if __name__ == '__main__':
    run()
