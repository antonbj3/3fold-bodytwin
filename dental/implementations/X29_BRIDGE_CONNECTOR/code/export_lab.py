import shutil, sys
import numpy as np, trimesh
from common import R, read, dump, sha

def run():
    exports = []
    for (design, tag) in [('A9', 'h4_b2.25_s0.18'), ('A9_short', 'h3_b3_s0.18'), ('A16', 'h4_b4_s0.18')]:
        dest = R / 'exports' / design
        dest.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(R / 'raw/fe' / tag / 'bridge.stl', dest / 'bridge.stl')
        a = np.load(R / 'inputs/D1_model.npz')
        for (i, s) in enumerate([-11.0, 11.0], 1):
            T = a['dT']
            F = np.vstack([T[:, [0, 2, 1]], T[:, [0, 1, 3]], T[:, [1, 2, 3]], T[:, [2, 0, 3]]])
            S = np.sort(F, axis=1)
            (_, idx, count) = np.unique(S, axis=0, return_index=True, return_counts=True)
            F = F[idx[count == 1]]
            surface = trimesh.Trimesh(a['dV'] + [s, 0, 0], F, process=True)
            surface.export(dest / f'die{i}.stl')
        files = [dict(path=str(p.relative_to(R)), sha256=sha(p), bytes=p.stat().st_size) for p in sorted(dest.glob('*.stl'))]
        exports.append(dict(design=design, units='mm (STL does not encode units)', coordinate_frame='X1BD1 local toothframe plus ±11mm translations', bridge_surface_gate_p99_mm=read(f'raw/fe/{tag}/FE.json')['surface_p99_error_mm'], files=files, status='CADfixture exported; CAM access/sintercompensation/asbuilt not validated', resolution='PER_POINT'))
    dump('EXPORT_MANIFEST.json', exports)
    return exports
if __name__ == '__main__':
    run()
