"""Independent trimesh triangle closest-point recomputation on actual cut faces."""
import json, time
import numpy as np
import trimesh
from identity import ROOT, write

def main():
    started = time.monotonic()
    checks = []
    rows = json.loads((ROOT / 'raw/R1_ROWS.json').read_text())
    for row in rows:
        if row.get('status') != 'MEASURED_DIGITAL_LABEL_GEOMETRY':
            continue
        a = np.load(row['spatial_export'])
        h = float(a['spacing'][0])
        c = a['cut_face_centers_zyx_voxel']
        e = a['cut_face_extents_voxel']
        triangles = []
        for axis in range(3):
            use = e[:, axis] == 0
            centers = c[use]
            u = (axis + 1) % 3
            v = (axis + 2) % 3
            corners = []
            for (s, t) in [(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)]:
                q = centers.copy()
                q[:, u] += s
                q[:, v] += t
                corners.append(q * h)
            quad = np.stack(corners, axis=1)
            triangles.extend([quad[:, [0, 1, 2]], quad[:, [0, 2, 3]]])
        tri = np.concatenate(triangles)
        q = a['pulp_local_zyx_voxel']
        indices = np.unique(np.linspace(0, len(q) - 1, min(10, len(q)), dtype=int))
        indices = np.unique(np.r_[indices, np.argmin(a['exact_point_cut_distance_mm'])])
        values = []
        for i in indices:
            point = q[i] * h
            closest = trimesh.triangles.closest_point(tri, np.broadcast_to(point, (len(tri), 3)))
            values.append(float(np.min(np.linalg.norm(closest - point, axis=1))))
        ref = a['exact_point_cut_distance_mm'][indices]
        err = float(np.max(np.abs(np.array(values) - ref)))
        poison = float(np.max(np.abs(np.array(values) - (ref + 1.0)))) > 1e-06
        checks.append({'fdi': row['fdi'], 'operation': row['operation'], 'independent_points': len(indices), 'max_error_mm': err, 'valid_pass': err <= 1e-06, 'distance_plus1_rejected': poison})
    assert all((r['valid_pass'] and r['distance_plus1_rejected'] for r in checks))
    write(ROOT / 'raw/INDEPENDENT_FACE_CHECKS.json', checks)
    write(ROOT / 'raw/INDEPENDENT_FACE_OUTCOME.json', {'rows': len(checks), 'points': sum((r['independent_points'] for r in checks)), 'max_error_mm': max((r['max_error_mm'] for r in checks)), 'all_plus1_faults_rejected': True, 'implementation': 'trimesh.triangles.closest_point actual cut-face triangles, everyR1solid; sampled pulp points plus each controlling point', 'kind': 'independent numerical implementation; not independent scientist or physical measurement', 'wall_s': time.monotonic() - started})
    print('Independent triangle distance checks PASS:', len(checks), 'actual solids')
if __name__ == '__main__':
    main()
