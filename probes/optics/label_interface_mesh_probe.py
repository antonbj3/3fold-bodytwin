"""Independent regional volume and edge-incidence checks, including a negative."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
from bodytwin.geometry.label_interfaces_v1 import label_interfaces
from label_interface_mechanism import fixtures, measure


def leg():
    rows = []
    for name, labels in fixtures().items():
        mesh = label_interfaces(labels)
        regions = []
        for region in np.unique(labels):
            if not region:
                continue
            faces = np.concatenate((mesh['faces'][mesh['back'] == region],
                                    mesh['faces'][mesh['front'] == region][:, ::-1]))
            triangles = mesh['vertices'][faces]
            six_volume = np.einsum('ij,ij->i', triangles[:, 0],
                                  np.cross(triangles[:, 1], triangles[:, 2])).sum()
            edges = faces[:, [0, 1, 1, 2, 2, 0]].reshape(-1, 2)
            unique, inverse, count = np.unique(np.sort(edges, axis=1), axis=0,
                                              return_inverse=True, return_counts=True)
            balance = np.zeros(len(unique), np.int64)
            np.add.at(balance, inverse, np.where(edges[:, 0] < edges[:, 1], 1, -1))
            regions.append(dict(region=int(region), six_volume=int(six_volume),
                                expected_six_volume=6 * int(np.count_nonzero(labels == region)),
                                nonmanifold_edges=int(np.count_nonzero(count != 2)),
                                unbalanced_edges=int(np.count_nonzero(balance)),
                                maximum_edge_incidence=int(count.max())))
        rows.append(dict(case=name, triangles=len(mesh['faces']), regions=regions,
                         hashes={k: hashlib.sha256(v.tobytes()).hexdigest() for k, v in mesh.items()}))
    return rows


if __name__ == '__main__':
    a, b = leg(), leg()
    expected = {r['case']: 2 * sum(r['pair_faces'].values()) for r in measure()}
    regions = [r for row in a for r in row['regions']]
    gates = dict(full_repeat=a == b,
                 exact_triangle_count=all(r['triangles'] == expected[r['case']] for r in a),
                 exact_region_volume=all(r['six_volume'] == r['expected_six_volume'] for r in regions),
                 manifold_edges=all(r['nonmanifold_edges'] == 0 for r in regions),
                 directed_edge_balance=all(r['unbalanced_edges'] == 0 for r in regions))
    report = dict(rows=a, gates=gates, scope='Synthetic voxel boundaries; unrestricted manifold claim requires every gate.')
    (ROOT / 'reports/label_interface_mesh_probe.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report))
    raise SystemExit(0 if all(gates.values()) else 2)
