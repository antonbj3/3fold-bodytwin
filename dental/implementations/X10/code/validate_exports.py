from common import *
from collections import Counter
import math

def inspect(path):
    verts = []
    for line in path.read_text().splitlines():
        if line.startswith('vertex '):
            verts.append(tuple((float(v) for v in line.split()[1:])))
    tris = [verts[k:k + 3] for k in range(0, len(verts), 3)]
    edges = Counter()
    directed = Counter()
    volume = 0.0
    for (a, b, c) in tris:
        for (u, v) in [(a, b), (b, c), (c, a)]:
            edges[tuple(sorted([u, v]))] += 1
            directed[u, v] += 1
        volume += (a[0] * (b[1] * c[2] - b[2] * c[1]) + a[1] * (b[2] * c[0] - b[0] * c[2]) + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6
    return {'triangles': len(tris), 'all_edges_have_two_faces': all((n == 2 for n in edges.values())), 'all_edges_opposite_orientation': all((n == directed[v, u] for ((u, v), n) in directed.items())), 'signed_volume_mm3': volume, 'max_z_mm': max((v[2] for v in verts))}
checks = {}
for item in load('EXPORT_MANIFEST.json')['files']:
    path = Path(item['path'])
    checks[path.name] = inspect(path)
    assert checks[path.name]['all_edges_have_two_faces'] and checks[path.name]['all_edges_opposite_orientation']
    assert checks[path.name]['signed_volume_mm3'] > 0
write('EXPORT_VALIDATION.json', {'mesh_checks': checks, 'facit_scope': 'our own analytic witness geometry, not an independent physical gap measurement', 'contact_residual_check': 'core run_r4 independent interval vertices and nonpenetration', 'finite_wall_overlap': 'contracted lift334um <5mm wall height; margin-plane projected gap only'})
print('STL topology, orientation and positive volume checks pass.')
