from dental_release.paths import expand as _release_expand
from common import *
import trimesh
from trimesh.ray.ray_triangle import RayMeshIntersector

def brute(tri, p, d):
    e1 = tri[:, 1] - tri[:, 0]
    e2 = tri[:, 2] - tri[:, 0]
    h = np.cross(np.broadcast_to(d, e2.shape), e2)
    det = np.einsum('ij,ij->i', e1, h)
    valid = np.abs(det) > 1e-12
    inv = np.zeros(len(det))
    inv[valid] = 1 / det[valid]
    s = p - tri[:, 0]
    u = np.einsum('ij,ij->i', s, h) * inv
    q = np.cross(s, e1)
    v = q @ d * inv
    t = np.einsum('ij,ij->i', e2, q) * inv
    ok = valid & (u >= -1e-10) & (v >= -1e-10) & (u + v <= 1 + 1e-10) & (t > 1e-07)
    return float(max(0, np.max(t[ok]) if np.any(ok) else 0))
start = time.perf_counter()
inputs = json.loads((ROOT / 'sources/INPUTS.json').read_text())
rows = []
for source in inputs['test_files']:
    case = source['name'].split('/')[-1].split('_')[0]
    inputpath = DATA / (case + '_R1.npz')
    if not inputpath.exists():
        inputpath = pathlib.Path(_release_expand('@DENTAL_WORK_ROOT@/LANE_XBREAK_HUNT_5')) / (case + '_R1.npz')
    z = np.load(inputpath)
    v = z['v']
    f = z['f'][z['lab'] == 1]
    tri = v[f]
    cent = tri.mean(axis=1)
    mesh = trimesh.Trimesh(v, f, process=False)
    rays = RayMeshIntersector(mesh)
    (topo, b) = topology(z['v'], z['f'], z['lab'])
    sub = []
    for (name, d) in [('centroid_axis', z['axis']), ('R1_cone_axis', z['d'])]:
        eps = 1e-07
        dist = np.zeros(len(cent))
        for first in range(0, len(cent), 32):
            current = cent[first:first + 32]
            (loc, ids, face) = rays.intersects_location(current + eps * d, np.broadcast_to(d, current.shape), multiple_hits=True)
            if len(ids):
                np.maximum.at(dist, first + ids, (loc - current[ids]) @ d)
        dist[dist < 2e-07] = 0
        witness = int(np.argmax(dist))
        check = sorted(set([witness] + np.linspace(0, len(cent) - 1, 16, dtype=int).tolist()))
        control = np.array([brute(tri, cent[i], d) for i in check])
        error = float(np.max(np.abs(control - dist[check])))
        sub.append({'name': name, 'direction': d.tolist(), 'rays': len(cent), 'nonzero_rays': int(np.sum(dist > 0)), 'max_lift_source_unit': float(dist.max()), 'p95_lift_source_unit': float(np.quantile(dist, 0.95)), 'area_fraction_over_0_05': float(np.sum(mesh.area_faces[dist > 0.05]) / mesh.area), 'witness_face': witness, 'witness_origin': cent[witness].tolist(), 'witness_envelope_point': (cent[witness] + dist[witness] * d).tolist(), 'control_error': error, 'control_samples': len(check), 'sample_scenario_pass': bool(dist.max() <= 0.05 and topo['disk_gate'])})
        np.savez_compressed(DATA / (case + '_R2_' + name + '.npz'), lift=dist, checked_indices=check, control=control)
    rows.append({'case': case, 'topology_gate': topo['disk_gate'], 'directions': sub, 'physical_gate': 'UNKNOWN'})
    print(case, [(r['name'], r['max_lift_source_unit'], r['control_error']) for r in sub], flush=True)
out = {'id': 'XBREAK5-R2', 'claim_type': 'capability', 'rows': rows, 'sample_pass_n': sum((any((q['sample_scenario_pass'] for q in r['directions'])) for r in rows)), 'validation_gate': 'PASS' if all((q['control_error'] <= 1e-05 for r in rows for q in r['directions'])) else 'FAIL', 'primary_gate': 'PASS_SAMPLED_SCENARIO' if any((q['sample_scenario_pass'] for r in rows for q in r['directions'])) else 'FAIL', 'continuous_gate': 'UNKNOWN_NOT_SAMPLED', 'physical_gate': 'UNKNOWN', 'wall_seconds': time.perf_counter() - start}
write('raw/R2.json', out)
print(json.dumps({k: v for (k, v) in out.items() if k != 'rows'}))
write('CURRENT_WORK_STATE.json', {'lane': _release_expand('XBREAK_HUNT_5'), 'status': 'R2_COMPLETE', 'latest_gate': out['primary_gate'], 'next_operation': 'R3 continuous projected-facet overlaps, no denser unchanged sampling', 'physical_gate': 'UNKNOWN'})
