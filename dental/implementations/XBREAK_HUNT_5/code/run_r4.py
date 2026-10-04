from dental_release.paths import expand as _release_expand
from common import *
from columns import coords
from shapely.geometry import Polygon
from shapely import points as shape_points, distance as shape_distance

def independent(p, ring):
    a = ring
    b = np.roll(ring, -1, axis=0)
    edge = b - a
    t = np.clip(np.einsum('ij,ij->i', p - a, edge) / np.einsum('ij,ij->i', edge, edge), 0, 1)
    dist = float(np.linalg.norm(p - (a + t[:, None] * edge), axis=1).min())
    crosses = (a[:, 1] > p[1]) != (b[:, 1] > p[1])
    den = b[:, 1] - a[:, 1]
    xx = np.full(len(a), np.inf)
    xx[crosses] = (b[crosses, 0] - a[crosses, 0]) * (p[1] - a[crosses, 1]) / den[crosses] + a[crosses, 0]
    inside = np.sum(crosses & (p[0] < xx)) % 2 == 1
    return 0.0 if inside or dist < 1e-12 else dist

def gap(poly, tri):
    active = tri.copy()
    lower = 0.0
    upper = 0.0
    witness = None
    iterations = []
    for depth in range(9):
        cent = active.mean(axis=1)
        r = np.linalg.norm(active - cent[:, None, :], axis=2).max(axis=1)
        dist = shape_distance(shape_points(cent), poly)
        k = int(np.argmax(dist))
        if dist[k] >= lower:
            lower = float(dist[k])
            witness = cent[k].tolist()
        up = dist + r
        clear = up <= 0.12
        if np.any(clear):
            upper = max(upper, float(up[clear].max()))
        iterations.append({'depth': depth, 'triangles': len(active), 'lower': float(dist.max()), 'upper_active': float(up.max()), 'cleared': int(clear.sum())})
        if np.any(dist > 0.12):
            return {'gate': 'FAIL', 'lower_max': lower, 'upper_max': max(upper, float(up.max())), 'witness': witness, 'iterations': iterations}
        bad = active[~clear]
        if len(bad) == 0:
            return {'gate': 'CONDITIONAL_PASS', 'lower_max': lower, 'upper_max': upper, 'witness': witness, 'iterations': iterations}
        if depth == 8:
            return {'gate': 'UNKNOWN', 'lower_max': lower, 'upper_max': max(upper, float(up.max())), 'witness': witness, 'iterations': iterations}
        (a, b, c) = (bad[:, 0], bad[:, 1], bad[:, 2])
        ab = (a + b) / 2
        bc = (b + c) / 2
        ca = (c + a) / 2
        active = np.concatenate([np.stack([a, ab, ca], axis=1), np.stack([ab, b, bc], axis=1), np.stack([ca, bc, c], axis=1), np.stack([ab, bc, ca], axis=1)])
start = time.perf_counter()
r3 = json.loads((ROOT / 'raw/R3.json').read_text())
rows = []
for prior in r3['rows']:
    case = prior['case']
    z = np.load(DATA / (case + '_R1.npz'))
    v = z['v']
    f = z['f'][z['lab'] == 1]
    seq = orientloop(z['boundary'])
    d = np.array(prior['direction'])
    (q, basis) = coords(v, d)
    ring = q[seq, :2]
    poly = Polygon(ring)
    r = {'case': case, 'original_margin_vertices': len(seq), 'simple_projected_margin': bool(poly.is_valid and poly.exterior.is_simple), 'source_coordinate_roundtrip_error': float(np.abs(q @ basis.T - v).max())}
    np.savez_compressed(DATA / (case + '_SOURCE_FINISH_LINE.npz'), curve=v[seq], source_vertex_indices=seq, projected_curve=ring, basis=basis)
    if r['simple_projected_margin']:
        outcome = gap(poly, q[f][:, :, :2])
        cent = q[f].mean(axis=1)[:, :2]
        selected = cent[np.linspace(0, len(cent) - 1, 16, dtype=int)]
        pts = np.vstack([selected, np.array(outcome['witness'])])
        control = np.array([independent(p, ring) for p in pts])
        external = shape_distance(shape_points(pts), poly)
        r['distance_control_error'] = float(np.abs(external - control).max())
        r['opening'] = outcome
    else:
        r['opening'] = {'gate': 'UNKNOWN_SELF_INTERSECTING_PROJECTED_MARGIN'}
        r['distance_control_error'] = None
    r['inherited_cavity_gate'] = prior['grids'][-1]['gap']['qualified_gate']
    r['linked_gate'] = 'CONDITIONAL_PASS' if r['opening']['gate'] == 'CONDITIONAL_PASS' and r['inherited_cavity_gate'] == 'CONDITIONAL_PASS' and (r['distance_control_error'] <= 1e-05) else 'FAIL_OR_UNKNOWN'
    r['physical_gate'] = 'UNKNOWN'
    rows.append(r)
out = {'id': 'XBREAK5-R4', 'claim_type': 'capability', 'rows': rows, 'primary_gate': 'CONDITIONAL_PASS' if any((r['linked_gate'] == 'CONDITIONAL_PASS' for r in rows)) else 'FAIL_OR_UNKNOWN', 'physical_gate': 'UNKNOWN', 'wall_seconds': time.perf_counter() - start}
write('raw/R4.json', out)
print(json.dumps(out))
write('CURRENT_WORK_STATE.json', {'lane': _release_expand('XBREAK_HUNT_5'), 'status': 'R4_COMPLETE', 'latest_gate': out['primary_gate'], 'next_operation': 'Consumer guards and independent replay, then same-object measured scale/margin/intaglio acquisition', 'physical_gate': 'UNKNOWN'})
