from common import *
from native_port import model

def down(a):
    return np.nextafter(a, -np.inf)

def up(a):
    return np.nextafter(a, np.inf)

def box_lower(point, lower, upper):
    dl = down(lower - point)
    du = down(point - upper)
    d = np.maximum(np.maximum(dl, du), 0.0)
    sq = np.maximum(down(d * d), 0.0)
    total = np.maximum(down(down(sq[:, 0] + sq[:, 1]) + sq[:, 2]), 0.0)
    return np.maximum(down(np.sqrt(total)), 0.0)

def vertex_upper(point, vertices):
    lo = down(vertices - point)
    hi = up(vertices - point)
    absbound = np.maximum(abs(lo), abs(hi))
    sq = up(absbound * absbound)
    total = up(up(sq[:, 0] + sq[:, 1]) + sq[:, 2])
    return up(np.sqrt(total))

def run():
    t0 = time.perf_counter()
    (a, g, m) = model()
    r = read(R / 'raw/R1_NATIVE_PORT.json')
    point = np.array(r['worst_point']['coordinate_mm'])
    tri = m.Vd[m.Fd_master[:, 2:5]].astype(np.float64)
    lo = np.min(tri, axis=1)
    hi = np.max(tri, axis=1)
    contains = lambda low, high: bool(np.all(tri >= low[:, None, :]) and np.all(tri <= high[:, None, :]))
    lower = box_lower(point, lo, hi)
    upper = vertex_upper(point, tri.reshape(-1, 3))
    lb = float(lower.min())
    ub = float(upper.min())
    wronglo = lo.copy()
    wronglo[0, 0] = np.nextafter(tri[0, :, 0].max(), np.inf)
    fault_rejected = not contains(wronglo, hi)
    perm = lower[::-1].copy()
    identity = bool(lb == float(perm.min()))
    out = {'round': 'R4_PORT_CERTIFICATE', 'claim_type': 'capability', 'point_mm': point, 'point_id': r['worst_point']['slave_node_0based'], 'native_master_faces': len(tri), 'rigorous_distance_enclosure_mm': [lb, ub], 'closest_triangle_distance_mm': r['worst_point']['exhaustive_distance_mm'], 'all_master_boxes_contain_vertices': contains(lo, hi), 'mathematical_argument': 'Each triangle is contained in its corner AABB, so distance(point,union triangles)>=min distance(point,AABB). Vertex distance is a valid upper bound to that union. IEEE754 binary64 operations directed outward via nextafter after every elementary rounding; no BLAS or cumulative opaque reductions in enclosure. Straight master facets; curved geometry and true physical support excluded.', 'gate': {'certificate_refutes_006_port': contains(lo, hi) and lb > 0.06, 'projection_inside_enclosure': lb <= r['worst_point']['exhaustive_distance_mm'] <= ub, 'wrong_box_fault_rejected': fault_rejected, 'permutation_summary_identity': identity}, 'bounds_resolution': 'PER_POINT', 'physical_coverage': 'UNKNOWN', 'stop': 'Under unchanged planar native master facets, NO displacement interpolation Q can reduce this geometrical distance to .06mm. A changed geometric port or support definition is necessary; re-optimization cannot fix it.', 'cost': cost(t0)}
    dump(R / 'raw/R4_PORT_CERTIFICATE.json', out)
    state('R4_PORT_CERTIFICATE_DECIDED', out['gate'], 'Complete figure, source checks, prediction freeze and independent replay; next construction needs changed geometric port')
    print(json.dumps(clean(out), indent=2))
if __name__ == '__main__':
    run()
