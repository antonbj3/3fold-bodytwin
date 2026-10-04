from joint_solid import *
from geometry import curve_error
import warnings
from scipy.optimize import linprog

def region_classification(v, f, S, A, r, prep=False):
    pts = v[f].mean(1)
    ds = fast_nearest(S, pts)[1]
    da = fast_nearest(A, pts)[1]
    if prep:
        return np.where(da < abs(ds - r - (0.025 + 0.025 * np.clip(da, 0, 1))), 2, 1)
    return np.argmin(np.c_[ds, abs(ds - r), da], axis=1)

def intersections(v, f, key):
    p = D / (key + '_check.mesh')
    meshwrite(p, v, f)
    prog = D.parent / 'PROOF_LANE_FULL_CROWN_R6/intersections'
    x = subprocess.run([str(prog), str(p)], capture_output=True, text=True)
    if x.returncode:
        return dict(status='UNKNOWN', returncode=x.returncode, stderr=x.stderr[:300])
    y = json.loads(x.stdout)
    count = y['count']
    y['pairs'] = y['pairs'][:20]
    y.update(status='PASS' if count == 0 else 'FAIL', backend='CGAL exact predicates', mesh_sha256=sha(p))
    return y

def shape_metrics(a, b):
    da = fast_nearest(b, sample(a, 8192))[1]
    db = fast_nearest(a, sample(b, 8192))[1]
    return dict(p95_mm=float(max(np.quantile(da, 0.95), np.quantile(db, 0.95))), sampled_max_mm=float(max(da.max(), db.max())))

def insertion(v, f):
    tri = v[f]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    n /= np.maximum(np.linalg.norm(n, axis=1)[:, None], 1e-30)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        x = linprog([0, 0, 0, -1], A_ub=np.c_[-n, np.ones(len(n))], b_ub=np.zeros(len(n)), bounds=[(-1, 1)] * 3 + [(None, None)], method='highs', options={'threads': 1})
    return dict(status='NUMERICAL_STRICT_DIRECTION' if x.success and x.x[3] > 1e-06 else 'NO_STRICT_DIRECTION_OR_UNKNOWN', angular_margin=float(x.x[3]) if x.success else None, direction=x.x[:3].tolist() if x.success else None, scope='Necessary local infinitesimal normal test; a strict result alone does not prove whole-arch insertion')

def evaluate(rec):
    t0 = time.perf_counter()
    m = load(rec['mesh_path'])
    (v, f) = (m['vertices'], m['faces'])
    (pv, pf) = (m['prep_vertices'], m['prep_faces'])
    S = m['native_triangles']
    A = m['closure_triangles']
    r = 0.6
    roles = region_classification(v, f, S, A, r)
    proles = region_classification(pv, pf, S, A, r, True)
    ext = v[f[roles == 0]].copy()
    inn = v[f[roles == 1]].copy()
    prep = pv[pf[proles == 1]].copy()
    pr = read(R / 'PREREG_B.json')['metrics']
    shape = shape_metrics(ext, S)
    center = ext.mean(1)
    rad = np.linalg.norm(ext - center[:, None], axis=2).max(1)
    ds = fast_nearest(inn, center)[1]
    wall = dict(sampled_min_mm=float(ds.min()), continuous_lower_mm=float((ds - rad).min()), scope='1-Lipschitz geometric bound, floating and source enclosure MISSING')
    pp = np.r_[v[np.unique(f[roles == 1])], inn.mean(1)]
    (q, d, _) = fast_nearest(prep, pp)
    desired = 0.025 + 0.025 * np.clip(fast_nearest(A, q)[1], 0, 1)
    gp = dict(actual_mm=stats(d), residual_mm=stats(abs(d - desired)), probe_count=len(pp), scope='Every cap vertex and facet centroid; continuous film enclosure MISSING')
    pp2 = np.r_[pv[np.unique(pf[proles == 1])], prep.mean(1)]
    cover = fast_nearest(inn, pp2)[1]
    coverage = dict(maximum_mm=float(cover.max()), probe_count=len(pp2))
    fullprep = trimesh.Trimesh(pv, pf, process=False)
    mm = trimesh.Trimesh(v, f, process=False)
    health = dict(watertight=bool(mm.is_watertight), winding_consistent=bool(mm.is_winding_consistent), volume_mm3=float(mm.volume), components=len(trimesh.graph.connected_components(mm.face_adjacency, nodes=np.arange(len(f)))))
    (native_v, native_f) = compact(S.reshape(-1, 3), np.arange(len(S) * 3).reshape(-1, 3))
    nm = trimesh.Trimesh(native_v, native_f, process=True)
    nl = loops(nm)
    source_curve = nm.vertices[nl[0]]
    outer = trimesh.Trimesh(v, f[roles == 0], process=False)
    try:
        ll = loops(outer)
        margin = dict(loop_count=len(ll), sampled_max_mm=curve_error(v[ll[0]], source_curve) if len(ll) == 1 else None)
    except Exception as e:
        margin = dict(error=repr(e), sampled_max_mm=None)
    ins = insertion(pv, pf[proles == 1])
    si = intersections(v, f, rec['key'] + '_B')
    points = sample(inn, 64)
    slow = closest(prep, points)[1]
    fast = fast_nearest(prep, points)[1]
    parity = float(abs(slow - fast).max())
    reductions = fast_nearest(prep, sample(S, 8192))[1]
    gates = dict(shape=shape['p95_mm'] <= pr['shape_p95_mm'][rec['family']], wall=wall['continuous_lower_mm'] >= 0.5, gap=gp['residual_mm']['max'] <= 0.01, coverage=coverage['maximum_mm'] <= 0.06, margin=margin['sampled_max_mm'] is not None and margin['sampled_max_mm'] <= 0.025, closed=health['watertight'] and health['winding_consistent'] and (health['volume_mm3'] > 0) and (health['components'] == 1), self_intersections=si['status'] == 'PASS', distance_parity=parity <= 1e-08)
    annotated = D / (rec['key'] + '_B_regions.npz')
    np.savez_compressed(annotated, vertices=v, faces=f, roles=roles, prep_vertices=pv, prep_faces=pf, prep_roles=proles)
    return dict(**rec, shape=shape, wall=wall, gap=gp, coverage=coverage, margin=margin, mesh=health, self_intersection=si, insertion=ins, preparation_reduction_mm=stats(reductions), reduction_scope='Nearest source-exterior to new prep-cap distance; not axial/occlusal directional reduction', gates=gates, geometric_conjunction=all(gates.values()), full_crown=False, physical_status='UNKNOWN_NO_TRUE_FINISH_LINE_PULP_BITE_OR_LAB', role_classification=dict(method='Nearest active surface, native/closure distinct', generator_label_disagreements=int(np.sum(roles != m['roles'])), canonical_path=str(annotated), canonical_sha256=sha(annotated)), seconds_scoring=time.perf_counter() - t0, resolution='PER_TOOTH')

def run():
    rows = []
    st = time.perf_counter()
    for rec in read(R / 'FROZEN_PREDICTIONS_B.json')['rows']:
        if rec['status'] != 'GENERATED':
            rows.append(rec)
            continue
        try:
            out = evaluate(rec)
        except Exception as e:
            out = dict(**rec, score_error=repr(e), geometric_conjunction=False, full_crown=False)
        rows.append(out)
        dump(R / 'raw/B_SCORE.json', rows)
        print(out['key'], out.get('gates'), out.get('score_error'), flush=True)
    dump(R / 'RESULTS_B.json', dict(claim_type='capability', rows=rows, seconds=time.perf_counter() - st, external_referent=read(R / 'PREREG_B.json')['external_referent']))
if __name__ == '__main__':
    run()
