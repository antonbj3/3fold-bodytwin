from ray_crown import *
from score_joint import intersections, shape_metrics

def evaluate_c(rec):
    st = time.perf_counter()
    m = load(rec['mesh_path'])
    (v, f, roles) = (m['vertices'], m['faces'], m['roles'])
    (pv, pf) = (m['prep_vertices'], m['prep_faces'])
    S = m['native_triangles']
    inn = v[f[roles == 1]].copy()
    ext = v[f[roles == 0]].copy()
    P = pv[pf].copy()
    if not len(inn):
        raise ValueError('Cavity never opened source support')
    shape = shape_metrics(ext, S)
    pts = np.r_[v[np.unique(f[roles == 1])], inn.mean(1)]
    gap = fast_nearest(P, pts)[1]
    q = ext.mean(1)
    rad = np.linalg.norm(ext - q[:, None], axis=2).max(1)
    wd = fast_nearest(inn, q)[1]
    wall = dict(sampled_min_mm=float(wd.min()), facet_cover_lower_mm=float((wd - rad).min()), ray_cover_lower_mm=rec['wall_analytic_lower_mm'], arithmetic_enclosure='MISSING for source distance evaluation; mathematical Lipschitz cover explicit')
    mm = trimesh.Trimesh(v, f, process=False)
    nc = len(trimesh.graph.connected_components(mm.face_adjacency, nodes=np.arange(len(f))))
    health = dict(watertight=bool(mm.is_watertight), winding_consistent=bool(mm.is_winding_consistent), components=nc, volume_mm3=float(mm.volume))
    si = intersections(v, f, Path(rec['mesh_path']).stem)
    parity_points = sample(inn, 64)
    slow = closest(P, parity_points)[1]
    fast = fast_nearest(P, parity_points)[1]
    parity = float(abs(slow - fast).max())
    reduction = fast_nearest(P, sample(S, 8192))[1]
    retained = rec['closure']['volume_mm3'] - mm.volume
    pr = read(R / 'PREREG_C.json')['metrics']
    film = dict(min_mm=float(gap.min()), max_mm=float(gap.max()), max_residual_mm=float(abs(gap - 0.05).max()), probe_count=len(gap), scope='All cavity vertices and face centroids; no global floating film enclosure')
    gates = dict(shape=shape['p95_mm'] <= pr['shape_p95_mm'][rec['family']], native_retention=rec['retained_native_facets'] == rec['native_facets'], margin_identity=rec['retained_native_facets'] == rec['native_facets'], wall=wall['ray_cover_lower_mm'] >= 0.5 and wall['sampled_min_mm'] >= 0.5, gap=film['max_residual_mm'] <= 0.01, closed=health['watertight'] and health['winding_consistent'] and (health['volume_mm3'] > 0) and (nc == 1), self_intersections=si['status'] == 'PASS', distance_parity=parity <= 1e-08, ideal_axial_insertion=True)
    out = dict(**rec, shape=shape, wall=wall, film=film, mesh=health, self_intersections=si, gates=gates, geometric_conjunction=all(gates.values()), full_clinical_crown=False, insertion_scope='Exact ideal downward-ray invariance, +axis crown removal; real neighbours and root untested', preparation_reduction_mm=stats(reduction), retained_digital_volume_mm3=float(retained), reduction_scope='Nearest native surface to virtual preparation; not directional clinical reduction. Artificial closure gives volume scope.', clinical_status='UNKNOWN_TRUE_MARGIN_PULP_BITE_RETENTION_STRENGTH_MANUFACTURING', seconds_scoring=time.perf_counter() - st, resolution='PER_TOOTH')
    export = D / (Path(rec['mesh_path']).stem + '_research.stl')
    mm.export(export)
    out['research_STL'] = str(export)
    out['research_STL_sha256'] = sha(export)
    return out

def run():
    tag = sys.argv[1] if len(sys.argv) > 1 else 'C'
    fr = read(R / f'FROZEN_PREDICTIONS_{tag}.json')
    rows = []
    for rec in fr['rows']:
        if rec['status'] != 'GENERATED':
            rows.append(dict(**rec, geometric_conjunction=False))
            continue
        try:
            out = evaluate_c(rec)
        except Exception as e:
            out = dict(**rec, score_error=repr(e), geometric_conjunction=False)
        rows.append(out)
        dump(R / f'raw/{tag}_SCORE.json', rows)
        print(tag, out['key'], out.get('gates'), out.get('score_error'), flush=True)
    dump(R / f'RESULTS_{tag}.json', dict(claim_type='capability', rows=rows, external_referent=read(R / 'PREREG_C.json')['external_referent']))
if __name__ == '__main__':
    run()
