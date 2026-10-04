from local import *
from check_intersections import intersect
import geometry as geometry_r4
sys.path.insert(0, str(BASE / 'LANE_X49_DESIGN_GATE/code'))
import design_gate as x49
from export import three_mf, load3mf

def ascii_stl(path, tri):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w') as h:
        h.write('solid R6_nominal_research_mm\n')
        n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        n /= np.maximum(np.linalg.norm(n, axis=1)[:, None], 1e-300)
        for (normal, t) in zip(n, tri):
            h.write('facet normal ' + ' '.join((format(float(x), '.17g') for x in normal)) + '\nouter loop\n')
            for p in t:
                h.write('vertex ' + ' '.join((format(float(x), '.17g') for x in p)) + '\n')
            h.write('endloop\nendfacet\n')
        h.write('endsolid R6_nominal_research_mm\n')

def x49_report(m, p, key, full=True):
    folder = D6 / 'diagnostic_stl' / key
    folder.mkdir(parents=True, exist_ok=True)
    tris = dict(crown=m['vertices'][m['faces']], prep=m['prep_triangles'], antagonist=p['antagonist'])
    paths = {}
    for (name, tri) in tris.items():
        path = folder / (name + '.stl')
        ascii_stl(path, tri)
        paths[name] = str(path)
    data = {k: x49.G.load(path, 1.0) for (k, path) in paths.items()}
    health = x49.mesh_health(data)
    contract = dict(units='mm', common_frame_confirmed=False, physical_geometry_status='virtual preparation; unverified raw bite', ifu_profile='katana-ht', indication='anterior' if 'anterior' in key else 'posterior', input_sha256={k: sha(path) for (k, path) in paths.items()}, surface_error_mm={k: 0.0 for k in paths}, error_source='Encoded nominal geometry only; scanner, registration and process uncertainty UNKNOWN', regions={'exterior': np.flatnonzero(m['roles'] == 0), 'intaglio': np.flatnonzero(m['roles'] == 1), 'preparation': np.arange(len(m['prep_triangles']))})
    save(folder / 'contract.json', contract)
    report = x49.check(paths, contract=str(folder / 'contract.json'), units='mm') if full else dict(rules=dict(mesh_health=health))
    save(folder / 'REPORT.json', report)
    return dict(report_path=folder / 'REPORT.json', report_sha256=sha(folder / 'REPORT.json'), rules={k: v['status'] for (k, v) in report['rules'].items()}, crown_mesh_diagnostics=health['diagnostics'].get('crown'), input_diagnostics=health['diagnostics'], paths=paths, scope='X49 mesh check does not assess self intersections; separate CGAL check mandatory')

def score_mesh(rec, rr, tag):
    m = npz(rec['mesh_path'])
    assert sha(rec['mesh_path']) == rec['mesh_sha256']
    (v, f, roles) = (m['vertices'], m['faces'], m['roles'])
    p = npz(rr['public_path'])
    target = npz(rr['private_path'])['target']
    orig = next((a for (a, b) in inputs() if a['key'] == rec['key']))
    src = npz(orig['mesh_path'])
    (pv, pf) = compact(src['vertices'], src['faces'][src['roles'] == 1])
    rim = pv[loops(trimesh.Trimesh(pv, pf, process=False))[0]]
    ext = v[f[roles == 0]]
    inner = v[f[roles == 1]]
    prep = m['prep_triangles']
    geometry_r4.closest = fast_nearest
    scorer.closest = fast_nearest
    shape = scorer.metrics(ext, target)
    wall = scorer.wall(ext, inner)
    probes = np.r_[v[np.unique(f[roles == 1])], inner.mean(1), sample(inner, 8192)]
    (q, d, j) = fast_nearest(prep, probes)
    gres = abs(d - boundary_gap(q, rim))
    qp = np.r_[pv, sample(prep, 8192)]
    (_, coverage, _) = fast_nearest(inner, qp)
    parity_points = sample(inner, 128)
    (_, fast, _) = fast_nearest(prep, parity_points)
    (_, slow, _) = closest(prep, parity_points)
    parity = float(abs(fast - slow).max())
    mm = trimesh.Trimesh(v, f, process=False)
    om = trimesh.Trimesh(v, f[roles == 0], process=False)
    ll = loops(om)
    margin = curve_error(v[ll[0]], p['margin_curve']) if len(ll) == 1 else None
    ident = cKDTree(v[np.unique(f[roles == 0])]).query(p['margin_curve'], workers=1)[0]
    si = intersect(v, f, rec['key'] + '_' + tag + '_whole')
    pairs = np.array(si.pop('pairs', []), int).reshape(-1, 2)
    dest = D6 / (rec['key'] + '_' + tag + '_intersections.npz')
    np.savez_compressed(dest, pairs=pairs)
    si.update(pairs_path=dest, pairs_sha256=sha(dest), sample_pairs=pairs[:10])
    xr = x49_report(m, p, rec['key'] + '_' + tag)
    gates = dict(shape=shape['p95_mm'] <= read(R6 / 'PREREG_A.json')['metrics']['shape_p95_mm'][rec['family']], margin=margin is not None and margin <= 0.025, margin_identity=bool((ident == 0).all()), closed=bool(mm.is_watertight and mm.is_winding_consistent and (mm.volume > 0) and (len(trimesh.graph.connected_components(mm.face_adjacency, nodes=np.arange(len(f)))) == 1)), wall=wall['continuous_lower_mm'] >= 0.5, gap=bool(gres.max() <= 0.01), coverage=bool(coverage.max() <= 0.06), self_intersections=si['status'] == 'PASS', x49_mesh=xr['rules']['mesh_health'] == 'PASS')
    return dict(rec, shape=shape, wall=wall, gap=dict(probe_count=len(probes), minimum_mm=float(d.min()), maximum_mm=float(d.max()), residual_max_mm=float(gres.max()), scope='All vertices and facet centroids plus8192 area-centroid probes; no continuous gap certificate', resolution='PER_POINT'), coverage=dict(probe_count=len(qp), maximum_mm=float(coverage.max()), resolution='PER_POINT', scope='Every source vertex plus8192 area-centroid probes; no continuous coverage certificate'), margin_mm=margin, margin_identity_error_mm=float(ident.max()), distance_parity_max_mm=parity, self_intersection=si, x49=xr, gates=gates, technical_export_pass=all(gates.values()), physical_status='UNKNOWN_NO_MEASURED_PREPARATION_REGISTERED_BITE_OR_SEATED_FILM')

def run(tag='F_FIRST'):
    st = time.perf_counter()
    fr = read(R6 / f'FROZEN_PREDICTIONS_{tag}.json')
    idx = {b['key']: b for (a, b) in inputs()}
    rows = []
    for rec in fr['rows']:
        t = time.perf_counter()
        if rec['status'] != 'GENERATED':
            r = dict(rec, technical_export_pass=False)
        else:
            try:
                r = score_mesh(rec, idx[rec['key']], tag)
            except Exception as e:
                r = dict(rec, score_error=repr(e), technical_export_pass=False)
                print('SCORE_ERROR', e, flush=True)
        r['score_seconds'] = time.perf_counter() - t
        rows.append(r)
        save(R6 / f'raw/{tag}_SCORE.json', rows)
        print(r['key'], r.get('gates'), r.get('gap'), flush=True)
    out = dict(claim_type='capability', rows=rows, seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    save(R6 / f'RESULTS_{tag}.json', out)
    state6(tag + '_SCORED', dict(requested=len(rows), complete=sum((r['technical_export_pass'] for r in rows))), 'Preserve gate failures; next construction attacks joint wall and topology')
if __name__ == '__main__':
    run(sys.argv[1] if len(sys.argv) > 1 else 'F_FIRST')
