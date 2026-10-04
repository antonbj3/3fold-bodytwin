from construct_b import *
from score_a import exact_failure
sys.path.insert(0, str(BASE / 'LANE_X49_DESIGN_GATE/code'))
import design_gate as x49

def tile(x, z):
    v = np.array([[x, 0, z], [x + 1, 0, z], [x + 1, 1, z], [x, 1, z]], float)
    return v[np.array([[0, 1, 2], [0, 2, 3]])]

def run():
    st = time.perf_counter()
    ref = np.r_[tile(0, 0), tile(3, 0)]
    a = np.r_[tile(0, 0.625), tile(3, 0.875)]
    b = np.r_[tile(0, 0.375), tile(3, 1.125)]
    pa = scorer.sample(a, 128)
    pb = scorer.sample(b, 128)
    da = distance(ref, pa)
    db = distance(ref, pb)
    mean_a = float(da.mean())
    mean_b = float(db.mean())
    assert mean_a == mean_b == 0.75
    wa = x49.robust_unsigned(a, mesh(ref), 0.5, [0.0, 0.0])
    wb = x49.robust_unsigned(b, mesh(ref), 0.5, [0.0, 0.0])
    assert wa['status'] == 'PASS' and wb['status'] == 'FAIL'
    suff = dict(identity_error_mm=abs(mean_a - mean_b), summary_mean_mm=[mean_a, mean_b], downstream_min_mm=[float(da.min()), float(db.min())], downstream_difference_mm=float(da.min() - db.min()), decision=[wa['status'], wb['status']], minimal_extension='The minimum local wall separates this pair; retain regional field for general queries.', resolution={'summary': 'PER_TOOTH', 'downstream': 'PER_SURFACE_REGION'}, kind='our_own_fixture', external_check='Published trimesh nearest-facet implementation inside unchanged X49; analytic parallel plane distance')
    from score_a import exact_failure
    exact_good = exact_failure(tile(0, 0)[0], tile(0, 0.75)[0])
    exact_bad = exact_failure(tile(0, 0)[0], tile(0, 0.25)[0])
    assert not exact_good['fails'] and exact_bad['fails']
    (x, f) = scorer.grid(ref, 0.25)
    g = np.full(len(x), 0.05)
    pr = read(BASE / 'PROOF_LANE_FULL_CROWN_R4/PREREG_B.json')['metrics']['v6']
    cg = scorer.contact.gates(scorer.contact.compare(x, f, g, g), pr)
    cb = scorer.contact.gates(scorer.contact.compare(x, f, g + 1, g), pr)
    curve = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0.25], [0, 1, 0]], float)
    mg = curve_error(curve, curve)
    mb = curve_error(curve + [0, 0, 1], curve)
    box = trimesh.creation.box()
    openbox = trimesh.Trimesh(box.vertices, box.faces[:-1], process=False)
    pv = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0]])
    pf = np.array([[0, 1, 2]])
    pts = np.array([[0.25, 0.25, 0.05], [0.25, 0.25, 0.2]])
    (sd, _, _, _) = igl.signed_distance(pts, pv, pf, igl.SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
    controls = dict(wall_X49=dict(positive=wa['status'] == 'PASS', injected_rejected=wb['status'] == 'FAIL'), wall_GenCAD_V2=dict(positive=not exact_good['fails'], injected_rejected=exact_bad['fails']), contact_V6=dict(positive=all(cg.values()), injected_rejected=not all(cb.values())), margin=dict(positive=mg <= 0.025, injected_rejected=mb > 0.025), topology=dict(positive=box.is_watertight, injected_rejected=not openbox.is_watertight), offset=dict(positive=abs(sd[0] - 0.05) <= 0.01, injected_rejected=abs(sd[1] - 0.05) > 0.01), shape=dict(positive=scorer.metrics(ref, ref, 128)['p95_mm'] <= 0.24, injected_rejected=scorer.metrics(ref + [0, 0, 1], ref, 128)['p95_mm'] > 0.24))
    (q, ds, j) = fast_nearest(ref, pts)
    (qt, dt, jt) = closest(ref, pts)
    err = float(np.max(abs(ds - dt)))
    controls['nearest_backend'] = dict(positive=err < 1e-09, injected_rejected=float(np.max(abs(ds + 0.1 - dt))) > 1e-09, error_mm=err)
    assert all((c['positive'] and c['injected_rejected'] for c in controls.values()))
    path = D / 'sufficiency.npz'
    np.savez_compressed(path, reference=ref, state_a=a, state_b=b)
    suff.update(raw_path=path, raw_sha256=sha(path))
    out = dict(controls=controls, sufficiency=suff, wall_reports=[wa, wb], exact_controls=[exact_good, exact_bad], all_pass=True, seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    save(R / 'raw/CONTROLS.json', out)
    print('Controls pass', len(controls), 'summary identity', suff['identity_error_mm'], 'downstream', suff['downstream_difference_mm'])
    return out
if __name__ == '__main__':
    run()
