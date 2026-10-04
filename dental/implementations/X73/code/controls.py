"""Recompute semantics and source bindings before report consumption; poison probes."""
from dental_release.paths import expand as _release_expand
from common import *
import copy, math, resource
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
sys.path.insert(0, str(X8))
from full_geometry import voxel_cylinder_bracket
sys.path.insert(0, str(X58 / 'code'))
from sections import ray_exit
from calibrate import quantile95, guide_budget

def independent_gate(x):
    ids = x.get('annotator_ids', [])
    return x.get('exact_image_equal') is True and x.get('blinded') is True and (len(set(ids)) >= 2) and (len(set(ids)) == len(ids)) and (x.get('preconsensus_masks_retained') is True)

def validate(r, p, s, certsha, checkbytes=True, byid=None):
    assert r['patient'] == p['patient'] and r['case'] == p['case'] == s['case'] and (r['fdi'] == s['fdi'])
    assert r['source_bindings'] == p['source_bindings'] and r['point_artifact'] == p['point_artifact']
    assert r['reference_certificate_sha256'] == certsha
    assert r['pose'] == {k: s[k] for k in ['entry_zyx_mm', 'axis_zyx', 'length_mm', 'radius_mm']}
    if checkbytes:
        assert sha(r['point_artifact']['path']) == r['point_artifact']['sha256']
    if byid is not None:
        expected_aliases = [ap for ap in byid.values() if ap['patient'] != p['patient'] and ap['image_identity']['tf1']['sha256'] == p['image_identity']['tf1']['sha256']]
        actual = r.get('alias_source_bindings', [])
        assert {a['patient'] for a in actual} == {a['patient'] for a in expected_aliases}
        for a in actual:
            ap = byid[a['patient']]
            assert a['source_bindings'] == ap['source_bindings'] and a['point_artifact'] == ap['point_artifact']
            assert a['image_sha256'] == p['image_identity']['tf1']['sha256'] == ap['image_identity']['tf1']['sha256']
            if checkbytes:
                assert sha(a['point_artifact']['path']) == a['point_artifact']['sha256']
        expected_keys = {x['dataset'] for x in p['source_bindings']} | {a['patient'] + '::' + x['dataset'] for a in expected_aliases for x in a['source_bindings']}
        assert set(r['distances']) == expected_keys
    for (key, d) in r['distances'].items():
        assert math.isfinite(d['lower_mm']) and math.isfinite(d['upper_mm']) and (0 <= d['lower_mm'] <= d['upper_mm'])
        assert d['gap_mm'] == d['upper_mm'] - d['lower_mm'] and d['gap_mm'] <= 0.0001
        assert r['classes'][key] == classify(d['lower_mm'], d['upper_mm'])
    lo = min((d['lower_mm'] for d in r['distances'].values()))
    hi = min((d['upper_mm'] for d in r['distances'].values()))
    assert r['union'] == dict(lower_mm=lo, upper_mm=hi, class_2mm=classify(lo, hi))
    ref = r['distances']['tf2']
    assert r['loss_lower_mm'] == max(0.0, ref['lower_mm'] - hi) and r['loss_upper_mm'] == max(0.0, ref['upper_mm'] - lo)
    assert max((abs(ref[k] - s[k]) for k in ['lower_mm', 'upper_mm'])) <= 0.0001
    return True

def rejects(f):
    try:
        f()
    except (AssertionError, ValueError, KeyError):
        return True
    return False

def main():
    check_frozen()
    start = time.perf_counter()
    pat = json.loads((ROOT / 'raw/PATIENTS.json').read_text())
    sites = json.loads((ROOT / 'raw/LINEAGE_SITES.json').read_text())
    byid = {p['patient']: p for p in pat}
    certpath = Path(_release_expand('@DENTAL_WORK_ROOT@/X8-guide-nerve-risk/R5_FULL_GEOMETRY.jsonl'))
    certsha = sha(certpath)
    cs = {(s['case'], s['fdi']): s for s in map(json.loads, certpath.read_text().splitlines())}
    for r in sites:
        validate(r, byid[r['patient']], cs[r['case'], r['fdi']], certsha, byid=byid)
    r = copy.deepcopy(sites[0])
    p = byid[r['patient']]
    s = cs[r['case'], r['fdi']]
    mutations = {}
    for key in r['classes']:
        bad = copy.deepcopy(r)
        bad['classes'][key] = 'BELOW' if bad['classes'][key] == 'ABOVE' else 'ABOVE'
        mutations['false_cached_' + key + '_class_rejected'] = rejects(lambda : validate(bad, p, s, certsha))
    for (name, path, val) in [('union_distance', ('union', 'lower_mm'), r['union']['lower_mm'] + 0.1), ('source_hash', ('point_artifact', 'sha256'), '0' * 64), ('mask_member_hash', ('source_bindings', 0, 'member_sha256'), '0' * 64), ('pose', ('pose', 'radius_mm'), r['pose']['radius_mm'] + 0.1), ('loss', ('loss_upper_mm',), r['loss_upper_mm'] + 0.1), ('reference', ('reference_certificate_sha256',), '0' * 64)]:
        bad = copy.deepcopy(r)
        obj = bad
        for key in path[:-1]:
            obj = obj[key]
        obj[path[-1]] = val
        mutations[name + '_mutation_rejected'] = rejects(lambda : validate(bad, p, s, certsha))
    chosen = next((p for p in pat if len(p['source_bindings']) == 3))
    with np.load(chosen['point_artifact']['path']) as dat:
        x = dat['tf1_voxels_zyx'].astype(int)
        y = dat['tf2_voxels_zyx'].astype(int)
        low = np.minimum(x.min(0), y.min(0)) - 1
        hi = np.maximum(x.max(0), y.max(0)) + 2
        shape = hi - low
        m = np.zeros(shape, bool)
        m[tuple((y - low).T)] = True
        edt = ndi.distance_transform_edt(~m, sampling=SP)
        e = edt[tuple((x - low).T)]
        kd = cKDTree(y * SP).query(x * SP, workers=1)[0]
        err = float(np.max(np.abs(e - kd)))
        brute = cdist(x[:257] * SP, y * SP).min(1)
        be = float(np.max(np.abs(brute - kd[:257])))
        mutations['real_distance_0p1mm_poison_rejected'] = np.max(np.abs(e - (kd + 0.1))) > 1e-10
    cube = np.zeros((15, 15, 15), bool)
    cube[3:12, 4:11, 5:10] = True
    rng = np.random.default_rng(73)
    dirs = rng.normal(size=(127, 3))
    dirs /= np.linalg.norm(dirs, axis=1)[:, None]
    origin = np.tile(np.array([7, 7, 7]) * SP, (127, 1))
    got = ray_exit(cube, origin, dirs, max_length=20)
    lo = (np.array([3, 4, 5]) - 0.5) * SP
    hi = (np.array([11, 10, 9]) + 0.5) * SP
    exact = np.min(np.where(dirs > 0, (hi - origin) / dirs, (lo - origin) / dirs), axis=1)
    re = float(np.max(np.abs(got - exact)))
    mutations['ray_exit_0p3mm_poison_rejected'] = float(np.max(np.abs(got + 0.3 - exact))) > 1e-10
    cr = voxel_cylinder_bracket(np.array([[2.0, 0.0, 4.0]]), np.full(3, 0.3), np.zeros(3), np.array([1.0, 0.0, 0.0]), 5.0, 1.0)
    ce = max((abs(cr[k] - 2.85) for k in ['lower_mm', 'upper_mm']))
    mutations['cylinder_0p1mm_poison_rejected'] = abs(cr['upper_mm'] + 0.1 - 2.85) > 1e-10
    ref = np.array([[z, 0, 0] for z in range(-16, 17)], float)
    aa = ref + [0, 0, -1]
    bb = ref + [0, 0, 1]
    q = np.array([[0, 0, 2.5]])
    da = cKDTree(ref).query(aa)[0]
    db = cKDTree(ref).query(bb)[0]
    sa = np.array([len(aa), 0, np.median(da), np.quantile(da, 0.95), da.max()])
    sb = np.array([len(bb), 0, np.median(db), np.quantile(db, 0.95), db.max()])
    qa = float(cdist(q, aa).min())
    qb = float(cdist(q, bb).min())
    suff = dict(summary=['point_count', 'Dice_to_reference', 'median', 'P95', 'max'], state_a=sa.tolist(), state_b=sb.tolist(), identity_error=float(np.max(np.abs(sa - sb))), entire_absolute_distribution_equal=bool(np.array_equal(np.sort(da), np.sort(db))), downstream_a_mm=qa, downstream_b_mm=qb, downstream_difference_mm=qa - qb, decision_2mm_a=classify(qa, qa), decision_2mm_b=classify(qb, qb), resolution='PER_POINT', time_scale='SIMULTANEOUS', smallest_extension='Signed clearance decrement for this fixed query; wall position/direction plus support for arbitrary new queries', external_referent=dict(kind='our_own_fixture', locator='code/controls.py mirrored canal columns', compared_quantity='summary identity versus point clearance', refutes_us=True))
    assert suff['identity_error'] == 0 and suff['entire_absolute_distribution_equal'] and (qa - qb == 2)
    rm = np.zeros((9, 9, 13), bool)
    rm[3:6, 2:7, 2:7] = True
    am = rm.copy()
    am[3:6, 2:7, 7] = True
    bm = rm.copy()
    bm[3:6, 2:7, 1] = True
    us = np.array([[0, 0, 1], [0, 1, 0], [0, 0, -1], [0, -1, 0]], float)
    oo = np.tile(np.array([4, 4, 4]) * SP, (4, 1))
    rr = ray_exit(rm, oo, us)
    dr1 = np.abs(ray_exit(am, oo, us) - rr)
    dr2 = np.abs(ray_exit(bm, oo, us) - rr)
    ds1 = np.array([am.sum(), 2 * (am & rm).sum() / (am.sum() + rm.sum()), np.median(dr1), np.quantile(dr1, 0.95), dr1.max()])
    ds2 = np.array([bm.sum(), 2 * (bm & rm).sum() / (bm.sum() + rm.sum()), np.median(dr2), np.quantile(dr2, 0.95), dr2.max()])
    rad = dict(summary=['voxel_count', 'Dice', 'median_abs_delta', 'P95_abs_delta', 'max_abs_delta'], state_a=ds1.tolist(), state_b=ds2.tolist(), identity_error=float(np.max(np.abs(ds1 - ds2))), entire_absolute_distribution_equal=bool(np.array_equal(np.sort(dr1), np.sort(dr2))), downstream_a_mm=float(3.0 - (np.argwhere(am)[:, 2].max() + 0.5) * SP), downstream_b_mm=float(3.0 - (np.argwhere(bm)[:, 2].max() + 0.5) * SP), smallest_extension='Direction or query-specific signed loss', resolution='PER_POINT', time_scale='SIMULTANEOUS', external_referent=dict(kind='our_own_fixture', locator='code/controls.py mirrored cuboid expansions', compared_quantity='radial summary versus point clearance', refutes_us=True))
    rad['downstream_difference_mm'] = rad['downstream_b_mm'] - rad['downstream_a_mm']
    assert rad['identity_error'] == 0 and rad['entire_absolute_distribution_equal']
    scores = [float(i) for i in range(40)]
    (q95, rank) = quantile95(scores)
    scalar = next((x for (j, x) in enumerate(sorted(scores), 1) if j >= math.ceil(41 * 0.95)))
    qe = abs(q95 - scalar)
    (infq, _) = quantile95([0.0] * 10 + [float('inf')] * 30)
    mutations['quantile_0p1mm_poison_rejected'] = abs(q95 + 0.1 - scalar) > 1e-10
    mutations['missing_support_finite_budget_rejected'] = math.isinf(infq)
    mutations['small_sample_finite_budget_rejected'] = math.isinf(quantile95([0.0] * 4)[0])
    profile = json.loads((X8 / 'GUIDE_PROFILES.json').read_text())['profiles'][2]
    g = guide_budget(profile)
    (mu, sd) = profile['apex_mean_sd_mm']
    t = g['apex_bound_mm']
    a = 0.05 / 3
    lowval = mu - sd / math.sqrt((1 - a) / a)
    highval = t
    probs = np.array([1 - a, a])
    vals = np.array([lowval, highval])
    me = float(abs(probs @ vals - mu))
    ve = float(abs(probs @ (vals - mu) ** 2 - sd ** 2))
    mutations['cantelli_optimistic_bound_rejected'] = highval > t - 0.1
    mutations['duplicate_reader_false_independence_rejected'] = not independent_gate(dict(exact_image_equal=True, blinded=True, annotator_ids=['a', 'a'], preconsensus_masks_retained=True))
    mutations['revision_flag_false_independence_rejected'] = not independent_gate(dict(exact_image_equal=True, independent=True))
    mutations['image_one_voxel_poison_rejected'] = not np.array_equal(np.array([1, 2, 3]), np.array([2, 2, 3]))
    mutations['scale_tenfold_poison_rejected'] = not np.array_equal(np.array([0.3] * 3), np.array([3.0, 0.3, 0.3]))
    split = json.loads((ROOT / 'raw/PATIENT_SPLIT.json').read_text())
    assert not set(split['calibration_patients']) & set(split['test_patients'])
    calsh = {byid[p]['image_identity']['tf1']['sha256'] for p in split['calibration_dataset_ids']}
    testsh = {byid[p]['image_identity']['tf1']['sha256'] for p in split['test_dataset_ids']}
    assert not calsh & testsh
    dup = next((r for r in sites if r['alias_source_bindings']))
    bad = copy.deepcopy(dup)
    bad['alias_source_bindings'][0]['image_sha256'] = '0' * 64
    mutations['alias_image_hash_poison_rejected'] = rejects(lambda : validate(bad, byid[bad['patient']], cs[bad['case'], bad['fdi']], certsha, byid=byid))
    bad = copy.deepcopy(dup)
    bad['alias_source_bindings'] = []
    mutations['cross_alias_omission_rejected'] = rejects(lambda : validate(bad, byid[bad['patient']], cs[bad['case'], bad['fdi']], certsha, byid=byid))
    mutations['patient_leakage_poison_rejected'] = bool(set(split['calibration_patients']) & (set(split['test_patients']) | {split['calibration_patients'][0]}))
    controls = dict(claim_type='capability', source_and_cache_semantics_rows=len(sites), exact_image_split_overlap=0, distance=dict(patient=chosen['patient'], kd_edt_max_error_mm=err, kd_brute_max_error_mm=be, tolerance_mm=1e-10, pass_gate=err <= 1e-10 and be <= 1e-10), ray_exit=dict(cuboid_analytic_max_error_mm=re, pass_gate=re <= 1e-10), cylinder=dict(closed_form_mm=2.85, computed_lower_mm=cr['lower_mm'], computed_upper_mm=cr['upper_mm'], max_error_mm=ce, pass_gate=ce <= 1e-10, external_referent=dict(kind='closed_form', locator='axis z cylinder 0<=z<=5, r=1; 0.3-mm cube at(z,y,x)=(2,0,4)', compared_quantity='minimum set distance', refutes_us=False)), quantile=dict(rank=rank, independent_sort_error_mm=qe, pass_gate=qe == 0), cantelli=dict(mean_error=me, variance_error=ve, tail_mass=a, pass_gate=me <= 1e-12 and ve <= 1e-12, external_referent=dict(kind='closed_form', locator='two-point Cantelli extremal distribution given mu, sigma', compared_quantity='first two moments and bound attainability', refutes_us=False)), sufficiency=suff, radial_sufficiency=rad, injected_errors={k: bool(v) for (k, v) in mutations.items()}, formal_floating_point_enclosure='MISSING; analytic identities valid under their assumptions, numerical tests only', equal_information_control='Matches digital distances/rank; no algorithm superiority claim', cost=dict(wall_seconds=time.perf_counter() - start, maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    dump(ROOT / 'raw/CONTROLS.json', controls)
    assert all((controls[k]['pass_gate'] for k in ['distance', 'ray_exit', 'cylinder', 'quantile', 'cantelli'])) and all(controls['injected_errors'].values())
    print(json.dumps(controls, indent=2))
    state('CONTROLS_COMPLETE', 'PASS_WITH_PHYSICAL_INDEPENDENCE_UNKNOWN', 'Build reader demo, figure and final graph feedback')
if __name__ == '__main__':
    main()
