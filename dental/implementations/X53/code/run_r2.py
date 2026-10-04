import time, resource
from common import *
from collision import *
from run_r1 import lib

def closed_cap(v, f):
    m = trimesh.Trimesh(v, f, process=False)
    if not m.is_watertight:
        edges = m.edges_sorted
        (unique, count) = np.unique(edges, axis=0, return_counts=True)
        boundary = unique[count == 1]
        ids = np.unique(boundary)
        centre = np.asarray(v)[ids].mean(0)
        nv = np.r_[v, [centre]]
        nf = np.r_[f, np.c_[boundary, np.full(len(boundary), len(v))]]
        m = trimesh.Trimesh(nv, nf, process=True)
        m.fix_normals()
    if not m.is_watertight or m.volume <= 0:
        raise ValueError('die topology invalid')
    return m

def normal_gaps(die, points, normals):
    (loc, idx, _) = die.ray.intersects_location(points, normals, multiple_hits=True)
    out = np.full(len(points), np.nan)
    if len(loc):
        ts = np.einsum('ij,ij->i', loc - points[idx], normals[idx])
        for i in range(len(points)):
            t = ts[(idx == i) & (ts >= 0)]
            if len(t):
                out[i] = t.min()
    return out

def controls():
    h0 = np.tile(np.array([60.0, 100.0]), 32)
    stock_A = np.tile(np.array([0.0, 40.0]), 32)
    stock_B = stock_A[::-1]
    hA = h0 - stock_A
    hB = h0 - stock_B
    error = abs(hA.mean() - hB.mean())
    threshold = read(ROOT / 'PREREG_R2.json')['metrics']['obstruction_scenario_gap_um']
    fracA = float(np.mean(hA < threshold))
    fracB = float(np.mean(hB < threshold))
    assert error == 0.0 and abs(fracA - fracB) >= 0.1
    p = np.array([[0.0, 0.0, 0.0]])
    n = np.array([[0.0, 0.0, 1.0]])
    observed = first_union_ball_entry(p, n, [[0, 0, 0.55]], [0.5])[0] * 1000
    exact = (0.55 - 0.5) * 1000
    assert abs(observed - exact) < 1e-06 and abs(observed + 50 - exact) > 1e-06
    return dict(summary_identity_error_um=error, stock_mean_identity_error_um=float(abs(stock_A.mean() - stock_B.mean())), mean_gap_A_um=float(hA.mean()), mean_gap_B_um=float(hB.mean()), minimum_gap_A_um=float(hA.min()), minimum_gap_B_um=float(hB.min()), obstruction_fraction_A=fracA, obstruction_fraction_B=fracB, downstream_obstruction_fraction_difference=abs(fracA - fracB), minimum_extension='Joint PER_POINT registration of nominal film and signed stock, including normal directions; a mean or independent stock histogram is insufficient', sphere_ray_absolute_error_um=abs(observed - exact), injected_50um_rejected=True, constructed_fields=True, empirical=False, external_referent=dict(kind='closed_form', locator='https://mathworld.wolfram.com/Sphere.html; quadratic line/sphere substitution explicitly in DECOMPOSITION_R2.json', compared_quantity='sphere normal-ray entry in mm, not manufactured fit', refutes_us=True))

def run():
    pre = read(ROOT / 'PREREG_R2.json')
    assert sha(ROOT / 'PREREG_R2.json') == (ROOT / 'PREREG_R2.sha256').read_text().strip()
    assert sha(ROOT / 'FROZEN_PREDICTIONS.json') == pre['frozen_R1_sha256']
    assert sha(ROOT / 'LAB_CASES_20.json') == pre['lab_cases_sha256']
    ctl = controls()
    dump(ROOT / 'raw/SUFFICIENCY_R2.json', ctl)
    start = time.perf_counter()
    manifest = {r['id']: r for r in read(ROOT / 'INPUT_MANIFEST.json')['rows']}
    rows = []
    hashes = {}
    for case in read(ROOT / 'LAB_CASES_20.json')['cases']:
        dest = DATA / 'R2' / (case['design_id'] + '.json')
        if dest.exists():
            out = read(dest)
            assert out['prereg_sha256'] == sha(ROOT / 'PREREG_R2.json')
            rows.append(out['summary'])
            hashes[str(dest)] = sha(dest)
            continue
        row = manifest[case['design_id']]
        a = np.load(row['file'], allow_pickle=False)
        scene = Scene(a['vertices'], a['faces'])
        r1 = read(case['prediction_file'])
        variant = next((v for v in r1['points_and_tools'] if v['library'] == 'vhf' and v['axes'] == 5))
        rec = [r for r in variant['records'] if r['region'] == 'intaglio']
        P = np.array([r['point'] for r in rec])
        N = np.array([r['normal'] for r in rec])
        tools = lib('vhf', 5)
        selected = []
        centres = []
        radii = []
        st = time.perf_counter()
        for (rec0, p, n) in zip(rec, P, N):
            available = [r for r in rec0['tools'] if r['status'] == 'FOUND']
            if available:
                w = min(available, key=lambda r: r['diameter_design_mm'])
            else:
                w = None
                for off in pre['metrics']['offset_grid_mm']:
                    for t in tools:
                        v = search(scene, p, n, t, 5, offsets=(off,))
                        if v['status'] == 'FOUND':
                            w = dict(**v, tool_id=t['id'], diameter_design_mm=t['diameter_mm'])
                            break
                    if w is not None:
                        break
            if w is not None:
                centres.append(w['centre'])
                radii.append(w['diameter_design_mm'] / 2)
            selected.append(w)
        stock = first_union_ball_entry(P, N, centres, radii)
        forced = first_union_ball_entry(P, N, P + tools[0]['diameter_mm'] / 2 * N, np.full(len(P), tools[0]['diameter_mm'] / 2))
        if row.get('die'):
            if sha(row['die']) != row['die_sha256']:
                raise ValueError('die drift')
            die = trimesh.load_mesh(row['die'], process=True)
        else:
            die = closed_cap(a['prep_vertices'], a['prep_faces'])
        h0 = normal_gaps(die, P, N)
        h_safe = h0 - stock
        h_forced = h0 - forced
        records = []
        for i in range(len(P)):
            records.append(dict(point=P[i], normal=N[i], face_id=rec[i]['face_id'], normal_gap_nominal_um=h0[i] * 1000, stock_safe_ball_union_um=stock[i] * 1000, forced_signed_stock_um=forced[i] * 1000, normal_gap_after_safe_union_um=h_safe[i] * 1000, normal_gap_after_forced_union_um=h_forced[i] * 1000, safe_pose_witness=selected[i], stock_status='FINITE_BALL_WITNESS_SET' if np.isfinite(stock[i]) else 'UNKNOWN_NO_RAY_CUT', fit_status='CONDITIONAL_FIXED_NORMAL_RAY' if np.isfinite(h0[i]) else 'UNKNOWN_NO_DIE_RAY', forced_status='INFEASIBLE_OR_UNVERIFIED_COMPENSATION_DIAGNOSTIC', geometry_resolution='PER_POINT', normal_ray_window_truncated=bool(np.isfinite(forced[i]) and forced[i] <= -0.4 + EPS)))

        def frac(x):
            return float(np.mean(x))
        sums = dict(case_id=case['case_id'], id=row['id'], family=row['family'], n=len(P), stand_off_new_witnesses=sum((rec[i]['smallest_found_green_mm'] is None and selected[i] is not None for i in range(len(P)))), stock_gt25um_fraction=frac(stock > 0.025), safe_union_overcut_gt25um_fraction=frac(stock < -0.025), forced_overcut_gt25um_fraction=frac(forced < -0.025), unknown_stock_fraction=frac(~np.isfinite(stock)), unknown_nominal_gap_fraction=frac(~np.isfinite(h0)), safe_nonpositive_gap_fraction=frac(h_safe <= 0), nominal_gap_median_um=float(np.nanmedian(h0) * 1000), safe_gap_median_um=float(np.nanmedian(h_safe) * 1000), max_safe_standoff_um=max((w['offset_mm'] * 1000 for w in selected if w is not None), default=None), seconds=time.perf_counter() - st, empirical_fit='UNKNOWN', scope='Sparse cutting-ball union, fixed ray, prescribed neck/holder; not full CAM swept stock or seated gap')
        out = dict(id=row['id'], prereg_sha256=sha(ROOT / 'PREREG_R2.json'), source_sha256=row['sha256'], records=records, summary=sums)
        dump(dest, out)
        rows.append(sums)
        hashes[str(dest)] = sha(dest)
        dump(ROOT / 'raw/R2_PARTIAL.json', rows)
        state('R2_FIELD_HANDOVER_RUNNING', case['case_id'] + ' saved', 'Continue local material-change maps; no physical gap claim')
        print(case['case_id'], row['id'], 'new stand-offs', sums['stand_off_new_witnesses'], 'forced overcut fraction', sums['forced_overcut_gt25um_fraction'], flush=True)
    result = dict(round='R2', claim_type='capability', construction_gate='PASS' if ctl['summary_identity_error_um'] == 0 and len(rows) == 20 else 'FAIL', physical_gate='UNKNOWN_UNPAIRED_CAM_AND_GEOMETRY', controls=ctl, rows=rows, raw_files=hashes, seconds_this_invocation=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, external_referent=pre['external_referent'], affine_sensitivity='Not used; rigorously rounded collision/ray arithmetic remains missing')
    dump(ROOT / 'RESULTS_R2.json', result)
    freeze(ROOT / 'FROZEN_PREDICTIONS_R2.json', dict(prereg_sha256=sha(ROOT / 'PREREG_R2.json'), result_sha256=sha(ROOT / 'RESULTS_R2.json'), point_files=hashes, physical_measurement_performed=False, code_sha256=sha(Path(__file__))))
    return result
if __name__ == '__main__':
    run()
