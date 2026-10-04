from dental_release.paths import expand as _release_expand
from common import *
import common, math, csv
from scipy.stats import t as student_t, weibull_min
from scipy.sparse.linalg import spsolve
from scipy.sparse.csgraph import connected_components
from geometry import MET, preparation, field_values
import trimesh
common.R = D / _release_expand('X1B')
X1B = load_module('x1b_frozen_predictor', D / 'LANE_X1B_CROWN_LOOP/code/calibrate.py')
BTE = load_module('bte1_surface_readonly', D / 'LANE_NEXT_E_CEMENT_SQUEEZE/squeeze.py')
X13 = load_module('x13_strip_readonly', D / 'LANE_X13_CEMENT_GAP/cement_port.py')
CAL = read(D / 'LANE_X1B_CROWN_LOOP/raw/CALIBRATION_R1.json')
MODEL = {k: np.asarray(v) if k in ['beta', 'cov', 'basis'] else v for (k, v) in CAL['models']['mixed'].items()}
GROUPS = sorted([r for r in read(D / 'LANE_X1B_CROWN_LOOP/raw/MATCHED_R2.json')['groups'] if r['study'] == 'PMC8558575'], key=lambda r: r['thickness_mm'])

def q05(F, m):
    return float(F * (-np.log(0.95)) ** (1 / m))

def source_quantile(thickness):
    if not np.isfinite(thickness) or thickness < GROUPS[0]['thickness_mm'] or thickness > GROUPS[-1]['thickness_mm']:
        return dict(status='UNKNOWN_OUTSIDE_SOURCE_THICKNESS', force05_N=None, parameter_rectangle_N=None)
    j = min(max(np.searchsorted([r['thickness_mm'] for r in GROUPS], thickness, side='right') - 1, 0), len(GROUPS) - 2)
    (a, b) = GROUPS[j:j + 2]
    u = np.log(thickness / a['thickness_mm']) / np.log(b['thickness_mm'] / a['thickness_mm'])

    def interp(x, y):
        return np.exp((1 - u) * np.log(x) + u * np.log(y))
    f = interp(a['characteristic_load_N'], b['characteristic_load_N'])
    m = interp(a['Weibull_m'], b['Weibull_m'])
    band = [q05(interp(a['F0_90CI_N'][k], b['F0_90CI_N'][k]), interp(a['m_90CI'][k], b['m_90CI'][k])) for k in [0, 1]]
    return dict(status='CONDITIONAL_SOURCE_CURVE', force05_N=q05(f, m), parameter_rectangle_N=band, interpolated_F0_N=float(f), interpolated_m=float(m), anchor_rows=[a['row_id'], b['row_id']], not_joint_confidence_interval=True, physical_transfer_enclosure='MISSING', numerical_enclosure='Mathematical monotone extrema of positive parameter rectangle; directed-rounding enclosure MISSING')

def force_query(row, rec):
    f = load_npz(row['field_path'])
    (p, _) = preparation(rec)
    use = (f['normals'][:, 2] >= 0.5) & (f['sites_mm'][:, 2] >= float(p['margin_z']) + 0.5) & np.isfinite(f['thickness_mm']) & (f['thickness_mm'] > 0)
    if use.sum() < 8:
        return dict(status='UNKNOWN_OCCLUSAL_RAY_SUPPORT', sample_count=int(use.sum()), calibrated_force05_N=None, calibrated_mean_force_N=None)
    ts = f['thickness_mm'][use]
    (lo, t, hi) = np.quantile(ts, [0.05, 0.5, 0.95])
    query = dict(material='3Y', thickness_mm=float(t), angle_deg=0, cement='resin')
    pred = X1B.predict(MODEL, query, n=6)
    x = X1B.feature(query)
    variance = float(x @ MODEL['cov'] @ x + MODEL['tau2'] + MODEL['sigma2'] + MODEL['specimen_log_variance'] / 6)
    direct = float(np.exp(x @ MODEL['beta']))
    half = float(student_t.ppf(0.95, len(MODEL['studies']) - 1) * np.sqrt(variance))
    direct_interval = np.exp(x @ MODEL['beta'] + np.array([-half, half])).tolist()
    assert abs(direct / pred['point_N'] - 1) < 1e-12 and np.allclose(direct_interval, pred['interval_N'], rtol=1e-12)
    supported = [r for r in read(D / 'LANE_X1B_CROWN_LOOP/raw/LITERATURE_GROUPS.json') if r['eligible_R1'] and r['material'] == '3Y' and (r['angle_deg'] == 0) and (r['cement'] == 'resin')]
    span = [min((r['thickness_mm'] for r in supported)), max((r['thickness_mm'] for r in supported))]
    qs = [source_quantile(float(z)) for z in [lo, t, hi]]
    return dict(status='CONDITIONAL_GEOMETRY_QUERY', resolution='PER_TOOTH', sample_count=int(use.sum()), occlusal_thickness_p05_median_p95_mm=[float(lo), float(t), float(hi)], source_force05=qs[1], spatial_scenario_force05=[dict(thickness_mm=float(z), **q) for (z, q) in zip([lo, t, hi], qs)], X1B_mean_force_diagnostic=pred, X1B_thickness_support_mm=span, X1B_thickness_in_support=bool(span[0] <= t <= span[1]), direct_control_relative_error=abs(direct / pred['point_N'] - 1), calibrated_force05_N=None, calibrated_mean_force_N=None, physical_status='UNKNOWN_UNMATCHED_SPECIMEN_FAILED_TRANSFER', model_material_scenario='3Y, resin, axial, six-specimen group; not known manufactured participant material', spatial_summary_sufficiency='FAIL for general mechanics; median alone cannot identify crown strength', missing=['matched material/process/batch', 'contact and support', 'same-geometry fracture measurements', 'local flaw/stress field'], enclosure='No rigorous physical-transfer or full continuous geometry enclosure')

def solve_conductance(v, faces, h, mu=100.0):
    if not np.isfinite(h).all() or h.min() <= 0:
        raise ValueError('NONPOSITIVE_OR_MISSING_FILM')
    used = np.unique(faces)
    mapping = np.full(len(v), -1, int)
    mapping[used] = np.arange(len(used))
    f = mapping[faces]
    vv = v[used] * 0.001
    s = BTE.Surface(vv, f)
    K = s.stiffness(h ** 3 / (12 * mu))
    graph = K.copy()
    graph.data = np.ones_like(graph.data)
    (n, lab) = connected_components(graph, directed=False)
    low = s.boundary
    high = np.setdiff1d(np.flatnonzero(vv[:, 2] >= np.quantile(vv[:, 2], 0.75)), low)
    if not len(low) or not len(high):
        raise ValueError('MISSING_PRESSURE_BOUNDARY')
    for j in range(n):
        if not np.any(lab[low] == j):
            raise ValueError('DISCONNECTED_UNVENTED_FILM')
    fixed = np.union1d(low, high)
    free = np.setdiff1d(np.arange(len(vv)), fixed)
    pressure = np.zeros(len(vv))
    pressure[high] = 1
    pressure[free] = spsolve(K[free][:, free], -K[free][:, high] @ np.ones(len(high)))
    reaction = K @ pressure
    flow = float(reaction[high].sum())
    out = -float(reaction[low].sum())
    mass = abs(flow - out) / max(abs(flow), 1e-300)
    if not np.isfinite(flow) or flow <= 0:
        raise ValueError('NONPOSITIVE_OR_SINGULAR_FLOW')
    return (dict(conductance_m3_Pa_s=flow, mass_relative_error=mass, faces=len(f), vertices=len(vv), low_nodes=len(low), high_nodes=len(high), components=n, area_mm2=float(s.area.sum() * 1000000.0)), s)

def hydraulic_query(row, rec):
    m = load_npz(row['mesh_path'])
    v = m['vertices']
    f = m['faces'][m['face_roles'] == 1]
    (p, pm) = preparation(rec)
    sites = v[f].mean(1)
    (g, cp, n) = field_values(pm, sites)
    path = DATA / 'R7' / rec['uid'] / 'cement_field.npz'
    np.savez_compressed(path, vertices_mm=v, faces=f, sites_mm=sites, normals=n, gap_mm=g)
    base = dict(field_path=str(path), field_sha256=sha(path), resolution='PER_POINT', observed_seated_film_um=None, physical_status='UNKNOWN_UNMEASURED_SEATED_CEMENT', model_state='Nominal dry digital clearance at zero seating displacement', mu_Pa_s=100.0, boundary='0Pa open edge,1Pa top quartile nodes; constructed steady query, not seating', continuous_geometry_enclosure='MISSING', experimental_viscosity_and_seating_state='UNKNOWN')
    try:
        (res, s) = solve_conductance(v, f, g * 0.001)
        mean = float(np.average(g, weights=s.area))
        (control, _) = solve_conductance(v, f, np.full(len(g), mean * 0.001))
        return dict(base, status='CONDITIONAL_REYNOLDS_CONDUCTANCE', **res, uniform_mean_control_m3_Pa_s=control['conductance_m3_Pa_s'], spatial_over_mean_ratio=res['conductance_m3_Pa_s'] / control['conductance_m3_Pa_s'], area_weighted_gap_mean_um=mean * 1000, gap_min_um=float(g.min() * 1000), numerical_pass=res['mass_relative_error'] <= 1e-08)
    except ValueError as e:
        return dict(base, status='UNKNOWN_HYDRAULIC_INPUT', reason=str(e))

def observation_gate(r, observed):
    if not np.isfinite(observed):
        return False
    return abs(observed - r['predicted']) <= 20 if r['kind'] == 'cement' else observed > 0 and abs(np.log(observed / r['predicted'])) <= 0.2

def source_tests():
    frozen = read(V4 / 'FROZEN_LITERATURE_PREDICTIONS.json')
    curated = read(V4 / 'payload/literature/CROWN_CURATED.json')
    cement = {r['row_id']: r for r in csv.DictReader((V4 / 'payload/literature/CEMENT.csv').open())}
    rows = []
    faults = []
    for r in frozen['rows']:
        if r['kind'] == 'cement':
            aa = [cement[i] for i in r['calibration_rows']]
            X = np.array([[1, 1 / float(a['internal_spacer_um'])] for a in aa])
            y = np.array([float(a['measured_mean_um']) for a in aa])
            value = float(np.array([1, 1 / r['x']]) @ np.linalg.solve(X, y))
            bad = r['predicted'] + 1000
        else:
            aa = sorted([a for a in curated if a['study'] == r['study'] and a['protocol'] == r['protocol'] and (a['product'] == r['product']) and (a['t'] in r['calibration_thicknesses_mm'])], key=lambda a: a['t'])
            X = np.array([[1, np.log(a['t'])] for a in aa])
            y = np.log([a['mean'] for a in aa])
            value = float(np.exp(np.array([1, np.log(r['x'])]) @ np.linalg.solve(X, y)))
            bad = r['predicted'] * 2
        assert abs(value / r['predicted'] - 1) < 1e-10
        rows.append(dict(r, predicted_recomputed=value, gate=observation_gate(r, r['observed']), quantity='dry_marginal_gap_group_mean' if r['kind'] == 'cement' else 'crown_group_force_mean', resolution='POPULATION'))
        faults.append(dict(normal_pass=observation_gate(r, r['predicted']), name=r.get('row_id', r['study'] + str(r['x'])), control='actual observed-value gate', fault='observed+1000um or force=2*prediction', rejected=not observation_gate(r, bad)))
    qchecks = []
    for r in GROUPS:
        q = q05(r['characteristic_load_N'], r['Weibull_m'])
        cdf = float(weibull_min.cdf(q, r['Weibull_m'], scale=r['characteristic_load_N']))
        wrong = float(weibull_min.ppf(0.95, r['Weibull_m'], scale=r['characteristic_load_N']))
        qchecks.append(dict(row=r['row_id'], force05_N=q, cdf=cdf, correct=abs(cdf - 0.05) < 1e-12, wrong95_rejected=abs(weibull_min.cdf(wrong, r['Weibull_m'], scale=r['characteristic_load_N']) - 0.05) > 1e-12, locator=r['locator']))
    return dict(rows=rows, source_faults=faults, quantile_controls=qchecks, X1B_original_validation_gates=CAL['gates'], X1B_original_scores=CAL['scores'], external_referent=read(P / 'PREREG_R1.json')['external_referent'])

def sufficiency():
    a = np.array([0.5, 1.5])
    b = a[::-1]
    demand = np.array([1.0, 2.0])
    m = 7.2
    H = lambda t: float(np.mean((demand / t ** 2) ** m))
    fa = (-np.log(0.95) / H(a)) ** (1 / m)
    fb = (-np.log(0.95) / H(b)) ** (1 / m)
    (low, high) = (3.2e-05, 9.6e-05)
    runs = [X13.strip(BTE, 16, low, high, k) for k in ['series', 'parallel']]
    film1 = np.tile(np.array([32.0, 96.0]), 128)
    film2 = np.sort(film1)
    summary1 = float(np.mean(film1))
    summary2 = float(np.mean(film2))
    ratio = runs[1]['conductance_m3_Pa_s'] / runs[0]['conductance_m3_Pa_s']
    return dict(thickness=dict(summary='mean,median and sorted thickness multiset', summary_a=float(a.mean()), summary_b=float(b.mean()), identity_error=float(abs(a.mean() - b.mean())), downstream='normalized weakest-link force05 under same asymmetric demand', downstream_a=float(fa), downstream_b=float(fb), difference=float(abs(fb - fa)), ratio=float(fb / fa), minimum_extension='For this fixed demand/exponent, one load-weighted hazard scalar H=sum(w*(demand/t^2)^m). For unknown demand, retain joint position/thickness/stress field.', claim='Conditional mathematical counterexample, not measured crown'), film=dict(summary='unweighted equal-cell mean and thickness multiset in um', summary_a=summary1, summary_b=summary2, identity_error=abs(summary1 - summary2), downstream='steady conductance m3/(Pa s)', difference=runs[1]['conductance_m3_Pa_s'] - runs[0]['conductance_m3_Pa_s'], ratio=ratio, runs=runs, minimum_extension='Fixed-flow effective resistance (one scalar); to change boundaries retain film connectivity and positions', external_referent=dict(kind='closed_form', locator='doi:10.1098/rstl.1886.0005; X13 DECOMPOSITION_R4.json', compared_quantity='series/parallel Reynolds conductance', refutes_us=True)), controls=[dict(name=k['kind'], relative_error=k['relative_error'], pass_gate=k['relative_error'] <= 0.02, mass_gate=k['mass_relative_error'] <= 1e-08, injected_2x_rejected=abs(2 * k['conductance_m3_Pa_s'] / k['exact_m3_Pa_s'] - 1) > 0.02) for k in runs])

def run(rows):
    import function_axis
    records = {r['uid']: r for r in cohort()}
    out = []
    start = time.perf_counter()
    for (i, r) in enumerate(rows):
        if r['status'] != 'ADAPTED':
            out.append(dict(r, physics_status='UNKNOWN_NO_GEOMETRY'))
            continue
        rec = records[r['uid']]
        st = time.perf_counter()
        out.append(dict(r, force=force_query(r, rec), cement=hydraulic_query(r, rec), adapted_function=function_axis.query(r, rec), physics_seconds=time.perf_counter() - st))
        dump(P / 'raw/PHYSICS_CHECKPOINT.json', out)
        if i % 20 == 0:
            print('physics', i + 1, flush=True)
    predictions = dict(implementation_hashes={name: sha(P / 'code' / name) for name in ['common.py', 'geometry.py', 'remesh_r7.py', 'exact_margin_r7.py', 'physics.py', 'function_axis.py']}, claim_type='capability', source_hashes=dict(X1B=sha(D / 'LANE_X1B_CROWN_LOOP/raw/CALIBRATION_R1.json'), Weibull=sha(D / 'LANE_X1B_CROWN_LOOP/raw/MATCHED_R2.json')), prereg_sha256=sha(P / 'PREREG_R7.json'), rows=[{k: v for (k, v) in r.items() if k not in ['seconds', 'physics_seconds']} for r in out], measurement_status='NO_NEW_PHYSICAL_MEASUREMENT; geometry already observed retrospectively; source outcomes already read', physical_quantities_calibrated=False)
    freeze(P / 'FROZEN_PREDICTIONS_R7.json', predictions)
    dump(P / 'FROZEN_PREDICTIONS.json', dict(schema='immutable_prediction_versions', versions=[dict(path='FROZEN_PREDICTIONS_R4.json', sha256=sha(P / 'FROZEN_PREDICTIONS_R4.json'), status='R4_NUMERICAL_ONLY_EXPORT_FAILED'), dict(path='FROZEN_PREDICTIONS_R7.json', sha256=sha(P / 'FROZEN_PREDICTIONS_R7.json'), status='FINAL_R7_CONDITIONAL')], measurement_status='NO_MATCHED_PHYSICAL_MEASUREMENTS'))
    (P / 'FROZEN_PREDICTIONS.sha256').write_text(sha(P / 'FROZEN_PREDICTIONS.json') + '\n')
    dump(P / 'raw/PHYSICAL_ROWS.json', out)
    dump(P / 'raw/SOURCE_VALIDATION.json', source_tests())
    dump(P / 'raw/SUFFICIENCY.json', sufficiency())
    return out
