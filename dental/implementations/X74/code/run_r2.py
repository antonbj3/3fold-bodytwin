import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '4'
import copy, itertools, json, math, resource, time
import numpy as np
import mpmath as mp
import trimesh
from scipy.stats import t
from common import ROOT, dump, sha, stamp, check_frozen, inputs, state
from sources import read_sources, source_error
from optics import de00, interval_de00
from run_r1 import tensor_data

def classify(box, threshold, kind):
    (lo, hi) = box
    if kind == 'at_least':
        return 'PASS' if lo >= threshold else 'FAIL' if hi < threshold else 'UNCERTAIN'
    return 'PASS' if hi <= threshold else 'FAIL' if lo > threshold else 'UNCERTAIN'

def ratio_box(box, q):
    mp.iv.dps = 35
    v = mp.iv.mpf([str(box[0]), str(box[1])]) / mp.iv.mpf(str(q))
    return [float(np.nextafter(float(v.a), -np.inf)), float(np.nextafter(float(v.b), np.inf))]

def brute_ray(triangles, origin, direction):
    e1 = triangles[:, 1] - triangles[:, 0]
    e2 = triangles[:, 2] - triangles[:, 0]
    h = np.cross(direction, e2)
    det = np.einsum('ij,ij->i', e1, h)
    ok = abs(det) > 1e-12
    inv = np.zeros_like(det)
    inv[ok] = 1 / det[ok]
    s = origin - triangles[:, 0]
    u = inv * np.einsum('ij,ij->i', s, h)
    r = np.cross(s, e1)
    v = inv * (r @ direction)
    dist = inv * np.einsum('ij,ij->i', e2, r)
    ok &= (u >= -1e-10) & (v >= -1e-10) & (u + v <= 1 + 1e-10) & (dist > 1e-09)
    return float(dist[ok].min()) if ok.any() else None

def ray_checks(rays):
    mesh = trimesh.load_mesh(ROOT / 'raw/crown.stl', process=True)
    h = rays['thickness_mm']
    checks = []
    for k in np.linspace(0, len(h) - 1, 12).astype(int):
        origin = rays['points'][k] - 0.01 * rays['normals'][k]
        direction = -rays['normals'][k]
        control = brute_ray(mesh.triangles, origin, direction)
        control = None if control is None else control + 0.01
        err = abs(float(h[k]) - control) if control is not None and np.isfinite(h[k]) else None
        baderr = abs(float(h[k]) + 0.02 - control) if err is not None else None
        checks.append({'ray_id': int(k), 'h_mm': float(h[k]) if np.isfinite(h[k]) else None, 'independent_h_mm': control, 'error_mm': err, 'actual_injected_plus0_02_error_mm': baderr})
    assert any((c['error_mm'] is not None for c in checks))
    return checks

def query_optical(row, h_mm, product, measured_patient=False):
    if h_mm != 1.0:
        return {'status': 'UNKNOWN', 'reason': 'unsupported thickness: source only samples 1 mm'}
    if product != 'Aidite 3D Pro Zir A1':
        return {'status': 'UNKNOWN', 'reason': 'different zirconia product'}
    if not measured_patient:
        return {'status': 'UNKNOWN', 'reason': 'coupon mean is not measured local patient Lab; cement-film thickness unavailable'}
    return {'status': 'LOCAL_MEASUREMENT_REQUIRED', 'reason': 'interface supplied; model still requires matched measurement and review'}

def main():
    started = time.perf_counter()
    p = check_frozen('PREREG_R2.json')
    inputs()
    (bond, optical, target) = read_sources()
    suff = json.loads((ROOT / 'raw/SUFFICIENCY_R1.json').read_text())
    assert suff['bond']['identity_error'] == suff['optical']['identity_error'] == 0
    assert source_error(bond, optical, target) <= 1e-10
    bond_out = []
    for r in bond:
        se = r['SD_MPa'] / math.sqrt(r['n'])
        dt = float(t.ppf(0.975, r['n'] - 1)) * se
        interval = [r['mean_MPa'] - dt, r['mean_MPa'] + dt]
        bond_out.append({**r, 'pointwise_95pct_mean_CI_MPa': interval, 'lab_challenge_MPa': 5.0, 'lab_mean_screen': classify(interval, 5.0, 'at_least'), 'patient_screen': 'UNKNOWN: nominal coupon SBS is not a local mixed-mode cohesive law'})
    optical_out = []
    interval_checks = []
    for (i, r) in enumerate(optical):
        box = [[mu - 2 * sd, mu + 2 * sd] for (mu, sd) in zip(r['mean_Lab'], r['SD_Lab'])]
        bounds = interval_de00(box, target)
        centre = de00(r['mean_Lab'], target)
        probes = []
        for w in itertools.product([0.0, 0.5, 1.0], repeat=3):
            a = [lo + v * (hi - lo) for ((lo, hi), v) in zip(box, w)]
            probes.append(de00(a, target))
        contains = all((bounds[0] <= v <= bounds[1] for v in probes))
        bad_bounds = [centre, centre]
        injected_rejected = not all((bad_bounds[0] <= v <= bad_bounds[1] for v in probes))
        interval_checks.append({'row': i, 'contains_all_27': contains, 'actual_zero_width_injection_rejected': injected_rejected, 'probe_range': [min(probes), max(probes)], 'enclosure': bounds})
        optical_out.append({**r, 'de00_of_mean_Lab': centre, 'Lab_box_mean_plusminus2SD': box, 'rigorous_formula_enclosure_de00': bounds, 'source_specific_screen': classify(bounds, 3.6, 'at_most'), 'source_threshold': 3.6, 'note': 'de00(mean Lab) is not mean de00 of raw specimens; rectangular box is descriptive, no probability assigned'})
    assert all((c['contains_all_27'] and c['actual_zero_width_injection_rejected'] for c in interval_checks))
    (d, n, q, traction) = tensor_data()
    tied = d['tied']
    demands = {k: float(v[tied].max()) for (k, v) in q.items()}
    inverse_targets = []
    conditional = []
    for (lc, peak) in demands.items():
        arg = np.flatnonzero(tied)[np.argmax(q[lc][tied])]
        inverse_targets.append({'loadcase': lc, 'facet_id': int(arg), 'coordinate_mm': d['cf'][arg].astype(float).tolist(), 'normal': n[arg].tolist(), 'conditional_tau_needed_at_150N_MPa': 150 * peak, 'force_challenge_N': 150.0, 'resolution_level': 'PER_POINT', 'time_scale': 'SIMULTANEOUS', 'status': 'MEASUREMENT_TARGET_ONLY; tied FE traction not yet a physical prediction'})
        for r in bond_out:
            if r['with_cover']:
                cap = ratio_box(r['pointwise_95pct_mean_CI_MPa'], peak)
                conditional.append({'system': r['system'], 'loadcase': lc, 'conditional_quotient_N': cap, 'q_MPa_per_N': peak, 'status': 'PHENOMENOLOGICAL_SCREEN_ONLY; no physical FE/nominal-SBS transfer enclosure'})
    rays = np.load(ROOT / 'raw/rays.npz')
    checks = ray_checks(rays)
    maxray = max((c['error_mm'] for c in checks if c['error_mm'] is not None))
    actual_bad = all((c['actual_injected_plus0_02_error_mm'] > 1e-08 for c in checks if c['error_mm'] is not None))
    assert maxray <= 1e-08 and actual_bad, (maxray, checks)
    h = rays['thickness_mm']
    finite = np.isfinite(h)
    visible = rays['visible']
    supported = finite & (h == 1.0)
    coverage = {'all_rays': len(h), 'finite_rays': int(finite.sum()), 'visible_rays': int(visible.sum()), 'exact_1mm_rays': int(supported.sum()), 'exact_1mm_visible_rays': int((supported & visible).sum()), 'patient_rejected': len(h), 'patient_rejection_fraction': 1.0, 'reason': 'all are missing matched optical product/cement-film/patient colour; exact single-thickness support is also absent', 'matched_bond_test_product': False, 'matched_optical_test_product': False, 'patient_result': 'UNKNOWN'}
    missing_faults = [query_optical(optical[0], 1.01, 'Aidite 3D Pro Zir A1'), query_optical(optical[0], 1.0, 'generic 5Y'), query_optical(optical[0], 1.0, 'Aidite 3D Pro Zir A1')]
    assert all((x['status'] == 'UNKNOWN' for x in missing_faults))
    product_join = {'bond_zirconia': '3M Lava Esthetic', 'optical_zirconia': 'Aidite 3D Pro Zir A1', 'bond_systems': ['Scotchbond Universal + RelyX Ultimate', 'Prime & Bond Universal + SmartCem2', 'Optibond Universal + Maxcem Elite Chroma'], 'optical_cement': 'GC G-CEM LinkForce', 'joinable': False, 'simultaneously_qualified_product_choices': 0, 'status': 'UNKNOWN: new paired measurement required; shared patient geometry does not imply shared material population'}
    payload = {'patient_case': 'STS ROI_L005 lo_M1_k6 D1', 'conditional_force_quotients_N': conditional, 'measurement_targets': inverse_targets, 'prospective_measurement_status': 'NOT_RUN', 'forecast_type': 'Conditional screening instrument, not calibrated physical forecast', 'missing_measurements': ['local mixed-mode cohesive law and stiffness with curing/aging on matching material', 'local stack Lab or reflectance spectra and cement-film thickness on same fabricated crown']}
    from run_r1 import freeze_predictions
    frozen = ROOT / 'FROZEN_PREDICTIONS.json'
    if frozen.exists():
        assert check_frozen('FROZEN_PREDICTIONS.json')['predictions'] == payload
    else:
        dump(frozen, {'frozen_utc': stamp(), 'prereg_sha256': sha(ROOT / 'PREREG_R2.json'), 'code_sha256': sha(__file__), 'predictions': payload})
        (ROOT / 'FROZEN_PREDICTIONS.json.sha256').write_text(sha(frozen) + '\n')
    with (ROOT / 'raw/LOCAL_BOND_DEMAND.csv').open('w') as f:
        f.write('facet,x_mm,y_mm,z_mm,area_mm2,tied,q_axial_MPa_per_N,q_offaxis30_MPa_per_N,resolution_level,time_scale\n')
        for i in range(len(n)):
            f.write(','.join(map(str, [i, *d['cf'][i].tolist(), float(d['af'][i]), int(tied[i]), float(q['axial'][i]), float(q['offaxis30'][i]), 'PER_POINT', 'SIMULTANEOUS'])) + '\n')
    forecasts = [{'substrate': s, 'cemented_nominal_pass_count': sum((r['de00_of_mean_Lab'] <= 3.6 for r in optical_out if r['substrate'] == s and r['cement'] != 'CG')), 'cemented_box_screens': {r['cement']: r['source_specific_screen'] for r in optical_out if r['substrate'] == s and r['cement'] != 'CG'}} for s in ['ZR', 'CR', 'MT']]
    out = {'round': 'R2', 'claim_type': 'information_link', 'bond_lab_screens': bond_out, 'optical_lab_screens': optical_out, 'rescue_queries': forecasts, 'conditional_force_quotients': conditional, 'inverse_measurement_targets': inverse_targets, 'patient_coverage': coverage, 'product_join': product_join, 'checks': {'intervals': interval_checks, 'rays': checks, 'max_ray_error_mm': maxray, 'actual_ray_plus0_02mm_rejected': actual_bad, 'missing_state_injections': missing_faults, 'unsupported_product_join_rejected': not product_join['joinable'], 'traction': traction, 'R1_same_summary_identity_still_zero': True}, 'external_referent': {'kind': 'independent_measurement', 'locator': ['doi:10.1016/j.jds.2023.05.011#Table5', 'doi:10.1016/j.heliyon.2023.e23046#Table2'], 'compared_quantity': 'Source-specific SBS and Lab, POPULATION; no matched patient output', 'refutes_us': True}, 'claim_limit': 'The stronger spatial/source contract is runnable. Empirical source decisions are valid in source regimes. Patient bond capacity/colour and joint product decision remain UNKNOWN.', 'cost': {'elapsed_s': time.perf_counter() - started, 'max_rss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'threads_max': 4}}
    dump(ROOT / 'rounds/R2.json', out)
    (ROOT / 'HANDOFF_R2.md').write_text("R2 completed: measured system/layer is maintained in local ports; rigorous colour inclusion applies to the declared Lab box. All local patient issues refrain from missing product/regime/measuring support.\nNext design: export the same patient's local keyframes as customizable optical and shear specimen with fixed product identity, and a scooter that requires these measurementss before UNKNOWN can be changed. The R1 faults and crude facits are preserved.\n")
    state(status='R2_DECIDED', latest_gate={'source_ports': 'PASS', 'patient_transfer': 'UNKNOWN', 'joint_product_join': 'UNKNOWN'}, next_operation='R3: runnable matched laboratory specimen export and measurement ingestion contract, with actual adverse input rejection')
    print(json.dumps({'covered_bond': [(r['system'], r['lab_mean_screen'], r['pointwise_95pct_mean_CI_MPa']) for r in bond_out if r['with_cover']], 'optical': forecasts, 'coverage': coverage, 'ray_max_error_mm': maxray}))
if __name__ == '__main__':
    main()
