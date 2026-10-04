import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '4'
from pathlib import Path
import datetime, itertools, json, math, time, resource
import numpy as np
from mpmath import iv
from freeze import ROOT, sha, dump
from run_r1 import validate_inputs
iv.dps = 30
OUT = Path(os.environ.get('X64_REPLAY_DIR', str(ROOT)))
(OUT / 'raw').mkdir(parents=True, exist_ok=True)

def iv_eval(h, a, b, y0, y1):
    (h, a, b, y0, y1) = (iv.mpf(float(x)) for x in (h, a, b, y0, y1))
    w = (h - a) / (b - a)
    return iv.exp((1 - w) * iv.ln(y0) + w * iv.ln(y1))

def endpoints(v):
    return (math.nextafter(float(v.a), -math.inf), math.nextafter(float(v.b), math.inf))

def corners(boxes):
    return list(itertools.product(*boxes))

def enclosure(boxes):
    vals = [endpoints(iv_eval(*c)) for c in corners(boxes)]
    return (min((v[0] for v in vals)), max((v[1] for v in vals)))

def direct(h, a, b, y0, y1):
    w = (h - a) / (b - a)
    return np.exp((1 - w) * np.log(y0) + w * np.log(y1))

def main():
    tic = time.perf_counter()
    (p1, manifest) = validate_inputs()
    p = json.loads((ROOT / 'PREREG_R2.json').read_text())
    assert sha(ROOT / 'PREREG_R2.json') == (ROOT / 'PREREG_R2.sha256').read_text().strip()
    assert sha(ROOT / 'raw/RESULTS_R1.json') == p['parent_results_sha256']
    source = json.loads((ROOT / 'raw/SOURCE_ROWS.json').read_text())
    old = json.loads((ROOT / 'raw/RESULTS_R1.json').read_text())
    rays = np.load(ROOT / 'raw/OPTICAL_RAYS.npz')
    h = rays['thickness_mm']
    visible = rays['visible']
    points = rays['points']
    delta = p['TM02']['geometry_scenario_halfwidth_mm']
    a = p['TM02']['knot0_mm_box']
    b = p['TM02']['knot1_mm_box']
    supported = np.isfinite(h) & (h - delta >= a[1]) & (h + delta <= b[0])
    rng = np.random.default_rng(6402)
    optical = []
    all_contained = True
    bad_rejected = True
    for route in ('ST', 'SP'):
        rr = [r for r in source['optical'] if r['route'] == route]
        y0 = [rr[0]['observed_min_TP'], rr[0]['observed_max_TP']]
        y1 = [rr[1]['observed_min_TP'], rr[1]['observed_max_TP']]
        assert y0[0] > y1[1] > 0, 'Required monotonic source boxes'
        data = []
        max_sample_violation = 0.0
        for i in np.flatnonzero(visible):
            row = {'ray_id': int(i), 'point_mm': points[i].tolist(), 'length_mm': float(h[i]) if np.isfinite(h[i]) else None, 'resolution_level': 'PER_POINT', 'time_scale': 'HANDOVER'}
            if supported[i]:
                boxes = [[h[i] - delta, h[i] + delta], a, b, y0, y1]
                (lo, hi) = enclosure(boxes)
                c = np.array(corners(boxes))
                sample = rng.uniform(np.array(boxes)[:, 0], np.array(boxes)[:, 1], size=(200, 5))
                v = np.concatenate([direct(*c.T), direct(*sample.T)])
                violation = max(0.0, lo - v.min(), v.max() - hi)
                max_sample_violation = max(max_sample_violation, float(violation))
                contained = bool((v >= lo - 1e-12).all() and (v <= hi + 1e-12).all())
                reject = bool(v.max() > lo + 0.5 * (hi - lo))
                all_contained &= contained
                bad_rejected &= reject
                state = 'MODEL_BOX_PASS' if lo >= 12 else 'MODEL_BOX_FAIL' if hi < 12 else 'MODEL_BOX_UNRESOLVED'
                row.update({'TP_model_box': [lo, hi], 'required_beta_lower_TP': 12 - lo, 'screen': state, 'interior_and_corner_containment': contained, 'injected_narrow_upper_rejected': reject})
            else:
                row.update({'TP_model_box': None, 'required_beta_lower_TP': None, 'screen': 'OUTSIDE_SOURCE_SUPPORT', 'reason': 'Whole +/-0.12mm scenario ray box not within every measured knot interval'})
            data.append(row)
        counts = {k: sum((r['screen'] == k for r in data)) for k in ('MODEL_BOX_PASS', 'MODEL_BOX_FAIL', 'MODEL_BOX_UNRESOLVED', 'OUTSIDE_SOURCE_SUPPORT')}
        optical.append({'route': route, 'duration_min': rr[0]['duration_min'], 'rows': data, 'counts': counts, 'max_sample_box_violation_TP': max_sample_violation, 'physical_D1_release': 'UNKNOWN', 'rigorous_scope': 'Declared log-affine surrogate/source boxes only; actual optical and geometry errors unbounded'})
    demands = np.load(ROOT / 'raw/DEMANDS.npz')
    mechanical = []
    frozen = []
    for name in ('axial', 'offaxis30'):
        d = demands[name + '_MPa']
        imax = int(np.argmax(d))
        peak = float(d[imax])
        rows = []
        for r in source['cure']:
            sbox = [max(0.0, r['mean_MPa'] - 2 * r['SD_MPa']), r['mean_MPa'] + 2 * r['SD_MPa']]
            (lo, hi) = endpoints(iv.mpf(sbox) / iv.mpf(peak))
            positive = d > 1e-12
            directlo = float((sbox[0] / d[positive]).min())
            directhi = float((sbox[1] / d[positive]).min())
            agree = lo <= directlo <= hi and lo <= directhi <= hi
            rejects = bool(2 * lo > directlo or 0.5 * hi < directhi)
            all_contained &= bool(agree)
            bad_rejected &= rejects
            state = 'MODEL_BOX_PASS' if lo >= 1 else 'MODEL_BOX_FAIL' if hi < 1 else 'MODEL_BOX_UNRESOLVED'
            predicted_N = r['mean_MPa'] / (peak / 150.0)
            forcebox = [sbox[0] / (peak / 150.0), sbox[1] / (peak / 150.0)]
            rows.append({'cure_min': r['cure_min'], 'strength_MPa_box': sbox, 'ratio_box': [lo, hi], 'alpha_required_at_mean': peak / r['mean_MPa'], 'screen': state, 'direct_facet_control_enclosed': bool(agree), 'injected_bad_ratio_box_rejected': rejects, 'resolution_level': 'PER_SURFACE_REGION', 'time_scale': 'HANDOVER'})
            frozen.append({'cure_min': r['cure_min'], 'loadcase': name, 'nominal_demand_MPa_at150N': peak, 'source_strength_mean_MPa': r['mean_MPa'], 'predicted_fracture_proxy_N': predicted_N, 'source_box_fracture_proxy_N': forcebox, 'critical_face_id': imax, 'critical_point_mm': demands['cf'][imax].tolist(), 'contact_zone_retained': True, 'comparison_endpoint': 'F_at_first_fracture on matched D1 geometry/loadcase; proxy validity UNKNOWN', 'physical_prediction_interval': None, 'resolution_level': 'PER_TOOTH', 'time_scale': 'HANDOVER'})
        passed = [r['cure_min'] for r in rows if r['screen'] == 'MODEL_BOX_PASS']
        mechanical.append({'loadcase': name, 'peak_demand_MPa': peak, 'critical_face_id': imax, 'critical_point_mm': demands['cf'][imax].tolist(), 'rows': rows, 'shortest_model_box_pass_cure_min': min(passed) if passed else None, 'physical_D1_release': 'UNKNOWN', 'rigorous_scope': 'Frozen FE stress field and declared strength scenario only'})
    fp = OUT / 'FROZEN_PREDICTIONS.json'
    content = {'frozen_utc': p['frozen_utc'], 'PREREG_R2_sha256': sha(ROOT / 'PREREG_R2.json'), 'purpose': 'Prospective matched D1 laboratory comparison; not independent validation', 'already_seen_sources': ['PMC10088447 Table2/3', 'PMC11012777 Table1', 'X1 D1 geometry/resin FE'], 'source_data_used_for_fitting': 'No fit. Source rows used as input; rechecking them is extraction verification, not held-out validation.', 'held_out_physical_measurements_seen': 0, 'mechanical': frozen, 'optical': {'routes': ['ST', 'SP'], 'predicted_surface_TP_file': 'raw/OPTICAL_INTERVALS.json', 'predicted_surface_TP_sha256': None, 'unsupported_prediction': 'UNKNOWN; require direct assay', 'matching': 'Katana HT same aged protocol; TP CIE metric on declared buccal sample sites'}}
    dump(OUT / 'raw/OPTICAL_INTERVALS.json', optical)
    content['optical']['predicted_surface_TP_sha256'] = sha(OUT / 'raw/OPTICAL_INTERVALS.json')
    if OUT != ROOT:
        assert json.loads((ROOT / 'FROZEN_PREDICTIONS.json').read_text()) == content, 'Replay prediction drift'
    if fp.exists():
        assert json.loads(fp.read_text()) == content, 'Frozen prediction drift'
    else:
        dump(fp, content)
    digest = sha(fp)
    (OUT / 'FROZEN_PREDICTIONS.sha256').write_text(digest + '\n')
    max_i = int(np.flatnonzero(visible & np.isfinite(h))[np.argmax(h[visible & np.isfinite(h)])])
    request = {'claim_type': 'capability', 'frozen_predictions_sha256': digest, 'optical': {'direct_surface_TP_target': 12.0, 'unrestricted_local_law_required_sites_per_route': int(visible.sum()), 'candidate_first_missing_measurement': {'ray_id': max_i, 'point_mm': points[max_i].tolist(), 'path_length_mm': float(h[max_i])}, 'one_site_sufficiency_condition': 'Only if the material is homogeneous and its TP is monotone with these ray lengths, neither verified for a curved crown', 'spectroscopy': 'Measure CIE TP with black/white background and identical illumination, batch and source aging protocol; retain location and repeated readings'}, 'mechanical': {'fixed_load_test_N': 150.0, 'candidate_cures_min': [16.0, 60.0], 'loadcases': ['axial', 'offaxis30'], 'minimum_distinct_matched_conditions': 4, 'protocol': 'Same D1 crown and die, resin batch, layer angle and cure chamber; force/displacement until fracture and fracture-origin image. Keep calibration and validation specimens separate.', 'replicate_count_for_population_reliability': 'UNKNOWN; four conditions do not mean four specimens validate reliability'}, 'consumer': 'Prototype manufacturing/test selection, physical release UNKNOWN', 'resolution_levels': ['PER_POINT optical sites', 'PER_TOOTH fracture measurement'], 'time_scale': 'HANDOVER'}
    dump(OUT / 'MEASUREMENT_REQUEST.json', request)
    result = {'round': 'R2', 'claim_type': 'capability', 'TM02': [{k: r[k] for k in ('route', 'duration_min', 'counts', 'physical_D1_release')} for r in optical], 'TM05': mechanical, 'verification': {'all_model_box_samples_enclosed': bool(all_contained), 'injected_bad_enclosures_rejected': bool(bad_rejected), 'R1_unchanged': sha(ROOT / 'raw/RESULTS_R1.json') == p['parent_results_sha256'], 'physical_release': 'UNKNOWN'}, 'PREREG_sha256': sha(ROOT / 'PREREG_R2.json'), 'FROZEN_PREDICTIONS_sha256': digest, 'external_referent': p['external_referent'], 'measurement_request': 'MEASUREMENT_REQUEST.json', 'numeric_enclosure_scope': p['numerical_limit'], 'cost': {'seconds': time.perf_counter() - tic, 'peak_RSS_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'new_FE_solves': 0, 'new_physical_measurements': 0, 'fit': 0, 'questions': 0}}
    dump(OUT / 'raw/RESULTS_R2.json', result)
    dump(OUT / 'CURRENT_WORK_STATE.json', {'lane': 'X64-manufacturing', 'phase': 'R2_COMPLETE_PENDING_REVIEW', 'last_gate': result['verification'], 'next_operation': 'R3 run the frozen matched optical/mechanical assays; quantify process-linked modulus/contact before claiming physical crown predictions'})
    assert all_contained and bad_rejected, 'Keep failed R2 and freeze changed construction before retry'
    print(json.dumps(result, indent=2))
if __name__ == '__main__':
    main()
