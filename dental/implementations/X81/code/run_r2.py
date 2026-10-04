import datetime, hashlib, itertools, math
import numpy as np
from scipy.stats import t
from common import read, save, state, verify_freeze, ROOT

def welch_hw(s1, s2, n1, n2, family=1):
    v1 = s1 * s1 / n1
    v2 = s2 * s2 / n2
    df = (v1 + v2) ** 2 / (v1 * v1 / (n1 - 1) + v2 * v2 / (n2 - 1))
    return t.ppf(1 - 0.05 / (2 * family), df) * np.sqrt(v1 + v2)

def allocation(delta, s1, s2, bias, family, limit=100):
    ns = np.arange(2, limit + 1)
    (n1, n2) = np.meshgrid(ns, ns, indexing='ij')
    hw = welch_hw(s1, s2, n1, n2, family) + bias + 0.01
    valid = hw < delta
    if not valid.any():
        return {'n': None, 'reason': 'BIAS_EXCEEDS_CONTRAST' if bias + 0.01 >= delta else 'NO_PLAN_WITHIN_FROZEN_SEARCH_LIMIT', 'candidate_count': int(hw.size)}
    (rr, cc) = np.where(valid)
    total = n1[rr, cc] + n2[rr, cc]
    vari = s1 * s1 / n1[rr, cc] + s2 * s2 / n2[rr, cc]
    z = int(np.lexsort((n2[rr, cc], n1[rr, cc], vari, total))[0])
    (i, j) = (int(rr[z]), int(cc[z]))
    (a, b) = (int(n1[i, j]), int(n2[i, j]))
    own = float(welch_hw(s1, s2, a, b, family) + bias + 0.01)
    v1 = s1 * s1 / a
    v2 = s2 * s2 / b
    df = (v1 + v2) ** 2 / (v1 * v1 / (a - 1) + v2 * v2 / (b - 1))
    control = float(t.isf(0.05 / (2 * family), df) * math.sqrt(v1 + v2) + bias + 0.01)
    smaller_valid = valid & (n1 + n2 < a + b)
    return {'n': [a, b], 'total_specimens': a + b, 'projected_lower': delta - own, 'halfwidth': own, 'candidate_count': int(hw.size), 'smaller_total_feasible_count': int(smaller_valid.sum()), 'control_difference': abs(own - control), 'sampling_model': 'Welch plug-in; projected interval, not power/physical bound'}

def paired_n(delta, s1, s2, family, limit, rounding):
    ns = np.arange(2, limit + 1)
    hw = t.ppf(1 - 0.05 / (2 * family), ns - 1) * (s1 + s2) / np.sqrt(ns) + rounding
    good = np.where(hw < abs(delta))[0]
    if not len(good):
        return {'n': None, 'reason': 'NO_PLAN_WITHIN_FROZEN_LIMIT'}
    i = int(good[0])
    n = int(ns[i])
    control = t.isf(0.05 / (2 * family), n - 1) * (s1 + s2) / math.sqrt(n) + rounding
    return {'n': n, 'projected_signed_interval': [delta - float(hw[i]), delta + float(hw[i])], 'worst_covariance_SD': s1 + s2, 'control_difference': abs(float(hw[i]) - float(control)), 'smaller_n_feasible_count': int((hw[:i] < abs(delta)).sum()), 'scope': 'Unknown covariance handled by worst rho=-1; plug-in normal sampling only'}

def interp_control(x, y, q):
    val = np.zeros_like(q, dtype=float)
    for i in range(len(x)):
        w = np.ones_like(q, dtype=float)
        for j in range(len(x)):
            if i != j:
                w *= (q - x[j]) / (x[i] - x[j])
        val += y[i] * w
    return val

def calibration_design():
    baseline = read('cbct_baseline')
    ref = next((q for q in baseline if q['device'] == 'Varian'))
    mats = [x['material'] for x in ref['inserts']]
    h = np.array([x['gray_median'] for x in ref['inserts']])
    devices = [q for q in baseline if q['device'] != 'Varian']
    candidates = []
    for k in [2, 3]:
        for inds in itertools.combinations(range(6), k):
            vals = []
            fits = {}
            monotonic = True
            parity = []
            for q in devices:
                g = np.array([x['gray_median'] for x in q['inserts']])
                center = (g.min() + g.max()) / 2
                scale = (g.max() - g.min()) / 2
                x = (g - center) / scale
                sel = list(inds)
                if len(set(x[sel])) < k:
                    monotonic = False
                    continue
                c = np.linalg.solve(np.vander(x[sel], N=k, increasing=True), h[sel])
                pred = np.polynomial.polynomial.polyval(x, c)
                derivative = [c[1]] if k == 2 else [c[1] - 2 * c[2], c[1] + 2 * c[2]]
                monotonic = monotonic and min(derivative) >= 0
                vals.append(float(abs(pred - h).max()))
                parity.append(float(abs(pred - interp_control(x[sel], h[sel], x)).max()))
                fits[q['device']] = {'coefficients': c.tolist(), 'center': float(center), 'scale': float(scale), 'signed_residual': (pred - h).tolist()}
            candidates.append({'anchors': [mats[i] for i in inds], 'k': k, 'indices': list(inds), 'max_baseline_error_HU_ref': max(vals) if vals else None, 'monotonic': bool(monotonic), 'fits': fits, 'control_max_difference': max(parity) if parity else None})
    valid = [x for x in candidates if x['monotonic']]
    best = min(valid, key=lambda x: (x['max_baseline_error_HU_ref'], x['k'], x['indices']))
    predictions = {'timestamp_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'claim_type': 'capability', 'training': 'Six insert medians per baseline scanner; all used in design discovery cost', 'validation': '18 already published acquisitions, retrospective replay; not unseen prospective experiment', 'selected': best, 'prereg_sha256': hashlib.sha256((ROOT / 'PREREG_R2.json').read_bytes()).hexdigest()}
    p = ROOT / 'FROZEN_PREDICTIONS_R2.json'
    if p.exists():
        old = verify_freeze('FROZEN_PREDICTIONS_R2.json')
        assert old['selected'] == best, 'Prediction drift'
    else:
        save(p.name, predictions)
        (ROOT / (p.name + '.sha256')).write_text(hashlib.sha256(p.read_bytes()).hexdigest() + '\n')
    held = read('cbct_held')
    rows = []
    for q in held:
        if q['device'] == 'Varian':
            continue
        g = np.array([x['gray_median'] for x in q['inserts']])
        fit = best['fits'][q['device']]
        x = (g - fit['center']) / fit['scale']
        pred = np.polynomial.polynomial.polyval(x, fit['coefficients'])
        err = pred - h
        isrepeat = 'test-retest2_' in q['zip_member'] or 'test-retest3_' in q['zip_member']
        rows.append({'device': q['device'], 'zip_member': q['zip_member'], 'stratum': 'REPEAT' if isrepeat else 'PROTOCOL_CHANGE', 'signed_error_HU_ref': err.tolist(), 'max_abs_HU_ref': float(abs(err).max()), 'passes40': bool(abs(err).max() <= 40), 'reference': 'Fixed Varian baseline HU_ref; absolute truth UNKNOWN', 'resolution': 'PER_SURFACE_REGION'})
    bounds = []
    for q in devices:
        g = np.array([x['gray_median'] for x in q['inserts']])
        for (i, j) in itertools.combinations(range(6), 2):
            if (g[i] - g[j]) * (h[i] - h[j]) < 0:
                bounds.append({'device': q['device'], 'materials': [mats[i], mats[j]], 'lower_bound_HU_ref': abs(float(h[i] - h[j])) / 2, 'kind': 'Exact minimax order bound on fixed observed medians; not physical uncertainty bound'})
    return {'all_candidates': candidates, 'selected': best, 'held_rows': rows, 'held_scans_evaluated': len(rows), 'published_held_scans_not_evaluated_reference_device': sum((q['device'] == 'Varian' for q in held)), 'held_max_HU_ref': max((q['max_abs_HU_ref'] for q in rows)), 'monotone_lower_bounds': bounds, 'outcome': 'SUFFICIENT_40HU_CALIBRATION' if all((q['passes40'] for q in rows)) else 'CALIBRATION_DESIGN_REFUTED_ON_EXISTING_HOLDOUT', 'new_exposures': 0, 'discovery_existing_ROIs': 24, 'selected_ROIs_per_scan': best['k'], 'cost_note': 'Anchor count is not scan count: all six inserts share one phantom image; cannot claim exposure saving', 'rigorous_physical_enclosure': 'MISSING'}

def thermal_design():
    rows = read('laser')
    pairs = []
    for freq in [2, 4, 6]:
        for energy in [250, 300, 350]:
            a = next((q for q in rows if q['frequency_Hz'] == freq and q['energy_mJ'] == energy and (not q['water'])))
            b = next((q for q in rows if q['frequency_Hz'] == freq and q['energy_mJ'] == energy and q['water']))
            delta = a['reported_rise_C'] - b['reported_rise_C']
            s1 = a['sample_SD_C']
            s2 = b['sample_SD_C']
            hw = float(welch_hw(s1, s2, 4, 4, 9)) + 0.01
            pairs.append({'frequency_Hz': freq, 'energy_mJ': energy, 'dry_minus_cooled_degC': delta, 'sample_SD_dry_degC': s1, 'sample_SD_cooled_degC': s2, 'current_n_per_group': 4, 'n_status': 'INFERRED_BALANCED_ALLOCATION_NOT_CONFIRMED', 'current_simultaneous95_interval_degC': [delta - hw, delta + hw], 'current_with_drift95_interval_degC': [delta - hw - 0.4, delta + hw + 0.4], 'minimum_no_drift': allocation(delta, s1, s2, 0, 9), 'minimum_drift_bound_0p4': allocation(delta, s1, s2, 0.4, 9), 'resolution': 'POPULATION', 'timescale': 'SIMULTANEOUS', 'locator': a['locator']})
    return {'rows': pairs, 'positive_current_no_drift': sum((x['current_simultaneous95_interval_degC'][0] > 0 for x in pairs)), 'positive_current_with_drift': sum((x['current_with_drift95_interval_degC'][0] > 0 for x in pairs)), 'optimally_allocated_specimens_no_drift': sum((x['minimum_no_drift'].get('total_specimens', 0) for x in pairs)), 'uniform_allocation_specimens_no_drift': 18 * max((max(x['minimum_no_drift']['n']) for x in pairs if x['minimum_no_drift']['n'])), 'drift_unresolved_count': sum((x['minimum_drift_bound_0p4']['n'] is None for x in pairs)), 'drift_structurally_impossible_count': sum((x['minimum_drift_bound_0p4'].get('reason') == 'BIAS_EXCEEDS_CONTRAST' for x in pairs)), 'drift_search_limit_count': sum((x['minimum_drift_bound_0p4'].get('reason') == 'NO_PLAN_WITHIN_FROZEN_SEARCH_LIMIT' for x in pairs)), 'scope': 'Laser source bovine teeth2mm dentin. Mean contrast only; no individual safety or rotary-preparation transfer', 'accuracy_source': 'Geraldo2005 p184 thermocouple+/-0.1degC. Worst four absolute readings => contrast0.4degC; no gain/contact bound', 'minimal_new_measurement': 'Baseline/reference channel bracketing each matched preparation, synchronized pulp sensor and exact remaining dentin; estimate differential drift separately from specimen scatter', 'rigorous_physical_enclosure': 'MISSING; normal plug-in simultaneous interval closure only'}

def paired_designs():
    replicas = []
    for r in read('replica'):
        delta = float(r['replica_mean_um']) - float(r['ct_mean_um'])
        s1 = float(r['replica_sd_um'])
        s2 = float(r['ct_sd_um'])
        n = r['n_specimens']
        hw = t.ppf(1 - 0.05 / 16, n - 1) * (s1 + s2) / math.sqrt(n) + 0.01
        replicas.append({'system': r['system'], 'region': r['region'], 'difference_um': delta, 'simultaneous95_worst_covariance_interval_um': [delta - float(hw), delta + float(hw)], 'minimum_projected_paired_n': paired_n(delta, s1, s2, 8, 500, 0.01), 'n': n, 'resolution': 'POPULATION', 'consumer_resolution': 'PER_SURFACE_REGION', 'observation_state': r['physical_state'], 'locator': r['source_locator'], 'does_not_identify': 'Sensor-only correction, cement-film agreement and same-state lab bias'})
    p = read('preload')
    delta = p['first']['mean_N'] - p['tenth']['mean_N']
    return {'replica': replicas, 'preload_minimum': paired_n(delta, p['first']['sd_N'], p['tenth']['sd_N'], 1, 200, 0.1), 'preload_scope': 'First-vs-tenth contrast only, same coated system and dry25Ncm protocol; projected design not new-system validation', 'minimum_replica_extension': 'Same region/specimen four-state contrasts: dry-seat, PVS seated, PVS removed/reseat, same-state independent film reference; held specimens retained. Crossed method/state avoids confounded correction.'}

def main():
    verify_freeze('PREREG_R2.json')
    out = {'claim_type': 'capability', 'calibration': calibration_design(), 'thermal': thermal_design(), 'paired': paired_designs(), 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'physical_measurement_execution': 'NOT_RUN'}
    save('rounds/R2_results.json', out)
    state('R2_COMPLETE', out['calibration']['outcome'], 'Freeze R3: complete R3 edge/group table, decision-specific lab frontier and adversarial controls')
    (ROOT / 'rounds/HANDOFF_R2.md').write_text('R2: minimum integer projected contrast designs delivered, with covariance/drift held explicit. Calibration selected on baseline then rejected on held protocols. Low thermal contrasts cannot overcome a drift budget by repetitions. Next operation: retain decision-specific uncertainty and finish all edge/group routing with quote-dependent lab frontier; do not call projected intervals physical certificates.\n')
    print('R2 anchors', out['calibration']['selected']['anchors'], 'held max', out['calibration']['held_max_HU_ref'])
    print('Thermal simultaneous positive', out['thermal']['positive_current_no_drift'], out['thermal']['positive_current_with_drift'], 'allocation', out['thermal']['optimally_allocated_specimens_no_drift'], 'uniform', out['thermal']['uniform_allocation_specimens_no_drift'], 'drift unresolved', out['thermal']['drift_unresolved_count'])
    print('Preload projected pairs', out['paired']['preload_minimum']['n'])
if __name__ == '__main__':
    main()
