"""Paired closure acquisition with common-mode cancellation before integration."""
from dental_release.paths import expand as _release_expand
import hashlib, json, resource, time, warnings
from pathlib import Path
import numpy as np
from scipy.optimize import linprog, OptimizeWarning
from run_r1 import ROOT, model, sha, write
from run_r2 import sample_ramp, work_bounds
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/LANE_XBREAK_HUNT_1_NOSEARCH'))

def paired_work(ts, baseline, shim, common_bound, residual_bound, guard):
    """Consumer sees paired measurements and bounds, never oracle stiffness."""
    integral = float(np.trapz(baseline - shim, ts))
    n = len(ts) - 1
    T = float(ts[-1])
    dt = T / n
    vbase = float(baseline[-1] - baseline[0] + 2 * common_bound + 2 * guard)
    vshim = float(shim[-1] - shim[0] + 2 * common_bound + 2 * residual_bound + 2 * guard)
    err = 0.5 * dt * (max(0.0, vbase) + max(0.0, vshim)) + T * (residual_bound + 2 * guard)
    return (integral - err, integral + err)

def control_work(ts, baseline, shim, common_bound, residual_bound, guard):
    dt = float(ts[1] - ts[0])
    T = float(ts[-1])

    def trapezoid(x):
        return dt * (float(np.sum(x)) - float(x[0] + x[-1]) / 2)
    vb = max(0.0, float(baseline[-1] - baseline[0]) + 2 * (common_bound + guard))
    vs = max(0.0, float(shim[-1] - shim[0]) + 2 * (common_bound + residual_bound + guard))
    radius = dt * (vb + vs) / 2 + T * (residual_bound + 2 * guard)
    center = trapezoid(baseline) - trapezoid(shim)
    return (center - radius, center + radius)

def simplex(upper, T):
    hi = np.minimum(T, np.maximum(0.0, np.asarray(upper)))
    if hi.sum() < T - 1e-05:
        return None
    return (np.maximum(0.0, T - hi.sum() + hi), hi)

def main():
    start = time.perf_counter()
    cpu = time.process_time()
    p = ROOT / 'PREREG_R4.json'
    cfg = json.loads(p.read_text())
    frozen = sha(p)
    write('raw/PREREG_R4_FROZEN.json', {'sha256': frozen, 'file': str(p)})
    worlds = json.loads((ROOT / 'raw/R1_WORLDS.json').read_text())
    obs = json.loads((ROOT / 'raw/R1_OBSERVATIONS.json').read_text())
    T = 205.0
    h = cfg['h_mm']
    maxn = max(cfg['force_intervals'])
    common = cfg['common_displacement_error_mm']
    guard = cfg['numerical_displacement_guard_mm']
    targets = list(range(11, 19)) + list(range(21, 29))
    records = []
    cost = []
    manifest = []
    work_coverage = []
    false_reference = []
    lp_calls = 0
    max_lp_diff = 0.0
    warnings.filterwarnings('ignore', category=OptimizeWarning, message='Unrecognized options detected.*')
    DATA.mkdir(parents=True, exist_ok=True)
    for (wi, w) in enumerate(worlds):
        (oracle, teeth, _, _) = model(w['case'], w['scenario'])
        (ts, base, stats) = sample_ramp(oracle, oracle.g, T, maxn)
        cost.append(stats)
        drift = common * (0.6 * np.sin(ts * 0.13 + wi) + 0.4 * np.cos(ts * 0.017))
        base_obs = base + drift
        arrays = {'force_N': ts, 'baseline_d_mm': base, 'common_drift_mm': drift}
        shims = {}
        for tooth in targets:
            direction = oracle.B[teeth.index(tooth)] if tooth in teeth else np.zeros(len(oracle.g))
            (_, d, stats) = sample_ramp(oracle, oracle.g - h * direction, T, maxn)
            cost.append(stats)
            shims[tooth] = d
            arrays[f'shim_{tooth}_d_mm'] = d
        datafile = DATA / f"R4_{w['case']}_{w['scenario']}.npz"
        np.savez_compressed(datafile, **arrays)
        manifest.append({'path': str(datafile), 'sha256': sha(datafile), 'bytes': datafile.stat().st_size, 'uncompressed_array_bytes': sum((a.nbytes for a in arrays.values()))})
        ideal = {t: w['native_V_Nmm'] - next((o['minus_work_Nmm'] for o in obs if o['case'] == w['case'] and o['scenario'] == w['scenario'] and (o['FDI'] == t) and (o['h_mm'] == h))) for t in targets}
        top = set(sorted(targets, key=lambda t: (-w['native_force_N'].get(str(t), 0), t))[:2])
        for n in cfg['force_intervals']:
            ix = np.arange(0, maxn + 1, maxn // n)
            tt = ts[ix]
            for eta in cfg['residual_differential_error_mm']:
                uppers = []
                independent_uppers = []
                details = []
                strongdiff = 0.0
                for tooth in targets:
                    residual = eta * np.sin(ts * 0.311 + tooth)
                    shim_obs = shims[tooth] + drift + residual
                    bounds = paired_work(tt, base_obs[ix], shim_obs[ix], common, eta, guard)
                    conventional = control_work(tt, base_obs[ix], shim_obs[ix], common, eta, guard)
                    strongdiff = max(strongdiff, max((abs(a - b) for (a, b) in zip(bounds, conventional))))
                    work_coverage.append(bounds[0] - 1e-05 <= ideal[tooth] <= bounds[1] + 1e-05)
                    uppers.append(bounds[1] / h)
                    wb = work_bounds(tt, base_obs[ix], common, guard)
                    ws = work_bounds(tt, shim_obs[ix], common + eta, guard)
                    independent_uppers.append((wb[1] - ws[0]) / h)
                    details.append({'work_difference_Nmm': list(bounds), 'heldout_ideal_work_difference_Nmm': ideal[tooth], 'control_work_difference_Nmm': list(conventional)})
                answer = simplex(uppers, T)
                ctrl = simplex(independent_uppers, T)
                if answer is None or ctrl is None:
                    raise ArithmeticError('Allowed noise caused infeasible simplex')
                (lo, hi) = answer
                (clo, chi) = ctrl
                for (i, tooth) in enumerate(targets):
                    q = w['native_force_N'].get(str(tooth), 0.0)
                    if n == maxn and eta == 1e-06:
                        objective = np.eye(len(targets))[i]
                        a = linprog(objective, A_eq=np.ones((1, len(targets))), b_eq=[T], bounds=list(zip(np.zeros(len(targets)), hi)), method='highs', options={'threads': 4})
                        b = linprog(-objective, A_eq=np.ones((1, len(targets))), b_eq=[T], bounds=list(zip(np.zeros(len(targets)), hi)), method='highs', options={'threads': 4})
                        lp_calls += 2
                        if not a.success or not b.success:
                            raise ArithmeticError('LP failure')
                        max_lp_diff = max(max_lp_diff, float(max(abs(a.fun - lo[i]), abs(-b.fun - hi[i]))))
                    records.append({'case': w['case'], 'scenario': w['scenario'], 'FDI': tooth, 'intervals': n, 'residual_error_mm': eta, 'facit_N': q, 'lower_N': float(lo[i]), 'upper_N': float(hi[i]), 'width_pp': float(100 * (hi[i] - lo[i]) / T), 'covered': bool(lo[i] - 1e-05 <= q <= hi[i] + 1e-05), 'top_two': tooth in top, 'independent_error_lower_N': float(clo[i]), 'independent_error_upper_N': float(chi[i]), 'independent_error_width_pp': float(100 * (chi[i] - clo[i]) / T), 'same_information_work_control_max_difference_Nmm': strongdiff, **details[i]})
                if n == maxn and eta == 1e-06:
                    tooth = sorted(top, key=lambda t: (-w['native_force_N'].get(str(t), 0), t))[0]
                    i = targets.index(tooth)
                    wrong = paired_work(tt, base_obs[ix], (shims[tooth] + drift + 0.002)[ix], common, eta, guard)
                    altered = list(uppers)
                    altered[i] = wrong[1] / h
                    bad = simplex(altered, T)
                    false_reference.append({'case': w['case'], 'scenario': w['scenario'], 'FDI': tooth, 'unshared_offset_mm': 0.002, 'asserted_residual_error_mm': eta, 'false_cap_excludes_facit': wrong[1] / h + 1e-05 < w['native_force_N'].get(str(tooth), 0), 'simplex_guard_status': 'UNKNOWN_INCONSISTENT_OBSERVATION' if bad is None else 'FEASIBLE_BUT_UNVALIDATED_CALIBRATION', 'required_physical_status_without_calibration': 'UNKNOWN'})
        write('raw/R4_PROGRESS.json', {'worlds_completed': wi + 1, 'cost': cost, 'records': records, 'data_manifest': manifest})
        print(w['case'], w['scenario'], 'paired ramp panel completed', flush=True)
    summary = []
    for eta in cfg['residual_differential_error_mm']:
        for n in cfg['force_intervals']:
            rr = [r for r in records if r['residual_error_mm'] == eta and r['intervals'] == n]
            toprr = [r for r in rr if r['top_two']]
            summary.append({'residual_error_mm': eta, 'intervals': n, 'n': len(rr), 'covered': sum((r['covered'] for r in rr)), 'fraction_all_le5pp': float(np.mean([r['width_pp'] <= 5 for r in rr])), 'fraction_top_two_le5pp': float(np.mean([r['width_pp'] <= 5 for r in toprr])), 'median_top_two_width_pp': float(np.median([r['width_pp'] for r in toprr])), 'max_top_two_width_pp': max((r['width_pp'] for r in toprr)), 'independent_error_fraction_all_le5pp': float(np.mean([r['independent_error_width_pp'] <= 5 for r in rr])), 'independent_error_fraction_top_two_le5pp': float(np.mean([r['independent_error_width_pp'] <= 5 for r in toprr])), 'independent_error_median_top_two_width_pp': float(np.median([r['independent_error_width_pp'] for r in toprr]))})

    def useful(s):
        return s['fraction_all_le5pp'] >= 0.8 and s['fraction_top_two_le5pp'] >= 0.8

    def gate(eta):
        return any((useful(s) for s in summary if s['residual_error_mm'] == eta))
    out = {'round': 'R4', 'claim_type': 'information_link', 'status': 'PENDING_INDEPENDENT_REVIEW', 'prereg_sha256': frozen, 'coverage_gate': all((r['covered'] for r in records)), 'work_difference_coverage_gate': all(work_coverage), 'one_nm_instrument_gate': gate(1e-06), 'ten_nm_instrument_gate': gate(1e-05), 'reference_value_gate': any((useful(s) and s['independent_error_fraction_top_two_le5pp'] < 0.8 for s in summary)), 'false_reference_gate': any((r['false_cap_excludes_facit'] for r in false_reference)), 'false_reference_probe': false_reference, 'same_information_LP_calls': lp_calls, 'max_LP_difference_N': max_lp_diff, 'cost': {'ramp_count': len(cost), 'force_readings_at_max_budget': sum((c['readings'] for c in cost)), 'slow_QPs': sum((c['slow_QP'] for c in cost)), 'fast_KKT': sum((c['fast_KKT'] for c in cost)), 'sampling_seconds': sum((c['seconds'] for c in cost)), 'wall_seconds': time.perf_counter() - start, 'CPU_seconds': time.process_time() - cpu, 'max_RSS_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'physical_acquisition_time': 'UNKNOWN; NOT_RUN', 'upstream_measured_round_seconds': [4.156158004887402, 3.5929280628915876, 2.4491713501047343]}, 'data_manifest': manifest, 'summaries': summary, 'records': records, 'physical_validation': 'UNKNOWN; residual nm errors are hypothetical instrument requirements; common reference, gap step and biology NOT_MEASURED', 'numerical_scope': 'Floating KKT checked; declared 1e-7 mm guard is not a formal numerical proof', 'algorithm_novelty': 'NOT_CLAIMED; conventional equal-information work control and LP executed'}
    write('rounds/R4/results.json', out)
    assert sha(p) == frozen
    print(json.dumps({k: out[k] for k in ['coverage_gate', 'one_nm_instrument_gate', 'ten_nm_instrument_gate', 'reference_value_gate', 'cost', 'summaries']}, indent=2))
if __name__ == '__main__':
    main()
