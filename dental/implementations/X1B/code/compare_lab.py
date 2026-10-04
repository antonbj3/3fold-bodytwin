"""Score a prospective specimen CSV without fitting or modifying predictions."""
import argparse, csv, json, math, statistics, hashlib
from collections import defaultdict, Counter
from pathlib import Path
from common import R, sha, read
INCOMPATIBLE = {'contact_damage', 'cement_debonding', 'die_fracture'}

def compare(frozen, rows):
    if not rows:
        raise ValueError('No specimens')
    designs = {d['design']: d for d in frozen['designs']}
    seen = set()
    out = []
    groups = defaultdict(list)
    for row in rows:
        sid = row['specimen_id']
        if not sid or sid in seen:
            raise ValueError('Empty or duplicate specimen_id')
        seen.add(sid)
        d = designs[row['design_id']]
        angle = float(row['load_angle_deg'])
        if angle not in (0, 30):
            raise ValueError('Unregistered load angle')
        obs = float(row['fracture_force_N'])
        gap = float(row['MG_distance_mean_um'])
        if not math.isfinite(obs) or obs <= 0 or (not math.isfinite(gap)) or (gap < 0):
            raise ValueError('Invalid force or gap')
        mode = row.get('failure_mode', '').strip().lower()
        origin = row.get('fracture_origin', '').strip().lower()
        mechanism = 'INCOMPATIBLE' if mode in INCOMPATIBLE or origin in ('contact_zone', 'cement_interface', 'die') else 'PASS' if mode == 'crown_tensile_fracture' and origin == 'intaglio_tensile_zone' else 'UNKNOWN'
        prediction = d.get('load_predictions', {}).get(str(int(angle)))
        pred = prediction['point_N'] if prediction else d['fracture_load_proxy_axial_N' if angle == 0 else 'fracture_load_proxy_offaxis30_N']
        x = dict(specimen_id=sid, design=row['design_id'], angle_deg=angle, force_observed_N=obs, force_predicted_N=pred, log_force_error=math.log(pred / obs), gap_bias_um=d['gap']['marginal_mean_um'] - gap, failure_mode=mode, fracture_origin=origin, mechanism=mechanism, prediction_validity=prediction.get('validity', 'UNKNOWN') if prediction else 'LEGACY_SCENARIO')
        out.append(x)
        groups[x['design'], angle].append(x)
    strata = []
    for ((design, angle), rs) in sorted(groups.items()):
        n = len(rs)
        counts = Counter((r['mechanism'] for r in rs))
        bad = counts['INCOMPATIBLE'] / n
        mechanism = 'FAIL' if bad > 0.5 else 'UNKNOWN' if counts['UNKNOWN'] or counts['INCOMPATIBLE'] else 'PASS'
        pred = designs[design].get('load_predictions', {}).get(str(int(angle)))
        mean = statistics.mean((r['force_observed_N'] for r in rs))
        interval_status = 'UNKNOWN'
        if pred and pred.get('interval_N') and (n >= pred.get('target_n', 6)):
            interval_status = 'PASS' if pred['interval_N'][0] <= mean <= pred['interval_N'][1] else 'FAIL'
        strata.append(dict(design=design, angle_deg=angle, n=n, force_mean_N=mean, incompatible_fraction=bad, counts=dict(counts), mechanism_gate=mechanism, force_interval_gate=interval_status, interval=pred.get('interval_N') if pred else None))
    mechanism = 'FAIL' if any((s['mechanism_gate'] == 'FAIL' for s in strata)) else 'UNKNOWN' if any((s['mechanism_gate'] == 'UNKNOWN' for s in strata)) else 'PASS'
    error = math.sqrt(statistics.mean((r['log_force_error'] ** 2 for r in out)))
    bias = statistics.mean((r['gap_bias_um'] for r in out))
    force_gate = error <= 0.3
    gap_gate = all((abs(statistics.mean((r['gap_bias_um'] for r in rs))) <= 30 for rs in groups.values()))
    interval_gate = 'FAIL' if any((s['force_interval_gate'] == 'FAIL' for s in strata)) else 'UNKNOWN' if any((s['force_interval_gate'] == 'UNKNOWN' for s in strata)) else 'PASS'
    model_gate = 'PASS' if all((r['prediction_validity'] == 'SUPPORTED' for r in out)) else 'UNKNOWN'
    conclusion = 'FAIL' if not force_gate or not gap_gate or mechanism == 'FAIL' or (interval_gate == 'FAIL') else 'UNKNOWN' if 'UNKNOWN' in (mechanism, interval_gate, model_gate) else 'PASS'
    double_ratios = []
    reference_variance = None
    for target in frozen.get('double_ratio_predictions', []):
        design = target['design']
        ref = target['reference_design']
        keys = [(design, 30), (ref, 30), (design, 0), (ref, 0)]
        if any((k not in groups or len(groups[k]) < target['n_min_per_group'] for k in keys)):
            double_ratios.append(dict(design=design, gate='UNKNOWN', reason='Incomplete matched design-angle groups'))
            continue
        rr = [groups[k] for k in keys]
        means = [statistics.mean((r['force_observed_N'] for r in x)) for x in rr]
        vv = [statistics.variance((r['force_observed_N'] for r in x)) / (len(x) * m * m) for (x, m) in zip(rr, means)]
        q = means[0] / means[1] / (means[2] / means[3])
        logq = math.log(q)
        se = math.sqrt(sum(vv))
        ci = [math.exp(logq - 1.95996398454 * se), math.exp(logq + 1.95996398454 * se)]

        def gate(window):
            return 'FAIL' if ci[1] < window[0] or ci[0] > window[1] else 'PASS' if ci[0] >= window[0] and ci[1] <= window[1] else 'UNKNOWN'
        fegate = gate(target['acceptance_window_Q'])
        ctrlgate = gate(target['control_acceptance_window_Q'])
        localmode = [s['mechanism_gate'] for s in strata if (s['design'], s['angle_deg']) in keys]
        matched_input = [row for row in rows if (row['design_id'], float(row['load_angle_deg'])) in keys]
        protocol = 'PASS'
        for col in ['material_batch', 'cement_batch', 'die_material_batch']:
            vals = {r.get(col, '').strip() for r in matched_input}
            if '' in vals:
                protocol = 'UNKNOWN'
            elif len(vals) > 1:
                protocol = 'FAIL'
                break
        for (col, expected, tol) in [('die_E_MPa', 18000, 1800), ('indenter_diameter_mm', 5, 0.05), ('crosshead_mm_min', 0.5, 0.025)]:
            for row in matched_input:
                try:
                    value = float(row.get(col, ''))
                except ValueError:
                    if protocol != 'FAIL':
                        protocol = 'UNKNOWN'
                    continue
                if not math.isfinite(value) or abs(value - expected) > tol:
                    protocol = 'FAIL'
        local_mechanism = 'FAIL' if 'FAIL' in localmode else 'UNKNOWN' if 'UNKNOWN' in localmode else 'PASS'
        verdict = 'FAIL_MECHANISM' if local_mechanism == 'FAIL' else 'FAIL_PROTOCOL' if protocol == 'FAIL' else 'UNKNOWN' if 'UNKNOWN' in [local_mechanism, protocol] else 'FE_REJECTED' if fegate == 'FAIL' else 'FE_SUPPORTED_CONTROL_REJECTED' if fegate == 'PASS' and ctrlgate == 'FAIL' else 'TIE' if fegate == 'PASS' and ctrlgate == 'PASS' else 'UNKNOWN'
        if local_mechanism == 'PASS' and protocol == 'PASS' and (fegate == 'FAIL') and (ctrlgate == 'FAIL'):
            verdict = 'FAIL_BOTH'
        double_ratios.append(dict(design=design, observed_Q=q, measured_95CI_Q=ci, log_Q_standard_error=se, FE_Q=target['Q_FE'], FE_gate=fegate, separable_control_Q=1.0, control_gate=ctrlgate, mechanism_gate=local_mechanism, protocol_gate=protocol, decision=verdict, CI_kind='delta method on independent arithmetic group means, normal approximation'))
        reference_variance = vv[1] + vv[3]
    return dict(rows=out, strata=strata, log_force_RMSE=error, force_accuracy_gate_pass=force_gate, marginal_bias_um=bias, gap_gate_pass=gap_gate, mechanism_gate=mechanism, interval_gate=interval_gate, model_domain_gate=model_gate, conclusion=conclusion, double_ratio_comparison=double_ratios, shared_reference_logQ_covariance=reference_variance, calibration_performed=False, clinical_inference=False)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('csv')
    ap.add_argument('--predictions', default=str(R / 'FROZEN_PREDICTIONS.json'))
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    p = Path(a.predictions)
    if sha(p) != p.with_suffix('.sha256').read_text().strip():
        raise ValueError('Frozen prediction hash drift')
    result = compare(json.loads(p.read_text()), list(csv.DictReader(open(a.csv))))
    result.update(prediction_sha256=sha(p), source_csv_sha256=sha(a.csv), prediction_artifact_unchanged=True)
    with Path(a.output).open('x') as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write('\n')
    print(result['conclusion'])
if __name__ == '__main__':
    main()
