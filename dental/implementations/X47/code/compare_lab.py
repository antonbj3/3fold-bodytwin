"""Score prospective retention observations without fitting predictions.
Synthetic observations are used only in verify.py to show rejection ability.
"""
import csv, json, math, hashlib, argparse
from pathlib import Path
import numpy as np
from scipy.stats import t
from lab_port import allocation, ROOT, sha

def score(rows):
    frozen = json.loads((ROOT / 'FROZEN_PREDICTIONS.json').read_text())
    lock = (ROOT / 'FROZEN_PREDICTIONS.sha256').read_text().split()[0]
    if sha(ROOT / 'FROZEN_PREDICTIONS.json') != lock:
        raise ValueError('FROZEN_PREDICTIONS_HASH_CHANGED')
    if sha(ROOT / 'LAB_ALLOCATION.csv') != frozen['allocation_sha256']:
        raise ValueError('FROZEN_ALLOCATION_HASH_CHANGED')
    if len(rows) != 72 or len({r['specimen_id'] for r in rows}) != 72:
        raise ValueError('ULTIMATE_ENDPOINT_SPECIMEN_REUSE_OR_MISSING')
    by = {r['specimen_id']: r for r in rows}
    for r in allocation():
        actual = by[r['specimen_id']]
        for k in r:
            if str(actual[k]) != str(r[k]):
                raise ValueError('ALLOCATION_MISMATCH:' + k)
    if all((not str(r.get('force_N', '')).strip() for r in rows)):
        return {'status': 'PENDING_MEASUREMENT', 'physical_measurements': 0, 'prediction_hash': lock}
    for r in rows:
        try:
            f = float(r['force_N'])
        except (ValueError, KeyError):
            raise ValueError('MISSING_FORCE')
        if not math.isfinite(f) or f < 0:
            raise ValueError('INVALID_FORCE_N')
        if r['ultimate_endpoint'] == 'retention':
            if r.get('force_unit') != 'N':
                raise ValueError('WRONG_FORCE_UNIT')
            required = ['material_batch', 'cement_product', 'cement_batch', 'surface_treatment', 'die_material_batch', 'actual_TC_cycles', 'post_storage_hours', 'force_trace_sha256', 'image_sha256']
            if any((not str(r.get(k, '')).strip() for k in required)):
                raise ValueError('MISSING_SETUP_OR_TRACE')
            if float(r.get('measured_pull_axis_deg', math.nan)) != 0:
                raise ValueError('NONAXIAL_RETENTION')
            if float(r.get('pull_speed_mm_min', math.nan)) != 1:
                raise ValueError('WRONG_PULL_SPEED')
            if float(r['actual_TC_cycles']) != float(r['thermal_cycles']):
                raise ValueError('WRONG_AGING_DOSE')
            if r['age_state'] == 'TC10000' and any((float(r.get(k, math.nan)) != v for (k, v) in [('TC_low_C', 5), ('TC_high_C', 55), ('dwell_s', 30), ('transfer_s', 5)])):
                raise ValueError('WRONG_THERMAL_PROTOCOL')
            if r['failure_mode'] not in ['cement_debonding', 'premature_decementation', 'crown_fracture', 'die_fracture', 'grip_failure', 'unknown']:
                raise ValueError('WRONG_ENDPOINT_MODE')
            if r['failure_mode'] == 'premature_decementation' and f != 0:
                raise ValueError('PREMATURE_LOSS_MUST_BE_ZERO_TERMINAL_FORCE')
    out = []
    for d in ['D1', 'M1', 'M2']:
        groups = [[r for r in rows if r['design_id'] == d and r['ultimate_endpoint'] == 'retention' and (r['age_state'] == s)] for s in ['sham37', 'TC10000']]
        flat = sum(groups, [])
        if any((r['failure_mode'] not in ['cement_debonding', 'premature_decementation'] for r in flat)):
            out.append({'design_id': d, 'status': 'UNKNOWN_COMPETING_FAILURE_OR_MODE', 'counts': [len(g) for g in groups]})
            continue
        if any((len({r[k] for r in flat}) != 1 for k in ['material_batch', 'cement_product', 'cement_batch', 'surface_treatment', 'die_material_batch'])):
            raise ValueError('UNMATCHED_MATERIAL_OR_CEMENT_BATCH')
        hours = [np.mean([float(r['post_storage_hours']) for r in g]) for g in groups]
        if abs(hours[0] - hours[1]) > 0.1:
            raise ValueError('UNMATCHED_SHAM_EXPOSURE_TIME')
        values = [np.array([float(r['force_N']) for r in g]) for g in groups]
        means = [float(a.mean()) for a in values]
        sd = [float(a.std(ddof=1)) for a in values]
        entry = {'design_id': d, 'mean_N': dict(zip(['sham', 'aged'], means)), 'sd_N': dict(zip(['sham', 'aged'], sd)), 'n': [len(x) for x in values], 'premature_losses': [sum((r['failure_mode'] == 'premature_decementation' for r in g)) for g in groups], 'resolution': 'POPULATION summaries of PER_TOOTH endpoint'}
        if any((s == 0 for s in sd)):
            entry.update(status='UNKNOWN_MEASUREMENT_VARIANCE', observed_ratio=means[1] / means[0] if means[0] > 0 else None, interval=None)
        elif any((m <= 0 for m in means)) or any(((a == 0).any() for a in values)):
            entry.update(status='UNKNOWN_INTERVAL_ZERO_OR_CENSORED', observed_ratio=means[1] / means[0] if means[0] > 0 else None, interval=None)
        else:
            vv = [(s / m) ** 2 / len(a) for (s, m, a) in zip(sd, means, values)]
            var = sum(vv)
            df = var ** 2 / sum((v * v / (len(a) - 1) for (v, a) in zip(vv, values))) if var > 0 else math.inf
            h = t.ppf(0.975, df) * math.sqrt(var) if var > 0 else 0
            mu = math.log(means[1] / means[0])
            ci = [math.exp(mu - h), math.exp(mu + h)]
            (lo, hi) = frozen['practice_aging_null_window']
            status = 'REJECT_AGING_NULL' if ci[1] < lo or ci[0] > hi else 'WITHIN_NULL_TOLERANCE' if lo <= ci[0] and ci[1] <= hi else 'UNKNOWN_OVERLAP'
            entry.update(observed_ratio=math.exp(mu), ratio_95pct_CI_approx=ci, welch_df=df, status=status, interval_caveat='Small n=4 per state, approximate independent-group delta interval; no rigorous enclosure. Single fitted batch is calibration, not validation of transfer.')
        out.append(entry)
    return {'status': 'SCORED', 'prediction_hash': lock, 'designs': out, 'candidate_absolute_force_prediction': 'UNKNOWN', 'candidate_ratio_prediction': 'UNKNOWN', 'physical_measurements': 72, 'specimen_reuse': False}

def main():
    p = argparse.ArgumentParser()
    p.add_argument('csv_path')
    p.add_argument('--output', default='LAB_COMPARISON.json')
    args = p.parse_args()
    with open(args.csv_path) as f:
        rows = list(csv.DictReader(f))
    out = score(rows)
    Path(args.output).write_text(json.dumps(out, indent=2, allow_nan=False) + '\n')
    print(out['status'])
if __name__ == '__main__':
    main()
