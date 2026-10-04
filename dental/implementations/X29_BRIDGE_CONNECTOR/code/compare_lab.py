"""Accept future laboratory data without fitting/replacing a frozen prediction."""
import csv, json, sys
import numpy as np
from scipy.stats import t
from common import R, read, dump, sha

def compare(path, output):
    frozen = read('FROZEN_PREDICTIONS.json')
    frozen_hash = sha(R / 'FROZEN_PREDICTIONS.json')
    rows = list(csv.DictReader(open(path)))
    out = []
    for pred in frozen['paired_contrasts']:
        angle = 0 if pred['load_case'] == 'axial' else 30
        a = [r for r in rows if r['design'] == 'A9' and float(r['angle_deg']) == angle]
        b = [r for r in rows if r['design'] == 'A16' and float(r['angle_deg']) == angle]
        reason = []
        if not frozen['numerical_ratio_gate_pass']:
            reason.append('paired_numerical_gate_failed')
        if len(a) < 2 or len(b) < 2:
            reason.append('too_few_measurements')
        if a and any((r['origin_confirmed'].lower() != 'true' or r['origin'] != 'pontic_contact' for r in a)):
            reason.append('reference_pontic_origin_not_confirmed_per_specimen')
        if any((not r[k] for r in a + b for k in ['batch_id', 'support_id', 'aging_protocol', 'origin'])):
            reason.append('missing_batch_support_aging_or_origin')
        for k in ['batch_id', 'support_id', 'aging_protocol']:
            if len({r[k] for r in a + b}) > 1:
                reason.append('unmatched_' + k)
        if reason:
            out.append(dict(angle_deg=angle, status='UNKNOWN', reasons=reason))
            continue
        A = np.array([float(r['force_N']) for r in a])
        B = np.array([float(r['force_N']) for r in b])
        if np.any(A <= 0) or np.any(B <= 0):
            raise ValueError('Forces must be positiveN')
        va = (A.std(ddof=1) / A.mean()) ** 2 / len(A)
        vb = (B.std(ddof=1) / B.mean()) ** 2 / len(B)
        se = np.sqrt(va + vb)
        df = (va + vb) ** 2 / (va ** 2 / (len(A) - 1) + vb ** 2 / (len(B) - 1)) if se > 0 else 1000
        rat = B.mean() / A.mean()
        ci = [rat * np.exp(-t.ppf(0.975, df) * se), rat * np.exp(t.ppf(0.975, df) * se)]
        (lo, hi) = pred['contact_origin_conditional_capacity_ratio_interval']
        tol = frozen['operational_log_ratio_tolerance']
        window = [lo * np.exp(-tol), hi * np.exp(tol)]
        out.append(dict(angle_deg=angle, status='REJECTED_PROXY_MODEL' if ci[1] < window[0] or ci[0] > window[1] else 'NOT_REJECTED_NOT_VALIDATED', mean_force_ratio=rat, ratio_CI95=ci, frozen_operational_window=window, n_reference=len(A), n_modified=len(B), resolution='POPULATION', new_modes=sorted({r['origin'] for r in b})))
    report = dict(frozen_sha256=frozen_hash, measurements_sha256=sha(path), fit_performed=False, comparisons=out)
    dump(output, report)
    if sha(R / 'FROZEN_PREDICTIONS.json') != frozen_hash:
        raise ValueError('Prediction mutation')
    return report
if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise SystemExit('usage: compare_lab.py measurements.csv output.json')
    print(json.dumps(compare(sys.argv[1], sys.argv[2]), indent=2))
