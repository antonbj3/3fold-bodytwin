"""Consume physical force observations against locked predictions; never fit or rewrite them."""
from common import *
import csv, argparse

def compare(rows, frozen):
    output = []
    for angle in [0, 30]:
        rr = [r for r in rows if int(r['angle_deg']) == angle]
        base = [r for r in rr if r['design'] == 'reference']
        test = [r for r in rr if r['design'] == 'surface_tie_optimized']
        record = dict(angle_deg=angle, resolution='POPULATION', claim_type='capability', physical_calibrated_force05_N=None)
        if len(base) < 6 or len(test) < 6:
            record.update(status='UNKNOWN', reason='need at least six positive numeric fracture endpoints per group')
            output.append(record)
            continue
        fields = ['material_batch', 'die_batch', 'cement_batch', 'footprint_id']
        if any((len({r.get(f, '') for r in base + test}) != 1 or next(iter({r.get(f, '') for r in base + test})) in ['', 'UNKNOWN'] for f in fields)):
            record.update(status='UNKNOWN', reason='batch/support/footprint missing or unmatched')
            output.append(record)
            continue
        if any((r.get('fracture_origin') not in ['ceramic_bulk', 'intaglio'] for r in base + test)):
            record.update(status='REFUTED_OR_UNKNOWN_MECHANISM', reason='requires confirmed ceramic tensile fracture; contact cone, die fracture and unknown origin do not validate volume-flaw model')
            output.append(record)
            continue
        a = np.array([float(r['force_N']) for r in base])
        b = np.array([float(r['force_N']) for r in test])
        if not np.isfinite(a).all() or not np.isfinite(b).all() or min(a.min(), b.min()) <= 0:
            raise ValueError('nonpositive/nonfinite measured force')
        ratio = float(b.mean() / a.mean())
        expected = frozen['contrasts'][str(angle)]['optimized_ratio']
        error = abs(np.log(ratio / expected))
        rng = np.random.default_rng(20261003)
        boot = b[rng.integers(0, len(b), (2000, len(b)))].mean(1) / a[rng.integers(0, len(a), (2000, len(a)))].mean(1)
        record.update(status='WITHIN_OPERATIONAL_WINDOW' if error <= 0.3 else 'REFUTED_FIXED_PREDICTION', observed_arithmetic_mean_ratio=ratio, predicted_conditional_ratio=expected, log_error=float(error), fixed_log_tolerance=0.3, bootstrap95_ratio_interval=np.quantile(boot, [0.025, 0.975]), n_reference=len(a), n_optimized=len(b), scope='Mean-ratio check with common-shape Weibull closure; not calibrated 5% force quantile, robustness certificate or clinical survival')
        output.append(record)
    return output

def run():
    ap = argparse.ArgumentParser()
    ap.add_argument('csv')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    lock = read(ROOT / 'FROZEN_PREDICTIONS.json')
    p = ROOT / 'FROZEN_SURFACE_TIE_PREDICTIONS.json'
    if sha(p) != lock['surface_tie_prediction_sha256']:
        raise ValueError('prediction hash drift')
    rows = list(csv.DictReader(open(args.csv)))
    dump(args.output, dict(prediction_sha256=sha(p), input_sha256=sha(args.csv), comparisons=compare(rows, read(p)), fits_performed=0))
if __name__ == '__main__':
    run()
