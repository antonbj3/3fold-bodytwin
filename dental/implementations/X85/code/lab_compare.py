"""Validate frozen physical predictions without fitting. Units/identity must match.

Release repair: reject nonfinite, unordered or malformed prediction intervals
before comparing calibrated observations; no numerical threshold changed.
"""
from common import *
import argparse

def compare(frozen, measurement):
    missing = [k for k in ['case', 'geometry_sha256', 'frame', 'unit', 'matched_load_rate_history', 'reaction_N', 'reaction_error_N'] if k not in measurement]
    if missing:
        return dict(status='UNKNOWN', reason='Missing measurement fields', missing=missing)
    if measurement['unit'] != 'N':
        return dict(status='UNKNOWN', reason='Absolute regional calibrated reactions required; relative T-Scan is a separate comparison')
    for key in ['case', 'geometry_sha256', 'frame']:
        if measurement[key] != frozen[key]:
            raise ValueError('Heldout identity mismatch: ' + key)
    if measurement['matched_load_rate_history'] is not True:
        return dict(status='UNKNOWN', reason='Load/rate/history mismatch')
    y = np.asarray(measurement['reaction_N'], float)
    e = np.asarray(measurement['reaction_error_N'], float)
    iv = np.asarray(frozen['force_interval_N'], float)
    if iv.ndim != 2 or iv.shape[1] != 2 or iv.shape[0] == 0 or not np.isfinite(iv).all() or np.any(iv[:, 0] > iv[:, 1]):
        raise ValueError('Frozen force intervals must be a nonempty finite ordered (n, 2) array')
    if y.shape != (len(iv),) or e.shape != y.shape or (not np.isfinite(np.r_[y, e]).all()) or (np.min(e) < 0):
        raise ValueError('Invalid force measurement/error')
    passall = bool(np.all(y - e >= iv[:, 0]) and np.all(y + e <= iv[:, 1]))
    return dict(status='PASS_FROZEN_OBSERVATION_GATE' if passall else 'FAIL_FROZEN_OBSERVATION_GATE', measured_interval_N=np.c_[y - e, y + e], frozen_interval_N=iv, fit_performed=False, physical_claim='One heldout check does not prove arbitrary nonlinear/historical transfer')

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('frozen')
    ap.add_argument('measurement')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    r = compare(read(args.frozen), read(args.measurement))
    dump(args.out, r)
    print(r['status'])
if __name__ == '__main__':
    main()
