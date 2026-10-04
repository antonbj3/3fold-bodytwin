"""Compare future independent geometry measurements to immutable digital targets.

python3 code/compare_metrology.py path/to/measurements.json
No measurements are fabricated or assumed. This does not assess anatomy or health.
"""
from common_local import *
import argparse

def compare(data):
    f = P / 'FROZEN_PREDICTIONS.json'
    assert sha(f) == f.with_suffix('.json.sha256').read_text().strip()
    frozen = load(f)
    if not data.get('measurements'):
        return {'status': 'UNKNOWN_NO_INDEPENDENT_MEASUREMENT', 'compared': 0}
    if data.get('kind') != 'independent_measurement':
        raise ValueError('Independent measurement required; fixture or digital replay is not physical validation')
    if data.get('frozen_prediction_sha256') != sha(f):
        raise ValueError('Frozen prediction hash mismatch')
    if data.get('units') != 'mm':
        raise ValueError('Explicit mm units required')
    if not data.get('instrument_record') or not data.get('registration_record'):
        raise ValueError('Instrument and registration evidence required')
    pred = {(r['case'], r['fdi'], r['advance_mm']): r for r in frozen['targets']}
    out = []
    for m in data['measurements']:
        k = (m['case'], m['fdi'], m['advance_mm'])
        if k not in pred:
            raise ValueError('Query absent from frozen targets')
        r = pred[k]
        if m['L_mm'] != r['L_mm'] or m['D_mm'] != r['D_mm']:
            raise ValueError('Dimensions differ from frozen specimen')
        value = float(m['distance_mm'])
        err = float(m['absolute_measurement_bound_mm'])
        if not np.isfinite(value + err) or value < 0 or err < 0:
            raise ValueError('Distance and absolute bound must be finite nonnegative')
        tol = frozen['absolute_acceptance_tolerance_mm']
        discrepancy = max(0.0, r['gap_lower_mm'] - (value + err), value - err - r['gap_upper_mm'])
        all_within = value - err >= r['gap_lower_mm'] - tol and value + err <= r['gap_upper_mm'] + tol
        out.append(dict(case=k[0], fdi=k[1], advance_mm=k[2], status='FAIL_DISJOINT' if discrepancy > tol else 'PASS_WITHIN_DECLARED_LAB_TOLERANCE' if all_within else 'UNKNOWN_MEASUREMENT_TOO_WIDE', distance_interval_mm=[value - err, value + err], frozen_gap_mm=[r['gap_lower_mm'], r['gap_upper_mm']], minimum_discrepancy_mm=discrepancy))
    return {'status': 'COMPARED_GEOMETRIC_LAB_TARGETS_ONLY', 'compared': len(out), 'rows': out, 'clinical_validation': 'UNKNOWN', 'anatomy_validation': 'UNKNOWN', 'no_automatic_refit': True}
if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input')
    args = ap.parse_args()
    print(json.dumps(compare(load(args.input)), indent=2))
