"""Minimum physical inputs for X31; fixtures never calibrate a measured tail."""
import argparse, csv, hashlib, importlib.util, json, math
from pathlib import Path
ROOT = Path(__file__).resolve().parent
FIELDS = ['specimen_id', 'fdi', 'measurement_kind', 'unit', 'source_locator', 'ray_frame_locator', 'horn_and_DEJ_verified', 'total_tissue_lower_mm', 'total_tissue_upper_mm', 'enamel_lower_mm', 'enamel_upper_mm', 'reduction_lower_mm', 'reduction_upper_mm']

def local_dentin_interval(row):
    for k in FIELDS:
        if k not in row or str(row[k]) == '':
            raise ValueError('Missing field: ' + k)
    if row['measurement_kind'] != 'INDEPENDENT_PHYSICAL':
        raise ValueError('A fixture is not a physical dentin measurement')
    if row['unit'] != 'mm':
        raise ValueError('Only verified millimeter observations accepted')
    if str(row['horn_and_DEJ_verified']).lower() != 'true':
        raise ValueError('Matched horn and DEJ not established')
    keys = ['total_tissue_lower_mm', 'total_tissue_upper_mm', 'enamel_lower_mm', 'enamel_upper_mm', 'reduction_lower_mm', 'reduction_upper_mm']
    (tl, tu, el, eu, rl, ru) = [float(row[k]) for k in keys]
    if not all((math.isfinite(v) and v >= 0 for v in [tl, tu, el, eu, rl, ru])):
        raise ValueError('Nonfinite/negative physical bounds')
    if tl > tu or el > eu or rl > ru or (el > tu):
        raise ValueError('Empty/invalid measurement interval')
    lo = max(0, tl - max(eu, ru))
    hi = max(0, tu - max(el, rl))

    def cls(b):
        return 'ABOVE_IN_ALL_WORLDS' if lo >= b else 'BELOW_IN_ALL_WORLDS' if hi < b else 'UNKNOWN_STRADDLES'
    return {'specimen_id': row['specimen_id'], 'fdi': int(row['fdi']), 'remaining_dentin_interval_mm': [lo, hi], 'threshold_0p5': cls(0.5), 'threshold_1': cls(1), 'resolution': 'PER_POINT', 'consumer_must_aggregate': 'One measured ray is not a whole horn/whole tooth bound', 'status': 'CONDITIONAL_ON_DECLARED_PHYSICAL_MEASUREMENT_CONTRACT', 'external_provenance_independently_verified': False, 'biological_risk': 'UNKNOWN', 'source_locator': row['source_locator']}

def self_check():
    base = {'specimen_id': 'TEST_ONLY', 'fdi': 36, 'measurement_kind': 'INDEPENDENT_PHYSICAL', 'unit': 'mm', 'source_locator': 'our_own_fixture:test-only', 'ray_frame_locator': 'our_own_fixture:ray', 'horn_and_DEJ_verified': 'true', 'total_tissue_lower_mm': 3, 'total_tissue_upper_mm': 3.2, 'enamel_lower_mm': 1.2, 'enamel_upper_mm': 1.4, 'reduction_lower_mm': 2, 'reduction_upper_mm': 2.1}
    a = local_dentin_interval(base)
    (lo, hi) = a['remaining_dentin_interval_mm']
    corners = []
    for t in [3, 3.2]:
        for e in [1.2, 1.4]:
            for r in [2, 2.1]:
                corners.append(max(0, t - max(e, r)))
    err = max(abs(lo - min(corners)), abs(hi - max(corners)))
    tests = [{'check': 'physical_interval_vs_independent_corners', 'valid_pass': err < 1e-07, 'error_mm': err, 'fault_rejected': abs(hi + 1 - max(corners)) > 1e-07, 'fixture_kind': 'our_own_fixture'}]
    enamel = base | {'reduction_lower_mm': 0.5, 'reduction_upper_mm': 0.5}
    val = local_dentin_interval(enamel)['remaining_dentin_interval_mm']
    expected = [1.6, 2.0]
    err = max((abs(val[i] - expected[i]) for i in [0, 1]))
    tests.append({'check': 'cut_before_DEJ', 'valid_pass': err < 1e-07, 'error_mm': err, 'fault_rejected': abs(3 - 0.5 - val[0]) > 1e-07, 'fixture_kind': 'our_own_fixture'})
    for (name, patch) in [('synthetic', {'measurement_kind': 'SYNTHETIC'}), ('unit', {'unit': 'm'}), ('DEJ', {'horn_and_DEJ_verified': 'false'}), ('locator', {'source_locator': ''}), ('NaN', {'enamel_lower_mm': 'nan'}), ('empty_interval', {'enamel_lower_mm': 4}), ('negative_reduction', {'reduction_lower_mm': -1})]:
        rejected = False
        try:
            local_dentin_interval(base | patch)
        except ValueError:
            rejected = True
        tests.append({'check': 'input_' + name, 'valid_pass': True, 'fault_rejected': rejected, 'fixture_kind': 'our_own_fixture'})
    spec = importlib.util.spec_from_file_location('cal', ROOT / 'inputs/calibrate_clearance.py')
    cal = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cal)
    sg = {'independent_block_id': 'TEST_ONLY', 'guide': 'fully_guided', 'planned_gap_mm': '4', 'measured_min_tool_gap_mm': '3', 'metrology_loss_upper_mm': '.1', 'source_locator': 'our_own_fixture:test-only', 'measurement_kind': 'INDEPENDENT_PHYSICAL', 'population_transfer_verified': 'true'}
    for (name, records) in [('synthetic', [sg | {'measurement_kind': 'SYNTHETIC'}]), ('duplicate', [sg, sg]), ('nonfinite', [sg | {'measured_min_tool_gap_mm': 'nan'}]), ('missing_locator', [sg | {'source_locator': ''}]), ('missing_transfer', [sg | {'population_transfer_verified': 'false'}])]:
        rejected = False
        try:
            cal.calibrate(records, x=0.05)
        except ValueError:
            rejected = True
        tests.append({'check': 'signed_loss_' + name, 'valid_pass': True, 'fault_rejected': rejected, 'fixture_kind': 'our_own_fixture'})
    if not all((t['valid_pass'] and t['fault_rejected'] for t in tests)):
        raise AssertionError('Measurement guard failed')
    (ROOT / 'raw/MEASUREMENT_CONTROLS.json').write_text(json.dumps(tests, indent=2) + '\n')
    print('Measured-input guard checks PASS; 0 physical observations supplied.')

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dentin-csv', type=Path)
    ap.add_argument('--output', type=Path)
    ap.add_argument('--self-check', action='store_true')
    a = ap.parse_args()
    if a.self_check:
        self_check()
        return
    rows = list(csv.DictReader(a.dentin_csv.open())) if a.dentin_csv else []
    seen = set()
    out = []
    for row in rows:
        key = (row['specimen_id'], row['fdi'], row['ray_frame_locator'])
        if key in seen:
            raise ValueError('Duplicate ray observation')
        seen.add(key)
        out.append(local_dentin_interval(row))
    obj = {'status': 'UNKNOWN_NO_PHYSICAL_MEASUREMENTS' if not rows else 'CONDITIONAL_LOCAL_MEASUREMENT_BOUNDS', 'physical_rows_n': len(rows), 'results': out, 'external_input_sha256': hashlib.sha256(a.dentin_csv.read_bytes()).hexdigest() if a.dentin_csv else None}
    text = json.dumps(obj, indent=2, allow_nan=False) + '\n'
    if a.output:
        a.output.write_text(text)
    else:
        print(text, end='')
if __name__ == '__main__':
    main()
