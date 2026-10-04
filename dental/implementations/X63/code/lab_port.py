"""Freeze measurable predictions; ingest a real CSV only when one exists."""
import argparse, csv, json, hashlib, sys
import numpy as np
from pathlib import Path
from assembly import beam_matrix, state
from demo import ROOT, write, sha, now
REQUIRED = ['prediction_key', 'contact_force_N', 'gap_mm', 'force_unit', 'gap_unit', 'source_locator', 'calibration_locator', 'fixture_geometry_sha256', 'condition_verified', 'independent_measurement', 'prediction_sha256']

def freeze():
    p = ROOT / 'FROZEN_LAB_PREDICTIONS.json'
    if p.exists():
        assert sha(p) == (ROOT / 'FROZEN_LAB_PREDICTIONS.sha256').read_text().strip()
        return json.loads(p.read_text())
    geo = json.loads((ROOT / 'raw/PATIENT_MODEL.json').read_text())
    L = geo['segment_lengths_quantized_mm']
    rows = []
    for h in [4, 5]:
        for g in [0, 0.02]:
            for (use, P) in [('first_force_scenario', 329.9), ('tenth_force_scenario', 253.7)]:
                for load in [0, 100]:
                    a = state(beam_matrix(L, h=h), np.array([g, -g, g, -g]), P, np.array([0, 0, 0, load]))
                    for support in range(4):
                        key = f'h{h}_g{g * 1000:.0f}_{use}_load{load}_support{support + 1}'
                        rows.append({'prediction_key': key, 'height_mm': h, 'misfit_vector_mm': [g, -g, g, -g], 'preload_target_each_bolt_N': P, 'lifting_load_on_support4_N': load, 'support': support + 1, 'contact_force_N': float(a['contact'][support]), 'bolt_tension_N': float(a['bolt'][support]), 'gap_mm': float(a['gap'][support]), 'resolution_level': 'PER_SURFACE_REGION', 'physical_status': 'CONDITIONAL_UNMEASURED'})
    d = {'frozen_at': now(), 'prereg_sha256': sha(ROOT / 'PREREG_R3.json'), 'patient_geometry_sha256': geo['member_sha256'], 'prediction_code_sha256': {'assembly.py': sha(ROOT / 'code/assembly.py'), 'lab_port.py': sha(ROOT / 'code/lab_port.py')}, 'measurement_status': 'NOT_MEASURED', 'rows': rows, 'comparison_tolerances': {'force': '10+0.10*abs(predicted_N)', 'gap_mm': 0.01}, 'claim_type': 'capability', 'time_scale': 'HANDOVER', 'regime': 'Force-controlled instrumented research fixture matching spatial geometry, springs, cross section and shims. Source torque alone does not meet force target.'}
    write(p, d)
    (ROOT / 'FROZEN_LAB_PREDICTIONS.sha256').write_text(sha(p) + '\n')
    with (ROOT / 'MEASUREMENTS_TEMPLATE.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=REQUIRED)
        w.writeheader()
        for row in rows:
            w.writerow({'prediction_key': row['prediction_key'], 'force_unit': 'N', 'gap_unit': 'mm', 'fixture_geometry_sha256': geo['member_sha256'], 'prediction_sha256': sha(p)})
    return d

def compare(d, measurements):
    predictions = {r['prediction_key']: r for r in d['rows']}
    seen = set()
    out = []
    for row in measurements:
        key = row.get('prediction_key')
        reasons = []
        if key not in predictions:
            reasons.append('UNKNOWN_PREDICTION_KEY')
        if key in seen:
            reasons.append('DUPLICATE_SUPPORT_RECORD')
        seen.add(key)
        for k in REQUIRED:
            if not str(row.get(k, '')).strip():
                reasons.append('MISSING_' + k)
        for (k, v) in [('force_unit', 'N'), ('gap_unit', 'mm'), ('condition_verified', 'true'), ('independent_measurement', 'true')]:
            if row.get(k) != v:
                reasons.append('INVALID_' + k)
        if row.get('prediction_sha256') != sha(ROOT / 'FROZEN_LAB_PREDICTIONS.json'):
            reasons.append('FROZEN_HASH_MISMATCH')
        if row.get('fixture_geometry_sha256') != d['patient_geometry_sha256']:
            reasons.append('GEOMETRY_MISMATCH')
        try:
            force = float(row['contact_force_N'])
            gap = float(row['gap_mm'])
            assert np.isfinite(force) and np.isfinite(gap)
        except (KeyError, ValueError, AssertionError):
            reasons.append('MISSING_OR_INVALID_MEASURED_VALUES')
        if reasons:
            out.append({'prediction_key': key, 'status': 'UNKNOWN', 'reasons': reasons})
            continue
        pred = predictions[key]
        fe = abs(force - pred['contact_force_N'])
        ge = abs(gap - pred['gap_mm'])
        out.append({'prediction_key': key, 'status': 'REJECT_PREDICTION' if fe > 10 + 0.1 * abs(pred['contact_force_N']) or ge > 0.01 else 'WITHIN_FROZEN_TOLERANCE', 'force_residual_N': fe, 'gap_residual_mm': ge, 'source_locator': row['source_locator']})
    missing = sorted(set(predictions) - seen)
    return {'measurement_status': 'NOT_MEASURED' if not measurements else 'SUPPLIED', 'rows': out, 'missing_keys': missing, 'campaign_status': 'FAIL' if any((r['status'] == 'REJECT_PREDICTION' for r in out)) else 'UNKNOWN' if missing or any((r['status'] == 'UNKNOWN' for r in out)) else 'WITHIN_SCOPED_TOLERANCE', 'scope': 'CSV metadata is an asserted measurement port, not automatically verified physical provenance; independent review remains required.'}

def selfcheck(d):
    r = d['rows'][0]
    base = {k: 'control' for k in REQUIRED}
    base.update({'prediction_key': r['prediction_key'], 'contact_force_N': str(r['contact_force_N']), 'gap_mm': str(r['gap_mm']), 'force_unit': 'N', 'gap_unit': 'mm', 'condition_verified': 'true', 'independent_measurement': 'true', 'source_locator': 'our_own_fixture_validator_probe', 'calibration_locator': 'our_own_fixture_validator_probe', 'fixture_geometry_sha256': d['patient_geometry_sha256'], 'prediction_sha256': sha(ROOT / 'FROZEN_LAB_PREDICTIONS.json')})
    bad = dict(base)
    bad['contact_force_N'] = str(r['contact_force_N'] + 1000)
    unit = dict(base)
    unit['force_unit'] = 'Ncm'
    missing = dict(base)
    missing['source_locator'] = ''
    result = {'synthetic_validator_checks_only': True, 'not_measured_data': True, 'reject_1000N': compare(d, [bad])['rows'][0]['status'] == 'REJECT_PREDICTION', 'reject_wrong_unit': compare(d, [unit])['rows'][0]['status'] == 'UNKNOWN', 'reject_missing_locator': compare(d, [missing])['rows'][0]['status'] == 'UNKNOWN', 'missing_real_data_abstains': compare(d, [])['campaign_status'] == 'UNKNOWN'}
    write(ROOT / 'raw/LAB_VALIDATOR_INJECTED_ERRORS.json', result)
    assert all((v for (k, v) in result.items() if k.startswith(('reject', 'missing'))))
    write(ROOT / 'raw/LAB_COMPARISON.json', compare(d, []))
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--measurements')
    a = p.parse_args()
    d = freeze()
    if a.measurements:
        with open(a.measurements) as f:
            res = compare(d, list(csv.DictReader(f)))
        write(ROOT / 'raw/SUPPLIED_LAB_COMPARISON.json', res)
        print(res['campaign_status'])
        sys.exit(1 if res['campaign_status'] == 'FAIL' else 0)
    else:
        selfcheck(d)
        print('64 frozen support predictions; actual measurements UNKNOWN')
