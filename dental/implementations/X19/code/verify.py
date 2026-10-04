"""Verify frozen scientific contracts and execute adverse-value controls."""
from dental_release.paths import expand as _release_expand
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[k] = '4'
from pathlib import Path
import hashlib, json, sys, zipfile, time
import numpy as np
from series import judge, fit, stiffness
from model import response
R = Path(__file__).resolve().parents[1]

def read(n):
    return json.load(open(R / n))

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as fh:
        for b in iter(lambda : fh.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()

def run():
    t = time.perf_counter()
    checks = {}
    for i in [1, 2, 3]:
        freeze = read(f'FROZEN_PREDICTIONS_R{i}.json')
        checks[f'R{i}_frozen_prediction_hash'] = sha(R / freeze['prediction_file']) == freeze['prediction_sha256']
        checks[f'R{i}_prereg_hash'] = sha(R / f'PREREG_R{i}.json') == freeze['prereg_sha256']
    checks['r1_input_hashes'] = all((sha(p) == s for (p, s) in [(R / 'inputs/ARCH_GEOMETRY.json', read('raw/PREDICTIONS_R1.json')['sha256_inputs']['geometry']), (Path(_release_expand('@DENTAL_INPUT_ROOT@/workspace/notes/aligner_intervals/ALIGNER_INTERVALS.json')), read('raw/PREDICTIONS_R1.json')['sha256_inputs']['intervals'])]))
    checks['r2_calibration_hash'] = sha(R / 'inputs/R2_CALIBRATION_ENDPOINTS.json') == read('raw/PREDICTIONS_R2.json')['calibration_sha256']
    for i in [1, 2]:
        pp = read(f'raw/PREDICTIONS_R{i}.json')
        checks[f'R{i}_frozen_code_hashes'] = all((sha(R / 'code' / n) == h for (n, h) in pp['code_sha256'].items()))
    zips = {}
    source_ok = True
    for s in read('inputs/ARCH_GEOMETRY.json')['source_manifest']:
        if 'zip' in s:
            if s['zip'] not in zips:
                zips[s['zip']] = zipfile.ZipFile(s['zip'])
            b = zips[s['zip']].read(s['member'])
            source_ok = source_ok and hashlib.sha256(b).hexdigest() == s['sha256']
        else:
            source_ok = source_ok and sha(Path(s['path'])) == s['sha256']
    for z in zips.values():
        z.close()
    checks['all_selected_source_hashes'] = source_ok
    a = read('inputs/ARCH_GEOMETRY.json')['arches'][0]
    k = int(next(iter(a['teeth'])))
    good = response(a, k, 0.2, 2189, 0.53746)
    ctrl = response(a, k, 0.2, 2189, 0.53746, True)
    bad = response(a, k, 0.2, 2189, 0.53746 * 1.2, True)

    def diff(x, y):
        return max((np.max(np.abs(np.array(x['wrenches'][q]) - y['wrenches'][q])) for q in x['wrenches']))
    checks['Hermite_correct_input_accepts'] = diff(good, ctrl) <= 1e-10
    checks['Hermite_wrong_h_rejects'] = diff(good, bad) > 1e-10
    F = sum((np.array(w)[:3] for w in good['wrenches'].values()))
    checks['force_balance_accepts'] = np.linalg.norm(F) <= 1e-09
    checks['force_plus1N_rejects'] = np.linalg.norm(F + np.array([1, 0, 0])) > 1e-09
    M = sum((np.array(w)[3:] + np.cross(a['teeth'][q]['centroid_mm'], np.array(w)[:3]) for (q, w) in good['wrenches'].items()))
    checks['torque_balance_accepts'] = np.linalg.norm(M) <= 1e-08
    checks['torque_plus1Nmm_rejects'] = np.linalg.norm(M + np.array([1, 0, 0])) > 1e-08
    rs = [r for r in read('raw/RESULTS_R2.json')['rows'] if r['primary']]
    p = [r['stiffness_N_per_mm'] for r in rs]
    v = [r['measured_K_N_per_mm'] for r in rs]
    checks['R2_real_primary_accepts'] = judge(p, v, 0.15)['gate'] == 'PASS'
    checks['R2_doubled_prediction_rejects'] = judge([2 * a for a in p], v, 0.15)['gate'] == 'FAIL'
    ep = read('inputs/R2_CALIBRATION_ENDPOINTS.json')['calibration']['labial']
    f = fit(ep)
    c = fit(ep, True)
    checks['bounded_fit_control_accepts'] = abs(stiffness(f, 0.625) - stiffness(c, 0.625)) / stiffness(f, 0.625) <= 1e-10
    checks['bounded_fit_corrupt_plus20percent_rejects'] = abs(stiffness(f, 0.625) - 1.2 * stiffness(c, 0.625)) / stiffness(f, 0.625) > 1e-10
    from calibrate_lab import calibrate
    checks['absent_local_bench_refuses'] = calibrate(read('inputs/MEASUREMENTS_TEMPLATE.json'))['status'] == 'NEEDS_MATCHED_MEASUREMENT'
    r = read('raw/RELAXATION_CASES.json')
    checks['relaxation_cases_not_pooled'] = r['pooled_curve'] is None and len({a['load_case'] for a in r['load_cases']}) == 6
    checks['unloaded_aging_not_called_relaxation'] = all((not a['is_relaxation_under_load'] for a in r['load_cases'] if 'no_applied_load' in a['load_case']))
    x = read('results.json')
    checks['all_physical_force_classifications_unknown'] = x['biological_reference']['patient_classification'] == 'UNKNOWN' and x['negative_result'] is True
    files = [p for p in R.rglob('*') if p.is_file()]
    size = sum((p.stat().st_size for p in files))
    checks['disk_cap'] = size <= 3000000000
    checks['large_arrays_on_data_disk'] = all((p.stat().st_size <= 50000000 for p in files))
    checks = {k: bool(v) for (k, v) in checks.items()}
    out = {'status': 'PASS' if all(checks.values()) else 'FAIL', 'checks': checks, 'wall_seconds': time.perf_counter() - t, 'lane_bytes': size, 'fault_injection_is_numerical_contract_not_external_facit': True}
    (R / 'raw/VERIFICATION.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out, indent=2))
    assert all(checks.values())
    state = read('CURRENT_WORK_STATE.json')
    state.update(phase='THREE_ROUNDS_ADJUDICATED_DEMO_VERIFIED_NEXT_BENCH_READY', latest_software_verification='PASS', verification_checks=len(checks), data_bytes=size)
    (R / 'CURRENT_WORK_STATE.json').write_text(json.dumps(state, indent=2) + '\n')
if __name__ == '__main__':
    run()
