"""Analytic binomial and independently frozen transport/result checks."""
from pathlib import Path
import json, math
P = Path(__file__).resolve().parent

def check_result(r):
    s = r['statistical_answers']
    e = r['export']
    return e['meshes'] == 11 and e['conversion_paths'] == 55 and (e['transport_passes'] == 55) and (e['maximum_correspondence_mm'] <= 1e-05) and (e['region_probe']['mismatches'] == 0) and (e['region_probe']['facets'] == 6704) and (e['screening']['rejected'] == 4) and (s['plan29']['n_per_condition'] == math.floor(math.log(0.05) / math.log(0.9)) + 1) and (s['plan42']['n_per_condition'] == math.floor(math.log(0.05 / 4) / math.log(0.9)) + 1) and (abs(s['actual100N']['population_upper_one_sided'] - (1 - 0.05 ** (1 / 3))) < 1e-12) and (abs(s['plan29']['upper_if_acceptance_boundary'] - (1 - 0.05 ** (1 / 29))) < 1e-12) and (abs(s['plan42']['upper_if_acceptance_boundary'] - (1 - (0.05 / 4) ** (1 / 42))) < 1e-12) and (s['all_observation_screening'] == {'total': 36, 'complete': 27, 'early_censored': 9, 'rejected': 0}) and (s['binary']['n_per_group'] == 69) and (s['continuous']['n_per_group'] == 26)

def main():
    r = json.loads((P / 'results.json').read_text())
    assert check_result(r)
    bad = json.loads(json.dumps(r))
    bad['statistical_answers']['actual100N']['population_upper_one_sided'] = 0.05
    assert not check_result(bad)
    bad = json.loads(json.dumps(r))
    bad['export']['maximum_correspondence_mm'] = 0.01
    assert not check_result(bad)
    log = (P / 'raw/tests_latest.log').read_text()
    assert 'Ran 15 tests' in log and '\nOK\n' in log
    assert not (P / 'raw/runs').is_symlink()
    (P / 'REVIEW_CONTROL_RESULT.json').write_text(json.dumps({'status': 'PASS', 'binomial_formula_checks': True, 'injected_risk_rejected': True, 'injected_geometry_error_rejected': True, 'fresh_test_log': 'raw/tests_latest.log', 'tests': 15, 'mutable_output_symlink': False, 'scope': 'Digital transport and conditional fixed-plan statistics; physical performance UNKNOWN'}, indent=2) + '\n')
    print('Review control PASS: analytic binomial checks and risk/geometry injections rejected.')
if __name__ == '__main__':
    main()
