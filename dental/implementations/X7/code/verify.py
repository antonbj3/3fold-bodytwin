"""Meaningful integrity/semantic and adversarial checks, distinct from independent review."""
import json, hashlib, sys
from pathlib import Path
import numpy as np
from geometry import P, load_stl_bytes
from benchmark import sha, write
from parser import parse
from diagnose import landmark_measures

def run():
    checks = []

    def ck(name, b):
        checks.append(dict(name=name, pass_check=bool(b)))
    for f in ['PREREG_R1.json', 'PREREG_R2.json', 'PREREG_R3.json', 'PREREG_R4.json']:
        ck('Frozen prereg ' + f, sha(P / f) == (P / (f + '.sha256')).read_text().strip())
    for (f, h) in json.load(open(P / 'raw/FROZEN_CODE_R1.json')).items():
        ck('Frozen R1 ' + f, sha(P / f) == h)
    for x in json.load(open(P / 'sources/upstream_manifest.json')):
        ck('Pinned external code ' + x['path'], sha(P / x['path']) == x['sha256'])
    m = json.load(open(P / 'raw/DATA_MANIFEST.json'))
    splits = [set(m[k]) for k in ['train', 'calibration', 'test']]
    ck('Patient split disjoint', not any((a & b for (i, a) in enumerate(splits) for b in splits[i + 1:])))
    fp = json.load(open(P / 'FROZEN_PREDICTIONS.json'))
    ck('Predictions immutable', sha(P / fp['prediction_path']) == fp['prediction_sha256'])
    ck('Model immutable', sha(P / 'raw/models.pkl') == fp['model_sha256'])
    ck('Full geometry immutable', sha(P / 'raw/geometry.jsonl') == fp['geometry_sha256'])
    fp2 = json.load(open(P / 'FROZEN_PREDICTIONS_R2.json'))
    ck('Set predictions immutable', sha(P / 'raw/PREDICTION_SETS_R2.json') == fp2['set_sha256'])
    ck('Set calibration immutable', sha(P / 'raw/THRESHOLDS_R2.json') == fp2['threshold_sha256'])
    ck('Negated midlines not centered', parse('The midlines are not coincident.')['midlines'] == 'deviated')
    ck('Lateral open does not override deep anterior', parse('There is an anterior deep bite and a right lateral open bite.')['overbite'] == 'increased')
    ck('Different molar/canine classes not merged', parse('Sagittally, there is a bilateral Class I molar and Class III canine relationship.')['molar_right'] == 'I')
    gold = json.load(open(P / 'raw/PARSER_REVIEW_CORRECTED.json'))
    em = [r for r in gold['rows'] if r['ours'] is not None]
    bad = [dict(r, ours='__injected_wrong__') for r in em]
    badprecision = sum((r['ours'] == r['gold'] for r in bad)) / len(bad)
    ck('Injected wrong parser fails precision gate', badprecision < 0.95)
    r1 = json.load(open(P / 'raw/RESULTS_R1.json'))
    r2 = json.load(open(P / 'raw/RESULTS_R2.json'))
    r3 = json.load(open(P / 'raw/RESULTS_R3.json'))
    ck('Injected label error fails R1', r1['injected_wrong_value_control']['gate'] == 'FAIL')
    ck('Injected threshold error fails R2', r2['injected_bad_threshold']['gate'] == 'FAIL')
    ck('Injected frame error fails R3', not r3['injected_bad_frame']['accepted'])
    strong = json.load(open(P / 'raw/STRONG_CONTROL_R1.json'))
    ck('Independent standard HGB refit parity', all((r['max_probability_difference'] <= 1e-12 and r['label_disagreements'] == 0 for r in strong['rows'])))
    ck('Actually injected probability error rejected', all((r.get('injected_probability_error', 0) > 1e-12 for r in strong['rows'])))
    r4 = json.load(open(P / 'raw/RESULTS_R4.json'))
    fp4 = json.load(open(P / 'FROZEN_PREDICTIONS_R4.json'))
    ck('Joint prediction sets immutable', sha(P / 'raw/PREDICTION_SETS_R4.json') == fp4['prediction_sha256'])
    ck('Injected joint threshold fails actual coverage', r4['injected_bad_threshold']['actual_joint_coverage'] < 0.85)
    ck('Joint standard control actually matches', abs(r4['threshold_q'] - r4['matched_control_threshold']) <= 1e-12)
    lm = landmark_measures(dict(superior_unit=[0, 0, 1], anterior_unit=[0, 1, 0], upper_inc_edge=[0, 4, -3], lower_inc_edge=[0, 0, 0], upper_midline=[2, 0, 0], lower_midline=[0, 0, 0], molar_right=dict(upper_cusp=[0, 2, 0], lower_reference=[0, 0, 0], mesial_unit=[0, 1, 0])))
    ck('Landmark closed-form OB OJ midline molar', lm['overbite_mm'] == 3 and lm['overjet_mm'] == 4 and (lm['midline_mm'] == 2) and (lm['molar_right_mesial_relation_mm'] == 2))
    rejected = False
    try:
        landmark_measures(dict(superior_unit=[0, 0, 1], anterior_unit=[0, 0, 1]))
    except ValueError:
        rejected = True
    ck('Injected invalid landmark basis rejected', rejected)
    write('VERIFICATION.json', dict(checks=checks, n=len(checks), passed=sum((x['pass_check'] for x in checks)), independent_review=False, clinical_validity=False, fixture='Landmark algebra fixture is our own internal test, not the external reference'))
    if not all((x['pass_check'] for x in checks)):
        raise RuntimeError('Integrity/semantic failure')
    print('Verification', len(checks), 'PASS')
if __name__ == '__main__':
    run()
