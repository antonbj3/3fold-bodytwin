from patient_geometry import P, sha, write
import json, numpy as np
from scipy.optimize import linprog
checks = []

def check(k, v):
    checks.append(dict(name=k, passed=bool(v)))
    if not v:
        raise AssertionError(k)
g = json.load(open(P / 'rounds/R1/results.json'))
c = json.load(open(P / 'rounds/R2/results.json'))
m = json.load(open(P / 'rounds/R3/results.json'))
x = json.load(open(P / 'rounds/R4/results.json'))
check('27 actual source crowns', g['source_crowns'] == 27)
check('physical Demo1 force remains unknown', g['physical_gate'].startswith('UNKNOWN') and m['physical_K16'].startswith('STILL_MISSING'))
check('all geometry witness controls pass', g['controls']['all_pass'])
check('external calibration failure preserved', c['n_fail'] == 27 and c['n_pass'] == 21 and c['outcome'].startswith('FAIL'))
check('published gains not transferred', c['patient_gain_transfer'].startswith('REJECTED'))
check('external wrong patient not admitted', m['published_force_transfer_rejected'])
check('independent contact LP and live contract faults', m['all_control_gates_pass'] and len(m['faults']) == 22)
check('external totals and fixed100N benchmark', all((v['balance_error_N'] == 0 for v in m['external_rows'])) and m['fixed100N_display_exceedances'] == 5)
check('exact contact equilibrium identities and corrupt compliance', x['all_pass'])
errors = []
for w in g['pairs']:
    U = np.array(w['upper_triangle_cbct_mm'])
    L = np.array(w['lower_triangle_cbct_mm'])
    A = np.zeros((4, 6))
    A[0, :3] = 1
    A[1, 3:] = 1
    A[2, :3] = U[:, 0]
    A[2, 3:] = -L[:, 0]
    A[3, :3] = U[:, 1]
    A[3, 3:] = -L[:, 1]
    ans = linprog(np.r_[U[:, 2], -L[:, 2]], A_eq=A, b_eq=[1, 1, 0, 0], bounds=(0, None), method='highs', options={'threads': 4})
    errors.append(abs(ans.fun - w['minimum_axial_gap_mm']))
check('27 independently replayed winning LPs', max(errors) <= 1e-07)
frozen = json.load(open(P / 'FROZEN_PREDICTIONS.json'))
check('first frozen result preserved', frozen['result_sha256'] == sha(P / 'rounds/R1/FIRST_RESULTS.json'))
check('frozen geometric prediction unchanged', frozen['predictions']['target36_first_touch_gap_mm'] == g['target36_first_touch_gap_mm'])
write(P / 'raw/VERIFICATION.json', dict(status='PASS', review_state='PENDING_INDEPENDENT_REVIEW', scope='Producer replay and adverse controls; not independent scientific review', checks=checks, maximum_independent_LP_error_mm=max(errors)))
print('PASS', len(checks), 'checks; LP max mm', max(errors))
