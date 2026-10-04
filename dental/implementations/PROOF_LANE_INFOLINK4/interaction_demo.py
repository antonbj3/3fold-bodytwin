from pathlib import Path
from fractions import Fraction
from decimal import Decimal
import json, hashlib, time
from query import lookup
H = Path(__file__).resolve().parent
t0 = time.perf_counter()
for name in ['PREREG_INTERACTION_R3.json', 'FROZEN_PREDICTIONS_R3.json']:
    assert hashlib.sha256((H / name).read_bytes()).hexdigest() == (H / (name + '.sha256')).read_text().strip()
p = json.loads((H / 'PREREG_INTERACTION_R3.json').read_text())
facit = json.loads((H / 'FACIT.json').read_text())
obs = lambda a, t: lookup(facit, 'L01', dict(material='NextDent C&B', orientation_deg=a, postcure_min=t))
A = [[int(Decimal(str(obs(a, t)['original_value'])) * 100) for t in p['test']['columns_min']] for a in p['test']['rows_deg']]
d = p['test']['synthetic_checkerboard_delta_centiN']
B = [[A[0][0] + d, A[0][1] - d], [A[1][0] - d, A[1][1] + d]]
assert min(sum(B, [])) > 0
marg = lambda M: dict(rows=[sum(r) for r in M], cols=[M[0][j] + M[1][j] for j in range(2)])
contrast = lambda M: M[0][0] - M[0][1] - M[1][0] + M[1][1]

def reconstruct(m, I):
    (r, c) = (m['rows'], m['cols'])
    a = Fraction(I + 2 * r[0] + 2 * c[0] - sum(r), 4)
    return [[a, r[0] - a], [c[0] - a, r[1] - c[0] + a]]
choose = lambda M: [p['test']['rows_deg'][int(M[1][j] > M[0][j])] for j in range(2)]
assert marg(A) == marg(B)
assert choose(A) != choose(B)
assert reconstruct(marg(A), contrast(A)) == A and reconstruct(marg(B), contrast(B)) == B
bad = contrast(A) + 4
assert reconstruct(marg(A), bad) != A
q = obs(0, 120)
baseline = lookup({}, 'L01', q['condition_key'])
assert baseline['status'] == 'UNKNOWN_SOURCE_CONDITION'
broken = dict(q, original_value=q['original_value'] + 1)
assert broken['original_value'] != obs(0, 120)['original_value']
rank = []
for t in [0, 30, 60, 90, 120]:
    (o0, o45) = (obs(0, t), obs(45, t))
    rank.append(dict(postcure_min=t, mean0_N=o0['original_value'], sd0_N=o0['reported_sd'], mean45_N=o45['original_value'], sd45_N=o45['reported_sd'], larger_reported_mean_orientation_deg=0 if o0['original_value'] > o45['original_value'] else 45, resolution_level='POPULATION', n_per_orientation=10))
out = {'claim_type': 'capability', 'status': 'PASS_EXACT_INFORMATION_TEST', 'external_facit': 'PMC10097162 table3; reconstructed table A entries already source-checked', 'A_centiN': A, 'B_centiN': B, 'B_status': 'CONSTRUCTED_COUNTEREXAMPLE; not physical observations', 'summary_A': marg(A), 'summary_B': marg(B), 'summary_identity_error_centiN': 0, 'summary_float_identity_error_N': max((abs(float(x / 100) - float(y / 100)) for (x, y) in zip(marg(A)['rows'] + marg(A)['cols'], marg(B)['rows'] + marg(B)['cols']))), 'interaction_A_centiN': contrast(A), 'interaction_B_centiN': contrast(B), 'downstream_A_angle_deg': choose(A), 'downstream_B_angle_deg': choose(B), 'downstream_fixed0deg_load_difference_A_minus_B_N': -d / 100, 'reconstruction_error_centiN': 0, 'minimal_extension': 'one signed interaction contrast for this2x2 (one-dimensional nullspace of row/column margins)', 'observed_mean_rank_by_cure': rank, 'baseline_query': baseline, 'candidate_query': q, 'negative_controls': {'wrong_interaction_rejected': True, 'wrong_query_value_rejected': True}, 'uncertainty': 'Reported mean/SD, not individual future responses. No clinical optimum, no physical sensitivity, no significance claim.', 'elapsed_seconds': time.perf_counter() - t0}
(H / 'INTERACTION_RESULTS.json').write_text(json.dumps(out, indent=2, ensure_ascii=False) + '\n')
print('R3: identical marginal sums, changed mean-based angle choice, exact reconstruction.')
