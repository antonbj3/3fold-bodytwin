"""R4 exact counterexample to a mean-preload sufficient state."""
import json, time
from fractions import Fraction as Q
import numpy as np
from assembly import beam_matrix, rational_solve, state, independent_state
from demo import ROOT, sha, write, now

def certificate(L, p):
    K = beam_matrix(L, exact=True)
    A = [r[:] for r in K]
    b = [Q(0)] * 8
    f = [Q(0)] * 8
    f[6] = Q(1)
    e = [Q(v) * (Q(1, 20000) + Q(1, 200000)) for v in p]
    for j in range(4):
        A[2 * j][2 * j] += 220000
        b[2 * j] = -20000 * e[j]
    x = rational_solve(A, b)
    v = rational_solve(A, f)
    assert all((x[2 * j] <= 0 and x[2 * j] + e[j] >= 0 for j in range(4)))
    roots = [(-x[2 * j] / v[2 * j], j) for j in range(4) if v[2 * j] > 0]
    roots += [(-(x[2 * j] + e[j]) / v[2 * j], j) for j in range(4) if v[2 * j] < 0]
    (onset, j) = min(roots)
    return {'mean_preload_fraction_N': str(sum((Q(z) for z in p)) / 4), 'onset_fraction_N': str(onset), 'onset_N': float(onset), 'support_index': j, 'initial_contact_fraction_N': [str(-200000 * x[2 * j]) for j in range(4)], 'exact_equilibrium_residual': str(max((abs(sum((A[i][k] * x[k] for k in range(8))) - b[i]) for i in range(8)))), 'affine_valid_interval_N': ['0', str(onset)], 'physical_error_enclosure': 'MISSING; this certificate covers the ideal quantized scenario only.'}

def main():
    t = time.perf_counter()
    p = json.loads((ROOT / 'PREREG_R4.json').read_text())
    geo = json.loads((ROOT / 'raw/PATIENT_MODEL.json').read_text())
    L = geo['segment_lengths_quantized_mm']
    fp = ROOT / 'FROZEN_PREDICTIONS_R4.json'
    if not fp.exists():
        write(fp, {'frozen_at': now(), 'prereg_sha256': sha(ROOT / 'PREREG_R4.json'), 'code_sha256': sha(__file__), 'states_N': p['two_states_N'], 'prediction': 'Identical mean and closed gaps; different guarded first event, physical source transport UNKNOWN.'})
        (ROOT / 'FROZEN_PREDICTIONS_R4.sha256').write_text(sha(fp) + '\n')
    else:
        assert sha(fp) == (ROOT / 'FROZEN_PREDICTIONS_R4.sha256').read_text().strip()
        assert sha(__file__) == json.loads(fp.read_text())['code_sha256']
    rows = []
    K = beam_matrix(L)
    for P in p['two_states_N']:
        a = state(K, np.zeros(4), np.array(P), np.zeros(4))
        b = independent_state(K, np.zeros(4), np.array(P), np.zeros(4))
        c = certificate(L, P)
        rows.append({'preload_vector_N': P, 'mean_preload_N': float(np.mean(P)), 'assembled_gaps_mm': a['gap'].tolist(), 'exact_certificate': c, 'independent_max_force_error_N': float(max(np.max(abs(a['contact'] - b['contact'])), np.max(abs(a['bolt'] - b['bolt']))))})
    (A, B) = rows
    dg = float(np.max(abs(np.array(A['assembled_gaps_mm']) - np.array(B['assembled_gaps_mm']))))
    dm = abs(A['mean_preload_N'] - B['mean_preload_N'])
    diff = abs(Q(A['exact_certificate']['onset_fraction_N']) - Q(B['exact_certificate']['onset_fraction_N']))
    passed = dg == dm == 0 and diff > 1 and (max((r['independent_max_force_error_N'] for r in rows)) < 1e-05)
    out = {'round': 'R4', 'claim_type': 'capability', 'external_referent': p['external_referent'], 'rows': rows, 'sufficiency': {'identity_error_mean_N': dm, 'identity_error_assembled_gap_mm': dg, 'downstream_difference_N': float(diff), 'downstream_difference_fraction_N': str(diff), 'outcome': 'MEAN_PRELOAD_SUMMARY_REFUTED' if passed else 'FAIL_GATE', 'minimum_extension_fixed_query': 'Per-support seat forces and their directional load response/validity guards; new arbitrary loads require the complete support transfer port.'}, 'physical_patient_validity': 'UNKNOWN; published population mean is not an observed4-support preload vector.', 'wall_s': time.perf_counter() - t, 'resolution_level': 'PHENOMENOLOGICAL', 'time_scale': 'HANDOVER'}
    write(ROOT / 'raw/RESULTS_R4.json', out)
    assert passed
    print(json.dumps({'round': 'R4', 'mean_identity_error_N': dm, 'gap_identity_error_mm': dg, 'onset_difference_N': float(diff)}))
if __name__ == '__main__':
    main()
