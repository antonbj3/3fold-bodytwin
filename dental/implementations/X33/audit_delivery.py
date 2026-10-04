"""Read-only integrity and deliberately corrupted decision audit."""
from pathlib import Path
import json, hashlib, numpy as np
from scipy.optimize import linprog
ROOT = Path(__file__).resolve().parent

def read(p):
    return json.loads((ROOT / p).read_text())

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def agree(actual, expected, tol):
    return abs(actual - expected) <= tol

def main():
    checks = {}
    for j in [1, 2, 3]:
        checks[f'PREREG_R{j}_hash'] = sha(ROOT / f'PREREG_R{j}.json') == (ROOT / f'PREREG_R{j}.sha256').read_text().split()[0]
    r1 = read('raw/R1_RESULTS.json')
    r2 = read('raw/R2_RESULTS.json')
    r3 = read('raw/R3_RESULTS.json')
    checks['pulse_matrix_hash'] = sha(r2['matrix_path']) == r2['matrix_sha256']
    checks['X26_abstains'] = r1['decision']['X26']['class'] == 2 and r1['decision']['X26']['physical_status'] == 'UNKNOWN'
    checks['R1_failed_gates_preserved'] = not r1['refinement_gate_pass'] and (not r1['external']['frozen_proxy_MAE_gate_pass'])
    checks['energy_and_scalar_replays'] = all((r['energy_gate_pass'] for r in r1['refinement'])) and all(r2['gates'].values())
    checks['R3_failed_consistency_preserved'] = not r3['resolution_consistency_gate_pass'] and r3['coarse_acceptances_refuted_by_finer_feasible_history'] == 3
    M = np.load(r2['matrix_path'], mmap_mode='r').reshape(-1, 42)
    row = r3['rows'][-1]
    B = row['bins']
    rel = row['relative_measurement_error']
    G = np.eye(42)
    fit = linprog(-np.asarray(M[row['worst_row']]), A_ub=np.vstack([G, -G, np.ones((1, 42)), -np.ones((1, 42))]), b_ub=np.r_[np.full(B, 2 / B * (1 + rel)), np.full(B, -2 / B * (1 - rel)), 2.02, -1.98], bounds=[(0, 0.1)] * 42, method='highs', options={'threads': 4})
    checks['real_LP_reference'] = fit.success and agree(row['worst_peak_pulp_C'], -fit.fun, 1e-07)
    checks['same_LP_validator_rejects_plus1'] = not agree(row['worst_peak_pulp_C'] + 1, -fit.fun, 1e-07)
    checks['same_replay_validator_rejects_reversed_pulses'] = not agree(r2['reverse_adverse_history_peak_C'], r2['maximum_possible_maximum_C'], 0.02)

    def below_gate(peak):
        return peak < 5.5
    checks['same_threshold_gate_rejects_coarse_label_on_fine'] = below_gate(r3['cross_resolution_uniform_history'][0]['maximum_C']) and (not below_gate(r3['cross_resolution_uniform_history'][-1]['maximum_C']))
    original = read('results.json')
    feedback = read('GRAPH_FEEDBACK.json')
    checks['feedback_hash_matches'] = feedback['sha256'] == sha(ROOT / 'results.json')
    checks['no_physical_certificate'] = original['physical_decision']['physical_certification'] == 'UNKNOWN' and original['physical_decision']['certifying_resolution_mm'] is None
    checks['isolated_replay'] = read('LATEST_REPLAY.json')['status'] == 'PASS_REPLAY'
    checks['mutations_original'] = all(r1['mutation_tests'].values()) and all(r2['mutations'].values()) and all(r3['mutations'].values())
    out = {'status': 'PASS' if all(checks.values()) else 'FAIL', 'checks': checks, 'scientific_admission': False, 'review_state': 'PENDING_INDEPENDENT_REVIEW'}
    (ROOT / 'DELIVERY_AUDIT.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out), flush=True)
    assert all(checks.values())
if __name__ == '__main__':
    main()
