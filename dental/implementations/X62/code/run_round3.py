import json, time
import numpy as np
from common import ROOT, DATA, save, utc, sha, state

def reciprocal_gate(S):
    return bool(np.max(np.abs(S - S.T)) <= 1e-10 and np.linalg.eigvalsh((S + S.T) / 2).min() >= -1e-08)

def run():
    start = time.perf_counter()
    port = np.load(DATA / 'patient_port_R2.npz')
    S = port['S']
    state('R3_RUNNING', 'R3 consumer-specific reciprocity experiment frozen', '3 measured-response queries versus full12-mode identification')
    rows = json.loads((ROOT / 'raw/PATIENT_GEOMETRY.json').read_text())['rows']
    i = 2
    directions = [np.array(rows[i][k]) for k in ['mesial_axis', 'long_axis', 'labial_axis']]
    W = np.zeros((3, 18))
    for (j, d) in enumerate(directions):
        W[j, 3 * i:3 * i + 3] = d
    responses = S @ W.T
    rng = np.random.default_rng(62)
    queries = rng.uniform(-1 / 1024, 1 / 1024, (12, 18))
    direct = np.array([W @ S @ u for u in queries])
    pulled = np.array([responses.T @ u for u in queries])
    maxerr = float(np.max(np.abs(direct - pulled)))
    sketch = dict(created_utc=utc(), claim_type='capability', consumer='FDI31 initial force3 local components', imposed_patterns=(W.reshape(3, 6, 3) / 1024).tolist(), predicted_reactions=(responses.T.reshape(3, 6, 3) / 1024).tolist(), note='measure reactions to each prescribed small displacement; divide by amplitude before reuse', prereg_sha256=sha(ROOT / 'PREREG_R3.json'), reference_port_sha256=sha(DATA / 'patient_port_R2.npz'), physical_measurement='NOT_RUN', validity='same geometry, supports and reciprocal linear law; no hysteresis or interface slip')
    p = ROOT / 'FROZEN_QUERY_SPECIFIC_PREDICTIONS.json'
    if p.exists():
        old = json.loads(p.read_text())
        assert old['predicted_reactions'] == sketch['predicted_reactions']
        sketch = old
    else:
        save(p.name, sketch)
    A = np.diag([1.0, 1.0])
    B = np.diag([1.0, 2.0])
    w = np.array([1.0, 0.0])
    v = np.array([0.0, 1.0])
    assert np.array_equal(A @ w, B @ w)
    exact = dict(summary_name='response to one observed deformation', summary_A=(A @ w).tolist(), summary_B=(B @ w).tolist(), identity_error=0.0, unobserved_response_A=(A @ v).tolist(), unobserved_response_B=(B @ v).tolist(), downstream_difference=1.0, minimal_extension='a response query spanning the additional consumer direction', witness_kind='our_own_fixture')
    bad = S.copy()
    bad[0, 5] += 10
    blocked = not reciprocal_gate(bad)
    r = dict(claim_type='capability', resolution='PER_TOOTH', timescale='SIMULTANEOUS', output_max_error_N=maxerr, gate_pass=maxerr <= 1e-09 and reciprocal_gate(S) and blocked, symmetry_residual=float(np.max(np.abs(S - S.T))), minimum_eigenvalue=float(np.linalg.eigvalsh(S).min()), required_response_queries=3, full_operator_queries=12, reduction_scope='counts of response experiments for restricted linear3-component consumer; not measured total cost advantage', numerical_cases=12, skew_operator_injection_rejected=blocked, sufficiency=exact, external_referent=json.loads((ROOT / 'PREREG_R3.json').read_text())['external_referent'], physical_gate='NOT_RUN; R1 failed physical closure retained; reciprocal contact/bond slip law cannot be assumed', frozen_prediction_sha256=sha(p), wall_seconds=time.perf_counter() - start)
    save('raw/R3.json', r)
    state('R3_COMPLETE', '3-query reciprocal identity verified; physical gate NOT_RUN', 'preserve negative physical gate and complete independent numerical/adversarial checks and lab package')
    print(json.dumps(r, indent=2))
if __name__ == '__main__':
    run()
