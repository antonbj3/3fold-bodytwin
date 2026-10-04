import numpy as np
from common import R, read, dump, state
from bridge_fe import solve

def capacity_interval(qc, qp):
    return [min(qc, qp), qp]

def run(solve_new=True):
    old = read('raw/BRIDGE_FE.json')
    new = []
    for spacing in [0.24, 0.18]:
        tag = f'h4_b4_s{spacing:g}'
        new.append(solve(4.0, 4.0, spacing, tag) if solve_new else read(f'raw/fe/{tag}/FE.json'))
    p = read('PREREG_R4.json')
    contrasts = []
    for spacing in [0.24, 0.18]:
        ref = next((r for r in old['runs'] if r['height_mm'] == 4.0 and r['field_spacing_mm'] == spacing))
        alt = next((r for r in new if r['field_spacing_mm'] == spacing))
        for (i, name) in enumerate(['axial', 'offaxis30']):
            A = ref['load_cases'][i]['regions']
            B = alt['load_cases'][i]['regions']
            kc_ref = max((A[k]['p99_MPa_per_N'] for k in ['connector_left', 'connector_right']))
            kc_new = max((B[k]['p99_MPa_per_N'] for k in ['connector_left', 'connector_right']))
            qp = A['pontic_contact']['p99_MPa_per_N'] / B['pontic_contact']['p99_MPa_per_N']
            qc = kc_ref / kc_new
            contrasts.append(dict(load_case=name, spacing_mm=spacing, q_connector=qc, q_contact=qp, contact_origin_conditional_capacity_ratio_interval=capacity_interval(qc, qp), crossover_x=qp / qc, interpretation='conditional proxy model only; specimen-level origin/batch/contact law unmeasured'))
    checks = []
    for lc in ['axial', 'offaxis30']:
        (a, b) = [c for c in contrasts if c['load_case'] == lc]
        for k in ['q_connector', 'q_contact']:
            e = abs(a[k] - b[k]) / b[k]
            checks.append(dict(load_case=lc, quantity=k, relative_mesh_difference=e, gate_pass=e <= 0.05))
    (qc, qp) = (2.0, 1.1)
    grid = np.geomspace(1, 100, 500)
    vals = np.minimum(qc * grid, qp)
    iv = capacity_interval(qc, qp)
    algebra = bool(vals.min() >= iv[0] and vals.max() <= iv[1] and np.isclose(vals[-1], iv[1]))
    wrong_inverse = capacity_interval(1 / qc, 1 / qp)
    faults = [dict(name='wrong_inverse_ratio', rejected=not (np.isclose(wrong_inverse[0], iv[0]) and np.isclose(wrong_inverse[1], iv[1]))), dict(name='swapped_mode_capacity', rejected=not np.isclose(min(1.2 * qc, qp), min(1.2 * qp, qc)))]
    out = dict(claim_type='capability', new_runs=new, contrasts=contrasts, paired_mesh_checks=checks, paired_mesh_gate_pass=all((c['gate_pass'] for c in checks)), algebra_control_pass=algebra, algebra_faults=faults, absolute_stress_parent_gate='FAILED_RETAINED; ratio gate does not repair absolute-stress gate', empirical_validation='UNKNOWN; no independent published16mm2pairedbridge geometry; no physical test performed', external_referent=p['external_referent'], capability_status='COMPUTATIONAL_MODE_ENVELOPE_ONLY', debt='x>=1 requires confirmed specimen-level reference pontic origin. Contact stress proxy must be physically calibrated. Mesh bracket is not statistical CI or rigorous discretization bound.')
    dump('raw/R4_PAIRED.json', out)
    state('R4_complete', 'paired geometry ratio PASS' if out['paired_mesh_gate_pass'] else 'paired geometry ratio FAIL', 'freeze prospective lab ratio envelope; assemble executable demo and missing graph coverage')
    print('R4 paired gate', out['paired_mesh_gate_pass'], contrasts, flush=True)
    return out
if __name__ == '__main__':
    run()
