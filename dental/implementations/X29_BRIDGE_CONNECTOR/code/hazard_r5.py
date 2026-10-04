import numpy as np
from common import R, read, dump, state

def hazard(s, area, m):
    return float(np.sum(area * np.maximum(s, 0.0) ** m))

def rates(tag, lc, m):
    base = R / 'raw/fe' / tag
    g = np.load(base / 'geometry.npz')
    s = np.load(base / (lc + '_stress.npz'))['s1']
    cf = g['cf']
    r = read(f'raw/fe/{tag}/FE.json')
    connector = np.abs(np.abs(cf[:, 0]) - 5.5) < 1.25
    contact = np.zeros(len(cf), bool)
    diag = next((x for x in r['load_cases'] if x['name'] == lc))['load']
    for p in diag['patches']:
        contact |= (np.linalg.norm(cf - p['p'], axis=1) < 3 * p['a_mm']) & (np.abs(cf[:, 0]) < 4.0)
    return [hazard(s[mask], g['af'][mask], m) for mask in [connector, contact]]

def run():
    p = read('PREREG_R5.json')
    contrasts = []
    for m in p['m_panel']:
        for lc in ['axial', 'offaxis30']:
            vals = []
            for spacing in [0.24, 0.18]:
                ref = rates(f'h4_b2.25_s{spacing:g}', lc, m)
                alt = rates(f'h4_b4_s{spacing:g}', lc, m)
                (qc, qp) = (np.array(ref) / np.array(alt)) ** (1 / m)
                vals.append(dict(spacing_mm=spacing, q_connector=float(qc), q_contact=float(qp), conditional_contact_origin_interval=[float(min(qc, qp)), float(qp)]))
            err = {k: abs(vals[0][k] - vals[1][k]) / vals[1][k] for k in ['q_connector', 'q_contact']}
            contrasts.append(dict(m=m, load_case=lc, mesh_values=vals, mesh_relative_errors=err, gate_pass=all((v <= 0.05 for v in err.values()))))
    s = np.array([2.0, 2.0])
    a = np.array([3.0, 5.0])
    m = 10
    expected = 8 * 2.0 ** m
    actual = hazard(s, a, m)
    fault1 = hazard(s, np.ones_like(a), m)
    fault2 = (actual / expected) ** (1 / (m + 1))
    q_correct = 32.0 ** (1 / m)
    q_wrong = 32.0 ** (1 / (m + 1))
    out = dict(claim_type='capability', contrasts=contrasts, all_mesh_gate_pass=all((c['gate_pass'] for c in contrasts)), analytic_integral_control=dict(expected=expected, actual=actual, pass_gate=actual == expected, missing_area_fault_rejected=fault1 != expected, wrong_root_fault_rejected=not np.isclose(q_correct, q_wrong)), external_referent=p['external_referent'], external_math_check=dict(kind='closed_form', locator='https://doi.org/10.1115/1.4010337; weakest-link equation4; constant-stress integral I=A*sigma^m is our derived reduction', compared_quantity='integral dimensions/area and inverse-mroot, not experimental strength', refutes_us=False), empirical_status='UNKNOWN_NO_MATCHED_GEOMETRY_NO_WEIBULL_BATCH_NO_CONTACTLAW', parent_failed_gates=['R1rank/studies', 'R2absolute regionalstress', 'R4pairedp99'], inherited_operator='X1B/crown_design_fe.py postprocess_arrays IS_m; no newalgorithm claim')
    dump('raw/R5_HAZARD.json', out)
    state('R5_complete', 'surfacehazard ratioPASS' if out['all_mesh_gate_pass'] else 'surfacehazard ratioFAIL', 'freeze diagnostic prospective output; make one-command demo; next physicalconstruction measuredroundedconnector andlocalstrain')
    print('R5 surfacehazard ratios', out['all_mesh_gate_pass'])
    print([(x['m'], x['load_case'], x['mesh_relative_errors']) for x in contrasts])
    return out
if __name__ == '__main__':
    run()
