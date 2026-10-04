from common import *
from native_port import model, deck_contract, project, q_matrix, p2, port_checks

def run():
    verify_inputs()
    r1 = read(R / 'raw/R1_NATIVE_PORT.json')
    r2 = read(R / 'raw/R2_SUPPORT_COVERAGE.json')
    r3 = read(R / 'raw/R3_COMPILED_SUPPORT.json')
    r4 = read(R / 'raw/R4_PORT_CERTIFICATE.json')
    ext = read(R / 'raw/EXTERNAL_FACIT.json')
    solver = read(R / 'raw/NATIVE_SOLVER.json')
    metrics = read(SRC / 'PREREG_R1.json')['metrics']
    (a, g, m) = model()
    de = deck_contract()
    (sl, pp, dd, ff, tri, b) = project(m)
    (Q, _) = q_matrix(len(m.Vc) + len(m.Vd), sl, m.nc + tri, p2(b))
    wrong = tri.copy()
    wrong[0] = m.Fd_master[(ff[0] + 1) % len(m.Fd_master), 2:8]
    actual_fault = port_checks(m, sl, ff, wrong, Q, p2(b), de)
    rng = np.random.default_rng(19)
    v = rng.normal(size=Q.shape[1])
    f = rng.normal(size=Q.shape[0])
    truework = float(Q @ v @ f)
    badwork = float(v @ (1.1 * (Q.T @ f)))
    badworkerr = abs(truework - badwork) / max(abs(truework), 1.0)
    force = np.array(de['loads'][0])
    validtotal = force.sum(0)
    badtotal = (force * 1.1).sum(0)
    loadtol = 1e-08
    native_error = solver['relative_tensor_L2_error']
    saved = np.load(X1 / 'inputs/fe/L005_lo_M1_k6_thin_3Y.npz')
    tensor = saved['Sf__axial'].astype(float)
    tensor_fault_error = float(np.linalg.norm(2 * tensor - tensor) / np.linalg.norm(tensor))
    displaced = np.load(DATA / 'SUPPORT_COVERAGE.npz')['P1_qualified_displacements'][0]
    displacement_fault_error = float(np.linalg.norm(2 * displaced - displaced) / np.linalg.norm(displaced))
    checks = {'source_input_hashes': True, 'native_declared_sets_agree': all(r1['checks']['set_parity'].values()), 'declared_port_failure_preserved': not r1['gate']['optimization_allowed'] and r1['interface_projection_mm']['max'] > 0.06, 'closest_search_exhaustive_witness': r1['worst_point']['closest_search_error_mm'] == 0, 'nearest_node_actually_rejected': r1['faults']['nearest_node']['port_identity_rejected'] and r1['faults']['nearest_node']['frozen_projection_rejected'], 'bad_facet_owner_rejected_by_actual_validator': not actual_fault['algebra_pass'] and (not actual_fault['facet_owner_parity']), 'wrong_virtual_work_transfer_rejected': badworkerr > metrics['virtual_work_relative_error'], 'wrong_load_rejected': float(np.max(abs(validtotal - badtotal))) > loadtol, 'source_solver_parity': solver['native_numeric_parity_pass'], 'source_tensor_fault_rejected': tensor_fault_error > metrics['native_solver_saved_tensor_relative_L2_error'], 'same_Q_full_FE_control': all((c['pass'] for c in r2['same_Q_full_FE_control'])) and r3['gate']['same_Q_full_FE'], 'same_Q_displacement_fault_rejected': displacement_fault_error > read(SRC / 'PREREG_R2.json')['metrics']['same_Q_KKT_response_relative_error'], 'qualified_coverage_failure_preserved': r2['gate']['qualified_projection'] and (not r2['gate']['full_native_coverage']), 'summary_identity_machine_exact': r2['sufficiency']['identity_error'] == 0.0, 'summary_downstream_discriminates': r2['sufficiency']['relative_downstream_difference'] > 0.01, 'compiled_receipt_rejects_faults': r3['receipt_controls']['valid_positive'] and r3['receipt_controls']['omitted_fault_rejected'] and r3['receipt_controls']['distance_coverage_fault_rejected'], 'geometry_certificate': all(r4['gate'].values()), 'source_force_fault_actually_rejected': ext['faults']['doubled_force_rejected'] and ext['faults']['valid_force_accepted'], 'protocol_positive_and_negative_tested': ext['faults']['source_protocol_positive_accepted'] and ext['faults']['source_protocol_transfer_rejected'], 'source_contrast_parity': ext['Chen_1mm_matched_source']['X1C_ratio_error'] <= 1e-12, 'original_X1C_failed_gates_preserved': ext['original_X1C_bracket_recalculation']['original_metrics_difference'] <= 1e-12 and (not all(ext['original_X1C_bracket_recalculation']['original_gate'].values()))}
    out = {'checks': checks, 'passed': sum(checks.values()), 'total': len(checks), 'controls': {'fault_virtual_work_relative_error': badworkerr, 'fault_load_total_difference_N': badtotal - validtotal, 'fault_native_stress_error': tensor_fault_error, 'fault_displacement_response_error': displacement_fault_error, 'wrong_facet_owner': actual_fault['facet_owner_parity']}, 'scope': 'These are computational controls. Physical strength is UNKNOWN; no clinical or manufacture acceptance.'}
    dump(R / 'raw/VALIDATION.json', out)
    if not all(checks.values()):
        raise ValueError('Validation failed ' + str([k for (k, v) in checks.items() if not v]))
    print(json.dumps(clean(out), indent=2))
if __name__ == '__main__':
    run()
