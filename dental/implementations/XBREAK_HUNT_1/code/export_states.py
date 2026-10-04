from bench import *
import numpy as np

def main():
    r4 = json.loads((ROOT / 'raw/R4.json').read_text())
    r5 = json.loads((ROOT / 'raw/R5.json').read_text())
    for row in r4['rows']:
        state = {'schema': 'conditional-edit-response-v1', 'case': row['case'], 'geometry_sha256': row['source']['sha256'], 'predicted_fdi': row['source']['predicted_fdi'], 'A': row['A'], 'basis_N': row['basis_N'], 'H_N_per_mm': row['observed_H_N_per_mm'], 'w_N': row['w_N'], 'baseline_force_N': row['observed_f0_N'], 'baseline_height_mm': [0.0] * row['n'], 'radius_model': 'full_response', 'spectral_error_N_per_mm': row['rho_N_per_mm'], 'baseline_error_N': row['eta0_N'], 'complete_active_calibration': row['acquisition_gate'], 'provenance_kind': 'SIMULATED_PROSPECTIVE_ACQUISITION', 'source_round': 'R4', 'source_result_sha256': sha(ROOT / 'raw/R4.json'), 'measurement_contract': 'Same fixed A/wrench, passive linear support, all candidate contacts initially loaded; raw regional force error <=0.01 N, exact height and wrench assumed'}
        write(ROOT / 'exports' / f"{row['case']}_R4_state.json", state)
    source = r5['rows'][0]
    world = source['worlds'][0]
    A = np.array(source['A'])
    from force_operator import balanced_basis
    state = {'schema': 'conditional-edit-response-v1', 'case': source['case'], 'geometry_sha256': source['source']['sha256'], 'predicted_fdi': source['source']['predicted_fdi'], 'A': source['A'], 'basis_N': balanced_basis(A).tolist(), 'H_N_per_mm': world['H_observed_N_per_mm'], 'w_N': source['w_N'], 'baseline_force_N': world['baseline_force_N'], 'baseline_height_mm': world['baseline_height_mm'], 'radius_model': 'activated_rank_one', 'spectral_error_N_per_mm': world['new_radius_N_per_mm'], 'old_spectral_error_N_per_mm': world['old_radius_N_per_mm'], 'activation_column_N_per_mm': world['v_N_per_mm'], 'baseline_error_N': 0.01 * np.sqrt(source['n']), 'complete_active_calibration': world['new_acquisition_gate'], 'provenance_kind': 'SIMULATED_PROSPECTIVE_ACQUISITION', 'source_round': 'R5/R6', 'source_result_sha256': sha(ROOT / 'raw/R5.json'), 'measurement_contract': 'R5 two directed 50/90 um activation levels, fixed A/wrench and same linear branch; raw regional force error <=0.01 N'}
    write(ROOT / 'exports' / f"{source['case']}_R5_state.json", state)
    edit = {'case': source['case'], 'geometry_sha256': source['source']['sha256'], 'unit': 'mm', 'height_change_mm': [-0.04] + [0.0] * (source['n'] - 1)}
    write(ROOT / 'exports' / f"{source['case']}_edit.json", edit)
    print('12 full-response states and one activated state exported; true supports/compliance excluded')
if __name__ == '__main__':
    main()
