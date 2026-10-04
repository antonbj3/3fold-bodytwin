"""Local reviewable adapters to Field TANDLAST and crown-design consumers.
No source checkout or source graph is changed. Forces stay at PER_TOOTH.
"""
from dental_release.paths import expand as _release_expand
import json
from pathlib import Path
R = Path(__file__).resolve().parents[1]

def design_request(s, a, fdi, geometry_sha256):
    if geometry_sha256 != s['geometry_sha256']:
        raise ValueError('Mechanics/geometry identity mismatch')
    if a.get('force_interval_N') is None:
        return {'status': 'UNKNOWN', 'reason': 'No bounded force response'}
    i = s['predicted_fdi'].index(fdi)
    return {'schema': 'crown-height-force-input-v1', 'status': 'CONDITIONAL_LAB_INPUT', 'case': s['case'], 'fdi': fdi, 'mechanics_geometry_sha256': geometry_sha256, 'height_change_mm': a['height_change_mm'][i], 'force_N': a['force_N'][i], 'force_interval_N': a['force_interval_N'][i], 'neighbor_force_intervals_N': {str(t): a['force_interval_N'][j] for (j, t) in enumerate(s['predicted_fdi']) if abs(t - fdi) == 1}, 'joint_force_state': {'A': s['A'], 'w_N': a['w_N'], 'f_N': a['force_N'], 'joint_l2_error_radius_N': a['joint_l2_error_radius_N']}, 'resolution': 'PER_TOOTH', 'timescale': 'SIMULTANEOUS', 'receiving_functions': ['cells/design/crown_design_pipeline.py:evaluate', 'cells/design/crown_design_fe.py:load_case', 'LANE_X18_CROWN_ANTAGONIST/code/round2.py:predict'], 'application_contract': 'Total tooth force supplies load magnitude. A separately bound spatial load shape is still required for stress/pressure; load_case returns per-1N node loads, F_ref controls the assumed Hertz patch. Do not substitute F_ref for an actual force scale.', 'generator_action': 'Apply the declared scalar height change along the calibrated axial action on the same crown surface. Recheck normals, force arms, fit/clearance and reference identity; a new surface patch requires a new measured mechanical action.', 'physical_validity': 'UNKNOWN; simulated leveled instrument, no clinical or manufacturing approval'}

def field_request(s, a, gap_convention):
    if gap_convention != 'MEASURED_LEVELED_REFERENCE':
        return {'status': 'UNKNOWN', 'reason': 'Native X21 gap classification and leveled instrument do not share a loaded contact state; no patient-force transfer', 'resolution': 'PER_TOOTH'}
    return {'status': 'CONDITIONAL_REGION_RESPONSE', 'case': s['case'], 'resolution': 'PER_TOOTH', 'timescale': 'SIMULTANEOUS', 'source': _release_expand('FALT_TANDLAST'), 'force_intervals_N': dict(zip(map(str, s['predicted_fdi']), a['force_interval_N'])), 'scope': 'Upper support regions only, same registered A/w. This replaces the arbitrary-positive-support premise with measured response information; it does not revise native MUST/CAN/NEVER or lower-arch forces.'}

def main():
    s = json.loads((R / 'inputs/F5367_R4_state.json').read_text())
    inv = json.loads((R / 'exports/F5367_inverse.json').read_text())
    a = inv['at_selected_height']
    out = {'design': design_request(s, a, 16, s['geometry_sha256']), 'field_native_guard': field_request(s, a, 'NATIVE_X21_GAPS'), 'field_conditional_port': field_request(s, a, 'MEASURED_LEVELED_REFERENCE')}
    try:
        design_request(s, a, 16, 'wrong')
        fault = False
    except ValueError:
        fault = True
    out['geometry_mismatch_fault_detected'] = fault
    (R / 'exports/CROWN_AND_FIELD_REQUEST.json').write_text(json.dumps(out, indent=2))
    print(json.dumps({'design': out['design']['status'], 'native_transfer': out['field_native_guard']['status'], 'geometry_fault': fault}))
if __name__ == '__main__':
    main()
