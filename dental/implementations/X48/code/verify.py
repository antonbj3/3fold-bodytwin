"""Recompute scientific gates from current values and exercise their rejection."""
from common import *
import time
from run_r2 import fraction_interval, rms
from run_r4 import control, projection
from sufficiency import run as sufficient

def run():
    start = time.monotonic()
    r2 = read(ROOT / 'raw/R2_RELATIVE_POSE.json')
    r4 = read(ROOT / 'raw/R4_TYPE_ENCLOSURES.json')
    checks = {}
    assert len(r2) == 1728 and len(r4) == 6912
    checks['r2_geodesic_tracking_current_values'] = sum((x['angular_tracking_gap_deg'] > 8 for x in r2)) == read(ROOT / 'raw/R2_RESULTS.json')['conditional_angular_tracking_rejections']
    checks['r2_lag_current_values'] = sum((x['observed_angular_deg'] + 4 < x['planned_angular_deg'] - 4 for x in r2)) == read(ROOT / 'raw/R2_RESULTS.json')['conditional_angular_lag_rows']
    checks['r2_direct_control_current_values'] = max((x['direct_point_error_mm'] for x in r2)) <= 1e-09 and max((x['direct_rotation_error_deg'] for x in r2)) <= 1e-09
    checks['r2_gauge_current_values'] = max((x['gauge_point_error_mm'] for x in r2)) <= 1e-09
    checks['r2_direct_control_injected_point_fails'] = not max((x['direct_point_error_mm'] for x in r2)) + 1 <= 1e-09
    checks['r2_direct_control_injected_angle_fails'] = not max((x['direct_rotation_error_deg'] for x in r2)) + 10 <= 1e-09
    checks['r4_ratio_box_current_values'] = all((x['conditional_ratio_interval'] == ratio_box(x['achieved'], x['planned'], x['achieved_radius'], x['planned_radius']) if x['achieved_radius'] is not None and x['planned_radius'] is not None else x['conditional_ratio_interval'] is None for x in r4))
    checks['r4_tracking_current_values'] = sum((x['achieved_radius'] is not None and x['planned_radius'] is not None and (abs(x['achieved'] - x['planned']) > x['achieved_radius'] + x['planned_radius']) for x in r4)) == read(ROOT / 'raw/R4_RESULTS.json')['conditional_tracking_rejections']
    n = control()
    checks['r4_finite_enclosure'] = n['finite_chord_disk_enclosure_pass'] and n['direct_quaternion_axis_error_deg'] <= 1e-09 and n['wrong_plus10deg_rejected']
    checks['summary_exact_identity_and_different_answer'] = sufficient()['scalar_motion']['identity_error'] == 0 and sufficient()['scalar_motion']['downstream_error_difference_mm'] == 0.5
    checks['angular_summary_exact_identity_and_different_answer'] = sufficient()['scalar_rotation']['identity_error'] == 0 and abs(sufficient()['scalar_rotation']['tracking_gap_B_deg'] - 20) < 1e-12
    rng = np.random.default_rng(48099)
    A = rng.normal(size=(16, 3))
    P = rng.normal(size=(16, 3))
    (ea, ep) = (0.1, 0.2)
    box = fraction_interval(A, P, ea, ep)
    inside = []
    for k in range(256):
        dA = rng.normal(size=A.shape)
        dA *= ea * rng.random() / rms(dA)
        dP = rng.normal(size=P.shape)
        dP *= ep * rng.random() / rms(dP)
        (a, p) = (A + dA, P + dP)
        alpha = float(np.mean(np.sum(a * p, axis=1)) / rms(p) ** 2)
        inside.append(box[0] - 1e-12 <= alpha <= box[1] + 1e-12)
    checks['fraction_nonlinear_ball_bound_covers'] = all(inside)
    checks['fraction_wrong_value_rejected'] = not box[0] <= box[1] + 1 <= box[1]
    checks['denominator_zero_refused'] = ratio_box(1.0, 0.0, 0.1, 0.1) is None
    frozen = read(ROOT / 'FROZEN_PREDICTIONS_R3.json')
    checks['frozen_forces_current_hash'] = sha(ROOT / frozen['prediction_file']) == frozen['sha256']
    force = read(ROOT / frozen['prediction_file'])
    ff = []
    wrong = []
    geom = {(a['patient'], a['jaw']): a for a in force['arches']}
    for row in force['rows']:
        if 'wrenches' not in row:
            continue
        a = geom[row['patient'], row['jaw']]
        F = sum((np.array(v[:3]) for v in row['wrenches'].values()), np.zeros(3))
        M = sum((np.array(v[3:]) + np.cross(a['teeth'][k]['centroid_mm'], v[:3]) for (k, v) in row['wrenches'].items()), np.zeros(3))
        ff.append(np.linalg.norm(F) <= 1e-08 and np.linalg.norm(M) <= 1e-08)
        wrong.append(np.linalg.norm(F + np.array([1.0, 0.0, 0.0])) > 1e-08 and np.linalg.norm(M + np.array([1.0, 0.0, 0.0])) > 1e-08)
    checks['wrench_current_balance'] = all(ff)
    checks['wrench_actual_injection_rejected'] = all(wrong)
    reg = read(ROOT / 'raw/REGISTRATION_AUDIT.json')
    checks['raw_external_STL_recomputed'] = len(reg['raw_members']) == 8 and all((x['points_recomputed_from_external_STL'] for x in reg['raw_members']))
    checks['independent_point_plane_closed_form'] = all((x['passed'] and x['wrong_1mm_rejected'] and x['wrong_10deg_rejected'] for x in reg['controls']))
    for i in range(1, 5):
        checks[f'prereg_R{i}_hash'] = sha(ROOT / f'PREREG_R{i}.json') == read(ROOT / f'PREREG_R{i}.sha256.json')['sha256']
    lit = read(ROOT / 'raw/LITERATURE_MOTION_REFERENTS.json')
    checks['wrong_quantity_bank_rejected'] = lit['supplied_bank']['usable_movement_records'] == 0 and lit['supplied_bank']['rejected_movement_records'] == 191
    checks['table_missing_cells_not_imputed'] = len(lit['source_unavailable_cells']) == 11 and len(lit['table_rows']) == 31 and (lit['direct_validation_eligible_rows'] == 0)
    checks['no_false_body_or_anatomic_claim'] = all((x['body'] == 'UNKNOWN_ROOT_NOT_OBSERVED' and x['fixed_reference'] is False and (x['input_radii_empirically_validated'] is False) for x in r4))
    result = dict(status='PASS' if all(checks.values()) else 'FAIL', checks=checks, scope='source integrity is recomputed separately at each run; these are numerical/gate checks only, not physical validation', wall_s=time.monotonic() - start)
    write(ROOT / 'raw/VERIFICATION.json', result)
    print(result)
    assert all(checks.values())
    return result
if __name__ == '__main__':
    run()
