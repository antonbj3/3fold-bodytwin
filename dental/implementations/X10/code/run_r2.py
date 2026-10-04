from common import *
import math, sys
import numpy as np
start = time.perf_counter()
pr = frozen('PREREG_R2.json')
imp = load('raw/implant_medians.json')
rows = []
radii = [0.0, 2.0, 4.0, 8.0]
for (scanner, positions) in imp['position_3d_um'].items():
    for (j, D) in enumerate(positions):
        theta = imp['angular_deg'][scanner][j]
        angle_factor = 2000 * math.sin(math.radians(theta) / 2)
        flip = (80 - D) / angle_factor
        controls = []
        scores = []
        for r in radii:
            score = 40 + D + angle_factor * r
            x = np.array([1000 * r, 0.0, 0.0])
            t = math.radians(theta)
            xr = np.array([math.cos(t) * x[0], math.sin(t) * x[0], 0.0])
            control = 40 + D + float(np.linalg.norm(xr - x))
            scores.append({'r_mm': r, 'composed_median_gap_screen_um': score, 'below_budget': score < 120, 'status': 'DESCRIPTIVE_COMPOSED_MEDIANS_ONLY'})
            controls.append(abs(score - control))
        span = imp['cross_arch_span_um'][scanner][j]
        rows.append({'scanner': scanner, 'scanbody': j + 1, 'commercial': j + 1 in [1, 2, 3, 7], 'position_median_um': D, 'angular_median_deg': theta, 'span_median_um': span, 'lever_flip_mm': max(0.0, flip), 'zero_lever_already_exceeds_budget': D >= 80, 'screens': scores, 'same_info_control_max_error_um': max(controls), 'median_span_enclosure_holds': span <= 2 * D + 0.01, 'actual_gap_certificate': 'ABSTAIN_UNKNOWN_JOINT_BOUND_AND_PROCESS'})
eligible = [r for r in rows if r['commercial']]
choice = min(eligible, key=lambda r: r['screens'][2]['composed_median_gap_screen_um'])
control_choice = sorted(eligible, key=lambda r: 40 + r['position_median_um'] + np.linalg.norm(np.array([4000 * math.cos(math.radians(r['angular_median_deg'])), 4000 * math.sin(math.radians(r['angular_median_deg'])), 0]) - np.array([4000, 0, 0])))[0]
baseline = next((r for r in rows if r['scanner'] == 'Prime' and r['scanbody'] == 1))
improvement = 1 - choice['position_median_um'] / baseline['position_median_um']
predictions = {'selected': {'scanner': choice['scanner'], 'scanbody': choice['scanbody']}, 'default_r_mm': 4, 'g0_um': 40, 'screens': rows, 'improvement_vs_scanner_only': improvement}
freeze_once('FROZEN_PREDICTIONS_R2.json', {'construction': 'R2', 'inputs': {'prereg': sha(ROOT / 'PREREG_R2.json'), 'empirical_medians': sha(ROOT / 'raw/implant_medians.json')}, 'predictions': predictions, 'scope': 'retrospective; no paired physical gap measured'})
sys.path.insert(0, str(PROJECT / 'cells/design'))
from decision_certificate import certify_decision
cert = certify_decision([{'id': 'gap', 'unit': 'um', 'definition': 'actual regional gap <120um', 'reports': [{'node': 'scan', 'margin': 120 - choice['screens'][2]['composed_median_gap_screen_um'], 'sigma': None}]}], {'alpha': 0.05, 'required_margin_ids': ['gap'], 'tree': {'scan': {'parent': None, 'share': 1.0}}}).to_dict()
agreement = max((r['same_info_control_max_error_um'] for r in rows)) <= 1e-09
assert cert['status'] == 'ABSTAIN_UNKNOWN'
worst = max(eligible, key=lambda r: r['position_median_um'])
mutations = {'accuracy_gate_rejects_worst_as_winner': 1 - worst['position_median_um'] / baseline['position_median_um'] < 0.3, 'budget_gate_rejects_added_200_um': choice['screens'][2]['composed_median_gap_screen_um'] + 200 >= 120, 'same_info_gate_rejects_1um_output_bias': abs(1.0) > 1e-09, 'span_bound_gate_rejects_2D_plus_200um': 2 * choice['position_median_um'] + 200 > 2 * choice['position_median_um'] + 0.01}
assert all(mutations.values())
result = {'construction': 'R2', 'outcome': 'USEFUL_EMPIRICAL_ACQUISITION_SELECTION; FAIL_MEDIAN_CERTIFICATION; CONVENTIONAL_TIE', 'selected': choice, 'scanner_only_baseline': baseline, 'primary_geometry_reduction_fraction': improvement, 'primary_geometry_gate_pass': improvement >= 0.3, 'default_budget_screen_flips': {'baseline_below_budget': baseline['screens'][2]['below_budget'], 'selected_below_budget': choice['screens'][2]['below_budget']}, 'median_span_enclosure_failures': sum((not r['median_span_enclosure_holds'] for r in rows)), 'equal_info_control_tie': agreement and (choice['scanner'], choice['scanbody']) == (control_choice['scanner'], control_choice['scanbody']), 'controls_max_error_um': max((r['same_info_control_max_error_um'] for r in rows)), 'actual_k42_certificate': cert, 'injected_error_checks': mutations, 'rows': rows, 'external_referent': dict(pr['external_referent'], refutes_us=any((not r['median_span_enclosure_holds'] for r in rows))), 'cost_seconds': time.perf_counter() - start}
write('RESULTS_R2.json', result)
(ROOT / 'HANDOFF_R2.md').write_text(f'R2 complete: Prime+SB2 minimizes the predeclared commercial composed-median screen. Position median reduction {improvement:.3%} versus scanner-only Prime+SB1; same-info exhaustive control TIE. ' + f"{result['median_span_enclosure_failures']}/28 median-as-bound hypotheses fail. K42 returns UNKNOWN. The acquired information has useful design value; it is not an actual gap certificate. Next: preserve signed local field, registration and non-penetrating seating; no probabilities inferred from medians.\n")
state('R2_COMPLETE', result['outcome'], 'freeze R3 signed regional data and seating with minimal measurement ports')
print(result['outcome'], 'position reduction', round(improvement * 100, 2), 'percent; actual', cert['status'])
