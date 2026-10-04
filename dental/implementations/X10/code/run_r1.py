from common import *
import math
start = time.perf_counter()
pr = frozen('PREREG_R1.json')
g = load('raw/global_scanners.json')['values']
targets = load('raw/clinical_gaps.json')
pred = [{'scanner': t['scanner'], 'predicted_group_gap_um': 40 + g[t['scanner']]['mean_um'], 'nominal_practice_um': 40.0, 'transferred_bound_claim': 'HYPOTHESIS_ONLY'} for t in targets]
freeze_once('FROZEN_PREDICTIONS.json', {'construction': 'R1', 'inputs': {'prereg': sha(ROOT / 'PREREG_R1.json'), 'global_values': sha(ROOT / 'raw/global_scanners.json')}, 'predictions': pred, 'scope': 'Retrospective known published group-summary comparison; frozen before evaluation, not before historical clinical measurements.'})
rows = []
for (p, t) in zip(pred, targets):
    control = sum([40.0, g[t['scanner']]['mean_um']])
    rows.append(dict(t, **p, error_um=abs(p['predicted_group_gap_um'] - t['value_um']), false_group_accept=p['predicted_group_gap_um'] < 120 <= t['value_um'], control_um=control, control_agreement=abs(control - p['predicted_group_gap_um']) <= 1e-09, actual_certificate='ABSTAIN_UNKNOWN: signed local error, process calibration, tails and version match absent'))
imp = load('raw/implant_medians.json')
mapname = {'Prime': 'Primescan', 'Omni': 'Omnicam 4.6', 'Trios 3': 'Trios 3', 'Trios 4': 'Trios 4'}
checks = []
for (name, values) in imp['position_3d_um'].items():
    bound = g[mapname[name]]['mean_um']
    for (j, v) in enumerate(values):
        checks.append({'scanner': name, 'scanbody': j + 1, 'proposed_bound_um': bound, 'measured_median_um': v, 'excess_um': v - bound, 'bound_holds': v <= bound + 0.01})
bad_target = targets[0]['value_um'] + 200
injected_accuracy_fail = abs(pred[0]['predicted_group_gap_um'] - bad_target) > 20

def within_bound(observed, bound):
    return math.isfinite(observed) and math.isfinite(bound) and (observed <= bound + 0.01)
bound = checks[0]['proposed_bound_um']
injected_bound_fail = within_bound(bound, bound) and (not within_bound(bound + 200, bound))
result = {'construction': 'R1', 'outcome': 'FAIL_GLOBAL_TO_GAP_TRANSFER; TIE_SCALAR_CONTROL', 'clinical_rows': rows, 'max_error_um': max((r['error_um'] for r in rows)), 'wrong_group_summary_accepts': sum((r['false_group_accept'] for r in rows)), 'implant_bound_rows': checks, 'implant_enclosure_failures': sum((not c['bound_holds'] for c in checks)), 'all_actual_certificates': 'ABSTAIN_UNKNOWN', 'equal_info_control_tie': all((r['control_agreement'] for r in rows)), 'injected_error_checks': {'clinical_plus_200_um_fails': injected_accuracy_fail, 'implant_plus_200_um_fails': injected_bound_fail}, 'external_referent': dict(pr['external_referent'], refutes_us=True), 'limitations': ['This refutes a universal cross-regime shortcut, not any manufacturer claim in its measured regime.', 'Global surface means and landmark medians differ in estimand: cannot establish same-observable metrological calibration.', 'Crown comparison uses cross-year nominal aliases; no version-equivalence claim; actual certificate refuses.', '2 clinical 2020 numeric rows remain secondary-table transcriptions.'], 'query_and_validation_seconds': time.perf_counter() - start}
write('RESULTS_R1.json', result)
(ROOT / 'HANDOFF_R1.md').write_text('R1 complete: global scanner summary cannot certify local seated gap. ' + f"Maximum retrospective clinical summary error {result['max_error_um']:.1f} um; {result['wrong_group_summary_accepts']} wrong summary acceptance; {result['implant_enclosure_failures']}/28 transferred implant bounds fail. " + 'Scalar control TIE. All actual certificates UNKNOWN. Next: use published scanner×scanbody geometry; do not invent a signed field from RMS.\n')
state('R1_COMPLETE', result['outcome'], 'freeze R2 scanner×scanbody operation with measured position and angular channels')
print(result['outcome'], 'max clinical error', result['max_error_um'], 'implant failures', result['implant_enclosure_failures'])
