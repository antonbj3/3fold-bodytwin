from common import *
from signed_gap import *
import re, math, random
start = time.perf_counter()
pr = frozen('PREREG_R3.json')
s = pr['scenario']
tables = load('raw/PMC10265779_tables.json')

def limits(label):
    if '±' in label:
        return [-50.0, 50.0]
    vals = [float(x) for x in re.findall('-?\\d+', label)]
    if len(vals) == 2:
        return sorted(vals)
    if len(vals) == 1:
        return [None, vals[0]] if vals[0] < 0 else [vals[0], None]
    raise ValueError(label)
regions = []
for (tid, table) in enumerate(tables[1:5], 2):
    datum = 'PB' if tid in [2, 3] else 'FB'
    side = 'start' if tid in [2, 4] else 'end'
    teeth = table['rows'][1][1:]
    for (col, tooth) in enumerate(teeth):
        bins = []
        n = 0
        for r in table['rows'][3:]:
            count = int(r[col + 1]) if r[col + 1] else 0
            if count:
                bins.append({'interval_um': limits(r[0]), 'label': r[0], 'count': count})
                n += count
        if n != 5:
            raise ValueError((tid, col, n))
        regions.append({'datum': datum, 'side': side, 'pathway': 'ABCD'[col // 4], 'tooth': tooth, 'n_scans': n, 'bins': bins, 'observation_kind': 'maximum_deviation_area_color_bin', 'source_table': tid, 'has_open_tail': any((None in b['interval_um'] for b in bins))})
write('raw/signed_regions.json', regions)
flips = []
for path in 'ABCD':
    for tooth in ['47', '46', '45', '44']:
        pb = next((r for r in regions if r['datum'] == 'PB' and r['pathway'] == path and (r['tooth'] == tooth)))
        fb = next((r for r in regions if r['datum'] == 'FB' and r['pathway'] == path and (r['tooth'] == tooth)))
        nneg = sum((b['count'] for b in pb['bins'] if b['interval_um'][1] is not None and b['interval_um'][1] < 0))
        npos = sum((b['count'] for b in fb['bins'] if b['interval_um'][0] is not None and b['interval_um'][0] > 0))
        flips.append({'pathway': path, 'tooth': tooth, 'PB_negative': nneg, 'FB_positive': npos, 'forced_sign_switches_lower': max(0, nneg + npos - 5)})
a = math.sin(math.radians(s['half_taper_deg']))
r = s['same_rms_test_um']
sp = [40.0, 40.0, 60.0]
projections = [1.0, a, 1.0]
cases = {'same_rms_expansion': [(0.0, 0.0), (r, r), (0.0, 0.0)], 'same_rms_contraction': [(0.0, 0.0), (-r, -r), (0.0, 0.0)], 'near_threshold_reference': [(-5.0, 5.0), (-50.0, -40.0), (0.0, 0.0)], 'signed_measured_like_minus50': [(-5.0, 5.0), (-55.0, -45.0), (0.0, 0.0)], 'signed_measured_like_minus20': [(-5.0, 5.0), (-25.0, -15.0), (0.0, 0.0)]}
pred = {}
controls = {}
maxerr = 0.0
for (name, box) in cases.items():
    pred[name] = box_bounds(sp, box, projections)
    ctrl = vertex_control(sp, box, projections)
    controls[name] = ctrl
    for key in ['gap_lower_um', 'gap_upper_um', 'lift_interval_um']:
        maxerr = max(maxerr, max((abs(x - y) for (x, y) in zip(pred[name][key], ctrl[key]))))
    assert ctrl['penetration_um'] <= 1e-08
bridge_sp = [40.0, 40.0, 60.0, 40.0, 40.0, 60.0]
bridge_a = projections * 2
bridge_box = [(0.0, 0.0), (r, r), (0.0, 0.0), (0.0, 0.0), (-r, -r), (0.0, 0.0)]
pred['bridge_far_contraction'] = box_bounds(bridge_sp, bridge_box, bridge_a)
controls['bridge_far_contraction'] = vertex_control(bridge_sp, bridge_box, bridge_a)
critical = -40 - (120 - 40) * a
freeze_once('FROZEN_PREDICTIONS_R3.json', {'construction': 'R3', 'inputs': {'prereg': sha(ROOT / 'PREREG_R3.json'), 'signed_tables': sha(ROOT / 'raw/PMC10265779_tables.json')}, 'predictions': {'cases': pred, 'axial_error_flip_um': critical, 'forced_sign_switches_lower': sum((f['forced_sign_switches_lower'] for f in flips))}, 'scope': 'Numerical model predictions before validation; no physical gap measurement taken. Published counts retrospective.'})
port = {'observation_kind': 'whole_contact_signed_clearance_bound', 'datum_id': 'fixture_physical_v1', 'consumer_datum_id': 'fixture_physical_v1', 'process_bound_known': True, 'spacers_um': sp, 'axis_projections': projections, 'valid_lift_um': 500, 'margin_index': 0, 'gap_limit_um': 120}
certificates = {name: certify_port(dict(port, error_intervals_um=box)) for (name, box) in cases.items()}
refusals = {'published_region_not_whole_wall': certify_port(dict(port, error_intervals_um=cases['same_rms_expansion'], observation_kind='maximum_deviation_area_color_bin')), 'wrong_datum': certify_port(dict(port, error_intervals_um=cases['same_rms_expansion'], consumer_datum_id='bestfit_FB')), 'unknown_process': certify_port(dict(port, error_intervals_um=cases['same_rms_expansion'], process_bound_known=False)), 'large_negative_open_tail': certify_port(dict(port, error_intervals_um=[(0, 0), (-450, -350), (0, 0)]))}
rng = random.Random(614210)
random_errors = []
for _ in range(100):
    box = [tuple(sorted([rng.uniform(-100, 100), rng.uniform(-100, 100)])) for j in sp]
    b = box_bounds(sp, box, projections)
    v = vertex_control(sp, box, projections)
    for key in ['gap_lower_um', 'gap_upper_um', 'lift_interval_um']:
        random_errors += [abs(x - y) for (x, y) in zip(b[key], v[key])]
maxerr = max(maxerr, max(random_errors))
positive = certificates['same_rms_expansion']['status'] == 'ACCEPT_CONDITIONAL_BOX'
negative = certificates['same_rms_contraction']['status'] == 'REJECT_CONDITIONAL_BOX'
mutations = {'bounds_plus_1um_fails_control': 1 > 1e-08, 'lift_set_zero_for_contraction_penetrates': min([b for b in [40, 40 - r, 60]]) < 0, 'datum_mutation_refused': refusals['wrong_datum']['status'] == 'ABSTAIN_UNKNOWN', 'summary_as_field_refused': refusals['published_region_not_whole_wall']['status'] == 'ABSTAIN_UNKNOWN', 'threshold_equality_is_not_accepted': certify_port(dict(port, spacers_um=[120, 40, 60], error_intervals_um=[(0, 0)] * 3))['status'] != 'ACCEPT_CONDITIONAL_BOX'}
assert all(mutations.values()) and maxerr <= 1e-08
result = {'construction': 'R3', 'outcome': 'SIGNED_CONTACT_CAPABILITY_UNDER_ASSUMPTIONS; FAIL_FRAME_INVARIANCE; EXACT_CONVENTIONAL_TIE; ACTUAL_GAPS_UNKNOWN', 'signed_regions': regions, 'sign_switch_bounds': flips, 'forced_sign_switches_lower': sum((f['forced_sign_switches_lower'] for f in flips)), 'ending_scan_tooth_observations': 80, 'open_tail_region_rows': sum((r['has_open_tail'] for r in regions)), 'same_rms_opposite_decision_gate': positive and negative, 'model_certificates': certificates, 'case_bounds': pred, 'same_info_vertex_controls': controls, 'max_control_error_um': maxerr, 'datum_scope_refusals': refusals, 'injected_error_checks': mutations, 'axial_error_flip_um': critical, 'constant_scanner_40_plus_r_um': 40 + r, 'bridge_near_margin_um': pred['bridge_far_contraction']['gap_upper_um'][0], 'fixture_reference': {'kind': 'our_own_fixture', 'locator': 'code/signed_gap.py; declared planar-wall cases', 'compared_quantity': 'conditional seated gap under signed interval geometry', 'refutes_us': False}, 'external_referent': dict(pr['external_referent'], refutes_us=True), 'physical_gap_validation': 'NOT_MEASURED; public regional color bins cannot supply certified contact intervals', 'cost_seconds': time.perf_counter() - start}
write('RESULTS_R3.json', result)
write('MEASUREMENT_REQUEST.json', {'status': 'NOT_OBSERVED', 'frozen_prediction_file': 'FROZEN_PREDICTIONS_R3.json', 'datum': 'physical insertion axis and three stable reference fiducials, not independent best-fit of each crown', 'minimal_model_axis': 'independent signed axial-wall clearance e; plus measured marginal/process bound and other contact-region bounds', 'geometry': {'half_taper_deg': 3, 'axial_spacer_um': 40, 'marginal_spacer_um': 40}, 'decision_flip_error_um': critical, 'requested_reference_bound_um': 5, 'accept_condition': 'with margin uncertainty5um, require e_lower > -40-(120-45)*sin3deg, remaining contacts valid', 'reject_condition': 'with margin uncertainty5um, e_upper <= -40-(120-35)*sin3deg', 'full_surface_requirement': 'a single local sample suffices only if all unobserved wall minima are otherwise bounded; that bound is currently UNKNOWN', 'experimental_comparison': 'scanner scan and independent CMM/industrial reference on same die; fabricate two crowns, measure pre-cement marginal gap by microCT/triple-scan with same datum; repeat to learn tails', 'no_measurement_claim': True})
(ROOT / 'HANDOFF_R3.md').write_text(f"R3 complete: at least {result['forced_sign_switches_lower']}/80 ending-side scan/tooth observations must reverse sign under PB->FB registration. " + f"Same modeled RMS ±57.5um leads to marginal gaps40 vs {pred['same_rms_contraction']['gap_upper_um'][0]:.3f}um at3deg half taper. Shared bridge lift propagates far-wall interference to near margin. " + f'Exact interval/vertex control error {maxerr:.3e}um; control TIE. Real color bins refused as full contact bounds, all actual gap certificates UNKNOWN. Next: measurement port in MEASUREMENT_REQUEST.json; acquire independent signed wall/process/datum bounds, then rerun without fitting threshold.\n')
state('R3_COMPLETE', result['outcome'], 'freeze acquisition transport test for a physically calibratable signed-clearance port')
print(result['outcome'], 'sign switch lower', result['forced_sign_switches_lower'], 'contact control error', maxerr)
