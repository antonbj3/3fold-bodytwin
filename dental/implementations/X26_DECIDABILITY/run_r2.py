"""Frozen-source dental table and typed reachability checks."""
import copy
import csv
import hashlib
import json
import math
import time
from collections import Counter
from pathlib import Path
from decidability import ScalarPort, interval_decision, refinement_sequence
ROOT = Path(__file__).resolve().parent

def read_json(name):
    return json.loads((ROOT / 'inputs' / name).read_text())

def read_csv(name):
    with (ROOT / 'inputs' / name).open() as f:
        return list(csv.DictReader(f))

def save(name, obj):
    (ROOT / name).write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n')

def hash_ok():
    for x in json.loads((ROOT / 'SOURCE_MANIFEST.json').read_text()):
        if hashlib.sha256((ROOT / x['local']).read_bytes()).hexdigest() != x['sha256']:
            return False
    for name in ['PREREG_R1', 'PREREG_R2', 'FROZEN_PREDICTIONS']:
        if hashlib.sha256((ROOT / (name + '.json')).read_bytes()).hexdigest() != (ROOT / (name + '.sha256')).read_text().strip():
            return False
    return True

def table_row(decision, quantity, level, sigma, klass, binding, certification='', **kwargs):
    return dict(decision=decision, quantity=quantity, resolution_level=level, dominant_sigma_or_bound=sigma, **{'class': klass}, binding_quantity=binding, certifying_resolution_or_budget=certification, physical_status='UNKNOWN', claim_type='capability', timescale='SIMULTANEOUS', **kwargs)

def validate(raw, parents, edge, zero, noseq, signed, frozen_hashes):
    return {'frozen_hashes': frozen_hashes, 'arithmetic_endpoint_match': all(((x['conditional_status'] != 'ABSTAIN') == x['direct_endpoint_certifies'] and abs(x['envelope'] - x['floor'] - x['numerical_error_bound']) < 1e-10 for x in raw)), 'finite_boundary_strict': edge['class'] == 2 and edge['pitch_limit_strict'] is None and (zero['class'] == 2), 'physical_refusals': all((x['physical_status'] == 'UNKNOWN' for x in raw)) and all((r['physical_status'] == 'UNKNOWN' for r in parents)), 'eight_parent_questions': len(parents) == 8, 'no_dental_divergence_inferred': noseq['class'] is None, 'signed_fixture_refutes_scalar_shortcut': signed['expansion_um'] < 120 < signed['contraction_um']}

def main():
    start = time.perf_counter()
    assert hash_ok(), 'frozen input/prereg hash mismatch'
    canal = read_json('GEOM_field.json')['canal']
    geom = read_json('GEOM_field.json')
    bone = geom['bone_outer']
    root = geom['root']
    bias = abs(canal['all']['bias'])
    sigmas = (canal['all']['sigma_w'], canal['tau'])
    floor = bias + 3 * sum(sigmas)
    input_rows = read_csv('X5_margins.csv')
    rows = []
    rejected = Counter()
    raw = []
    for (n, r) in enumerate(input_rows):
        if r['region'] != 'canal_body':
            rejected['outside_canal_question'] += 1
            continue
        try:
            m = float(r['signed_margin_mm'])
            if not math.isfinite(m):
                raise ValueError()
        except (KeyError, ValueError):
            rejected['missing_or_nonfinite_margin'] += 1
            continue
        port = ScalarPort(m, 'mm', bias, sigmas)
        q = port.evaluate(0.3)
        lo = m - bias - 3 * sum(sigmas) - math.sqrt(3) * 0.3
        hi = m + bias + 3 * sum(sigmas) + math.sqrt(3) * 0.3
        direct = lo > 0 or hi < 0
        q.update(record_index=n, case=r['case'], site=r['site'], selection=r['selection'], direct_endpoint_certifies=direct, conditional_transfer='POPULATION field -> regional minimum is an unvalidated closure; no surface confidence')
        raw.append(q)
        h = q['suggested_sufficient_pitch']
        rows.append(table_row(f"X5:{r['case']}:{r['site']}:{r['selection']}", 'canal_body clearance minus 2 mm engineering guard', 'PER_SURFACE_REGION', f'sigma_w={sigmas[0]:.12g} mm; tau={sigmas[1]:.12g} mm; bias_bound={bias:.12g} mm', q['class'], 'same-site image-to-anatomy and simultaneous surface error; conditional image-label floor', f"h < {q['pitch_limit_strict']:.12g} mm; suggested h={h:.12g} mm" if h else 'no finite resampling pitch in conditional model', certifying_scope='conditional scalar envelope ONLY; not anatomy or full-surface confidence', numerical_pitch_mm=0.3, conditional_floor_mm=floor, signed_margin_mm=m, conditional_current_status=q['conditional_status'], field_source_level='POPULATION', conditional_floor_resolution_level='PHENOMENOLOGICAL', numerical_pitch_resolution_level='PER_POINT', signed_margin_resolution_level='PER_SURFACE_REGION', external_locator='https://ditto.ing.unimore.it/toothfairy2/; inputs/X5_margins.csv', confidence='3 sigma conditional Gaussian scalar; surface coverage UNKNOWN'))
    x18 = read_json('X18_results.json')['force_uncertainty_capability']['rows']
    force = []
    for r in x18:
        lower_witness = r['independent_mixed_solve']['peak_MPa']
        upper = r['sharp_upper_peak_MPa']
        answer = 'ROBUST_BELOW' if upper < 100 else 'ABOVE_WITNESS' if lower_witness > 100 else 'UNRESOLVED_FORCE_DISTRIBUTION'
        force.append(dict(fdi=r['fdi'], patches=r['patch_count'], upper_MPa=upper, lower_witness_MPa=lower_witness, engineering_threshold_MPa=100, response=answer, ratio=r['upper_over_area_peak']))
    x10 = read_json('X10_R4.json')
    x13 = read_json('X13_R4.json')
    x1b = read_json('X1b_results.json')
    x24 = read_json('X24_results_R1.json')
    parents = [table_row('X5:physical_local_canal', 'true implant-to-canal minimum clearance', 'PER_SURFACE_REGION', f'image-label sigma {sigmas[0]:.6g} mm; anatomy error UNKNOWN', 2, 'matched local anatomy and surface-tail coverage', 'paired anatomical reference; resampling fixed 0.3mm acquisition is insufficient', class_evidence='unvalidated local/surface transfer'), table_row('X10:physical_gap', 'projected vertical marginal-plane separation', 'PER_SURFACE_REGION', 'signed registered scanner error UNKNOWN; group medians/RMS are not local bounds', 2, 'signed registered normal error and contact lift law', 'same-specimen signed field + seated gap measurement', class_evidence='independent signed study rejects universal scalar summary'), table_row('X18:physical_crown_stress', 'maximum tensile stress at known total force', 'PER_POINT', 'force response upper/area ratio ' + f"{min((r['upper_over_area_peak'] for r in x18)):.4g}..{max((r['upper_over_area_peak'] for r in x18)):.4g}; not sigma", 2, 'patch force fractions, contact registration, support/cement law', 'registered force map; mesh refinement sequence also required', class_evidence='same geometry and total force yield unequal stresses'), table_row('X1b:absolute_fracture_load', 'single-crown fracture load', 'PER_TOOTH', f"held-study logRMSE={x1b['literature']['R2_log_ratio_RMSE']:.9g}; discrepancy, not SD", 2, 'cement/support, flaw population and fracture origin matched to specimen', 'matched cement/defect/fracture-origin/load experiment', timescale_override='HANDOVER', class_evidence='held-study calibration fails; no specimen law'), table_row('X15:biological_envelope', 'dehiscence after a specified tooth trajectory', 'PER_SURFACE_REGION', f"bone/root image-label sigma_w={bone['all']['sigma_w']:.6g}/{root['all']['sigma_w']:.6g} mm; local anatomical/biological error UNKNOWN", 2, 'anatomical bone edge plus trajectory-to-remodeling response', 'same-site anatomical baseline and later bone outcome', timescale_override='HANDOVER', class_evidence='static bone/moved tooth sensitivity 0.53 in independent study; not our sensitivity'), table_row('X13:cement_inversion', 'CAD setting to regional assembled gap and hydraulic conductance', 'PER_SURFACE_REGION', f"R1 region errors 71.0/64.1/213.4 um; same-mean flow capacity ratio {x13['aliasing_conductance_ratio']:.9g}; not SD", 2, 'within-lab regional calibration and spatial cement film', 'two settings plus held-out batch; spatial film/viscosity for flow', class_evidence='external group means reject identity spacer; own topology witness is only analytic'), table_row('X20:physical_short_vs_sinus', 'full-body bone/sinus clearance and clinical choice', 'PER_SURFACE_REGION', 'height +/-0.6mm scenario box; physical sigma UNKNOWN', 2, 'local pose and bone/sinus truth; outcome response separate', 'fiducial-matched whole-body geometry; no clinical certificate', class_evidence='failed reread and unknown anatomy; cohort is not local truth'), table_row('X24:bone_stiffness', 'CBCT gray to mandibular modulus', 'PER_SURFACE_REGION', f"held phantom max error {x24['max_heldout_error_HU']:.9g} HU_ref; not modulus sigma", 2, 'scanner attenuation calibration plus bone density/fabric/modulus law', 'matched HA phantom + mandibular mechanical test', class_evidence='two-anchor heldout gate fails; polymer facit not bone modulus')]
    for r in parents:
        if 'timescale_override' in r:
            r['timescale'] = r.pop('timescale_override')
    rows = parents + rows
    fields = list(dict.fromkeys((k for row in rows for k in row)))
    with (ROOT / 'DECIDABILITY_TABLE.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    save('RAW_R2_CANAL.json', raw)
    save('R2_FORCE_SCREENS.json', force)
    signed_case = {'expansion_um': x10['case_bounds']['same_rms_expansion']['gap_upper_um'][0], 'contraction_um': x10['case_bounds']['same_rms_contraction']['gap_lower_um'][0], 'threshold_um': 120, 'same_information_unmeasured_sign': True, 'external_signed_observation': {'lower_switch_count': x10['forced_sign_switches_lower'], 'eligible': x10['eligible_ending_scan_tooth_observations']}, 'fixture_scope': 'own analytic contact fixture, not independent physical gap truth'}
    assert signed_case['expansion_um'] < 120 < signed_case['contraction_um']
    save('R2_SIGNED_WITNESS.json', signed_case)
    edge = ScalarPort(3.0, 'mm', 0, (1.0,), coefficient=1).evaluate(0.1)
    zero = ScalarPort(0, 'mm').evaluate(0.1)
    noseq = refinement_sequence([10.0], [0.3])
    tests = validate(raw, parents, edge, zero, noseq, signed_case, hash_ok())
    kwargs = dict(raw=raw, parents=parents, edge=edge, zero=zero, noseq=noseq, signed=signed_case, frozen_hashes=True)
    mutation = {}
    changes = {'frozen_hashes': {'frozen_hashes': False}, 'arithmetic_endpoint_match': {'raw': [dict(raw[0], envelope=raw[0]['envelope'] + 1)] + raw[1:]}, 'finite_boundary_strict': {'edge': dict(edge, **{'class': 1, 'pitch_limit_strict': 0.1})}, 'physical_refusals': {'parents': [dict(parents[0], physical_status='CERTIFIED')] + parents[1:]}, 'eight_parent_questions': {'parents': parents[:-1]}, 'no_dental_divergence_inferred': {'noseq': dict(noseq, **{'class': 3})}, 'signed_fixture_refutes_scalar_shortcut': {'signed': dict(signed_case, contraction_um=signed_case['expansion_um'])}}
    for (gate, change) in changes.items():
        mutation[gate] = not validate(**kwargs | change)[gate]
    mutation['legacy_naive_count'] = json.loads((ROOT / 'legacy/reports/decidability_abstention_atlas.json').read_text())['gates']['G3_naive_misfire_count'] != 2
    counts = Counter((x['class'] for x in raw))
    current = Counter((x['conditional_status'] for x in raw))
    result = {'round': 'R2', 'claim_type': 'capability', 'outcome': 'CONDITIONAL_RESOLUTION_BUDGETS_AND_PHYSICAL_BINDERS', 'gates': tests, 'mutation_rejections': mutation, 'canal_floor': {'bias_bound_mm': bias, 'sigma_w_mm': sigmas[0], 'tau_mm': sigmas[1], 'worst_covariance_sigma_envelope_mm': sum(sigmas), 'conditional_floor_mm': floor, 'resolution_level': 'PHENOMENOLOGICAL', 'replacement_measurement': 'matched local signed anatomy error; full-surface maximum distribution', 'point_confidence': 0.9973002039367398, 'physical_surface_confidence': 'UNKNOWN'}, 'canal_counts': dict(counts), 'canal_current_status': dict(current), 'parent_physical_class2_count': len(parents), 'dropout': {'input_margin_rows': len(input_rows), 'canal_kept': len(raw), 'rejected': dict(rejected), 'fraction_outside_scope_or_missing': sum(rejected.values()) / len(input_rows)}, 'force_screens': force, 'same_information_control': 'direct endpoints match, a consistency check; no algorithm superiority claim', 'divergence': 'UNKNOWN for actual dental meshes; no multiresolution sequence in X18', 'cost': {'wall_seconds': time.perf_counter() - start, 'fit_seconds': 0, 'source_preparation_and_human_work': 'UNKNOWN', 'physical_measurement': 'NOT_PERFORMED', 'independent_review': 'PENDING'}, 'external_referent': {'kind': 'published_dataset', 'locator': 'https://ditto.ing.unimore.it/toothfairy2/; inputs/X5_margins.csv and SOURCE_MANIFEST.json', 'compared_quantity': 'label-based per-region canal clearance margins; not physical local max error', 'refutes_us': False}, 'next_operation': 'R3 measured-force constraint polytope replaces simplex; compare to external convexity theorem and independent source FEM mixed solve'}
    save('RESULTS_R2.json', result)
    assert all(tests.values()) and all(mutation.values())
    print(json.dumps({'canal_counts': dict(counts), 'current': dict(current), 'floor_mm': floor, 'gates': tests, 'mutations': mutation}))
if __name__ == '__main__':
    main()
