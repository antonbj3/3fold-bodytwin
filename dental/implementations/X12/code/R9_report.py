from dental_release.paths import expand as _release_expand
import csv, hashlib, json, os, time
from collections import Counter
from pathlib import Path
import numpy as np
from consumer_ports import material_table, heat_table, compatible_anatomy
PACKAGE = Path(__file__).resolve().parents[1]
ROOT = Path(os.environ.get('X12_RUN_ROOT', str(PACKAGE)))
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X12'))

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def csvwrite(p, rows):
    if not rows:
        return
    with p.open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

def readrows(path):
    rows = list(csv.DictReader(path.open()))
    floats = ['pulp_voxels', 'tooth_voxels', 'containment', 'voxel_mm', 'pca_tilt_deg', 'crown_area_contrast', 'horn_boundary_min_mm', 'horn_boundary_p05_mm', 'horn_boundary_median_mm', 'occlusal_vertical_min_mm', 'axial_cardinal_horn_min_mm', 'axial_cardinal_1p5_below_min_mm', 'digital_two_surface_radius_mm', 'horn_boundary_lower_mm']
    for r in rows:
        r['fdi'] = int(r['fdi'])
        r['source_pulpy_fdi'] = int(r['source_pulpy_fdi'])
        for k in floats:
            r[k] = float(r[k]) if r.get(k) not in ('', None) else None
        r['domain_truncated'] = r['domain_truncated'] == 'True'
        r['calibration_anchor'] = r.get('calibration_anchor') == 'True'
    return rows

def summary(rows, key):
    out = []
    for group in sorted(set((r[key] for r in rows))):
        rr = [r for r in rows if r[key] == group]
        d = np.array([r['horn_boundary_min_mm'] for r in rr])
        o = np.array([r['occlusal_vertical_min_mm'] for r in rr if r['occlusal_vertical_min_mm'] is not None])
        a = np.array([r['axial_cardinal_horn_min_mm'] for r in rr if r['axial_cardinal_horn_min_mm'] is not None])
        out.append({key: group, 'n_teeth': len(rr), 'n_cases': len(set((r['case'] for r in rr))), 'hard_tissue_min_q05_mm': float(np.quantile(d, 0.05)), 'hard_tissue_min_median_mm': float(np.median(d)), 'hard_tissue_min_q95_mm': float(np.quantile(d, 0.95)), 'occlusal_projection_median_mm': float(np.median(o)) if len(o) else None, 'axial_cardinal_median_mm': float(np.median(a)) if len(a) else None, 'n_calibration_anchor': sum((r['calibration_anchor'] for r in rr))})
    return out

def main():
    t0 = time.time()
    rows = readrows(ROOT / 'raw/R8_teeth.csv')
    clean = [r for r in rows if not r['domain_truncated']]
    heldout = [r for r in clean if not r['calibration_anchor']]
    groups = summary(clean, 'tooth_type')
    fdi = summary(clean, 'fdi')
    csvwrite(ROOT / 'population_by_type.csv', groups)
    csvwrite(ROOT / 'population_by_FDI.csv', fdi)
    cal = json.loads((ROOT / 'raw/R8_calibration.json').read_text())
    r6 = json.loads((ROOT / 'raw/R6_pairs.json').read_text())
    c6 = json.loads((ROOT / 'raw/R6_cost.json').read_text())
    c8 = json.loads((ROOT / 'raw/R8_cost.json').read_text())
    cards = json.loads((PACKAGE / 'sources/product_cards.json').read_text())['cards']
    table = material_table(clean, cards)
    mat = DATA / ('material_uniform_budget.csv' if ROOT == PACKAGE else 'REPLAY_material_uniform_budget.csv')
    csvwrite(mat, table)
    htable = heat_table(clean)
    heat = DATA / ('heat_distance_scenarios.csv' if ROOT == PACKAGE else 'REPLAY_heat_distance_scenarios.csv')
    csvwrite(heat, htable)
    decisions = []
    for prod in [c['id'] for c in cards]:
        for x in [0, 0.25, 0.5, 1]:
            for delta in [0, 0.15, 0.3]:
                rr = [r for r in table if r['product'] == prod and r['buffer_mm'] == x and (r['annotation_delta_each_mm'] == delta)]
                decisions.append({'product': prod, 'buffer_mm': x, 'annotation_delta_each_mm': delta, 'n_teeth': len(rr), 'excluded_uniform_fraction': sum((r['decision'] == 'EXCLUDED_UNIFORM_HORN_OFFSET' for r in rr)) / max(1, len(rr))})
    csvwrite(ROOT / 'material_summary.csv', decisions)
    refs = json.loads((PACKAGE / 'sources/anatomical_referents.json').read_text())['referents']
    comparisons = []
    m1 = [r for r in clean if r['tooth_type'] == 'molar1']
    observed = float(np.median([r['occlusal_vertical_min_mm'] for r in m1])) if m1 else None
    for ref in refs:
        comparisons.append({**ref, 'atlas_quantity': 'topmost horn to outer tooth surface in same CT-axis column', 'atlas_molar1_vertical_median_mm': observed, 'nominal_difference_mm': observed - ref['mean_mm'] if observed is not None else None, 'comparison_status': 'UNKNOWN_UNMATCHED_LANDMARK_OR_TOOTH_TYPE', 'used_as_calibration': False})
    (ROOT / 'raw/ANATOMICAL_COMPARISON.json').write_text(json.dumps(comparisons, indent=2) + '\n')
    held = [r for r in r6 if int(r['case'][1:]) >= 300 and r.get('paired')]
    fraction = sum((r.get('aligned', False) for r in held)) / max(1, len(held))
    all_candidates = [r for r in r6 if int(r['case'][1:]) >= 300 and r.get('reason') != 'TRUNCATED_OR_MISSING_CANDIDATE']
    ct_fraction = sum((r.get('paired', False) for r in all_candidates)) / max(1, len(all_candidates))
    predictions = {'frozen_file': str(PACKAGE / 'FROZEN_PREDICTIONS.json'), 'sha256': sha(PACKAGE / 'FROZEN_PREDICTIONS.json'), 'heldout_case_ids_min': 300, 'heldout_CT_match_fraction': ct_fraction, 'CT_95pct_gate': ct_fraction >= 0.95, 'heldout_global_swap_fraction': fraction, 'global_swap_80pct_gate': fraction >= 0.8, 'n_exact_paired_heldout': len(held), 'digital_EDT_prediction': 'VERIFICATION_ALREADY_PERFORMED_BEFORE_FREEZE; not prospective physical prediction'}
    (ROOT / 'raw/FROZEN_PREDICTIONS_COMPARISON.json').write_text(json.dumps(predictions, indent=2) + '\n')
    valid = [r for r in cal if r.get('accepted')]
    mapcount = Counter((r['calibration']['selected_map'] for r in valid))
    control = json.loads((ROOT / 'raw/CONTROLS.json').read_text())
    thermal = json.loads((ROOT / 'raw/THERMAL_CONTROL.json').read_text())
    geom = json.loads((ROOT / 'raw/DEMO_GEOMETRY.json').read_text())
    manifest = [{'path': str(mat), 'sha256': sha(mat), 'bytes': mat.stat().st_size}, {'path': str(heat), 'sha256': sha(heat), 'bytes': heat.stat().st_size}] + json.loads((ROOT / 'raw/R8_manifest.json').read_text()) + geom['files']
    preregs = [{'path': str(p), 'sha256': sha(p)} for p in sorted(PACKAGE.glob('PREREG_*.json'))]
    code = [{'path': str(p), 'sha256': sha(p)} for p in sorted((PACKAGE / 'code').glob('*.py'))]
    gate = len(valid) >= 20 and sum((r['tooth_type'].startswith('molar') for r in clean)) >= 20
    result = {'lane': 'X12-pulpy3d', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'outcome': 'CALIBRATED_ANNOTATED_SAMPLE_ATLAS; CLINICAL_AND_LANDMARK_VALIDITY_UNKNOWN', 'population_gate_pass': gate, 'counts': {'visible_Pulpy_CT_records': 423, 'complete_CT_records': 422, 'R6_considered': len(r6), 'exact_paired_CT_cases': sum((r.get('paired', False) for r in r6)), 'R6_fixed_swap_aligned': sum((r.get('aligned', False) for r in r6)), 'R8_semantically_calibrated_cases': len(valid), 'R8_map_counts': dict(mapcount), 'measured_tooth_rows': len(rows), 'nontruncated_tooth_rows': len(clean), 'heldout_nontruncated_tooth_rows': len(heldout), 'nontruncated_cases': len(set((r['case'] for r in clean)))}, 'population_by_type': groups, 'population_by_FDI': fdi, 'material_decision_summary': decisions, 'prediction_comparison': predictions, 'external_referent': {'kind': 'published_dataset', 'locator': 'https://ditto.ing.unimore.it/toothfairy2/', 'compared_quantity': '95% FDI containment of non-anchor pulp voxels in independently provided whole-tooth labels, on exact same CT after two-anchor semantic calibration', 'refutes_us': not predictions['global_swap_80pct_gate'], 'scope': 'Refutes universal side encoding; sampled label congruence is not independent histological distance validation'}, 'external_anatomical_comparisons': comparisons, 'thermal_control': thermal, 'geometry_demo': geom, 'controls': control, 'equally_informed_control': {'method': 'independent KDTree against EDT on same boundaries', 'max_abs_error_mm': max((r['max_abs_error_mm'] for r in json.loads((ROOT / 'raw/R8_controls.json').read_text()))), 'outcome': 'TIE', 'methodological_novelty_claim': False}, 'cost': {'R6_full_CT_validation': c6, 'R8_semantic_and_mask_measurement': c8, 'consumer_wall_s': time.time() - t0, 'geometry_and_heat_wall_s': geom['wall_s'] + thermal['wall_s'], 'fit': '2discreteanchors percase; no geometric or continuous fit', 'preparation_and_failed_discovery_cost': 'Pilot timing logs retained; total agent preparation time and token use unmeasured', 'physical_measurement_and_fabrication': 'NOT_PERFORMED', 'fallback': 'rejected case rows and unverified boundaries remain UNKNOWN'}, 'uncertainty': {'digital_two_surface_envelope_mm': float(np.sqrt(3) * 0.3), 'annotation_delta_each_sensitivity_mm': [0, 0.15, 0.3], 'physical_accuracy': 'UNKNOWN', 'regional_orientation': 'area proxy, cusp and CEJ unvalidated', 'enamel_dentin_split': 'UNKNOWN', 'clinical_buffer': 'free scenario x; no recommended value'}, 'validity_limits': ['Annotated shared CT cohort, not random representative population', 'Two anchors are calibration, remaining labels held out; same-scan observations correlated', '150manual seed scans plus model-assisted remainder in published paper; per-case origin unknown', 'Total hard tissue, not dentin-only RDT', 'Top2mm pulp horn region proxy; per-cusp horns and axial dental planes unknown', 'No regional full-crown feasibility, no K2 global impossibility', 'No physically calibrated heat amplitude or clinical temperature'], 'preregistrations': preregs, 'code_manifest': code, 'data_manifest': manifest}
    (ROOT / 'results.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result['counts']))
    print('population_gate', gate, 'R6prospectiveglobalmap', predictions['global_swap_80pct_gate'])
if __name__ == '__main__':
    main()
