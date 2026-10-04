"""Replay actual controls and inject wrong values; expected scientific FAIL is preserved."""
import json, hashlib, zipfile, sys
from pathlib import Path
import numpy as np
from measure_r1 import HERE, sha, dump, tf2_io, exact_ray, side_measure
from measure_r2 import slab_control, edge
from measure_r3 import matrix_control
from motion_envelope import interval, query

def numeric_gate(candidate, reference, tolerance):
    return bool(np.max(abs(np.asarray(candidate) - np.asarray(reference))) <= tolerance)

def coverage_gate(teeth, scans, need_teeth, need_scans):
    return teeth >= need_teeth and scans >= need_scans

def edge_gate(discrepancy):
    return bool(np.quantile(discrepancy, 0.95) <= 0.3)

def survival_gate(survived, admitted):
    return admitted > 0 and survived / admitted >= 0.95

def biology_gate(sensitivity):
    return sensitivity >= 0.9

def main():
    tests = []
    for n in range(1, 6):
        f = json.load(open(HERE / f'FROZEN_PREDICTIONS_R{n}.json'))
        assert sha(HERE / f'PREREG_R{n}.json') == f['prereg_sha256']
        assert f['timestamp_utc']
    tests.append({'name': 'five_preregs_hash_frozen', 'status': 'PASS'})
    r = next((json.loads(l) for l in open(HERE / 'RAW_R1.jsonl') if json.loads(l).get('control_max_error_mm') is None and json.loads(l)['status'] == 'MEASURED_SECTION'))
    (lab, sp, _) = tf2_io.load(r['case'])
    sp = np.array(sp)
    o = np.array(r['origin_xyz_mm'])
    b = np.array(r['buccal_direction_xyz'])
    ex = exact_ray(lab, o, b, sp)
    sl = slab_control(lab, o, b, sp)
    exarr = np.array([x[:2] for x in ex])
    slarr = np.array([x[:2] for x in sl])
    assert numeric_gate(exarr, slarr, 1e-08)
    injected = exarr.copy()
    injected[len(injected) // 2, 1] += 1
    assert not numeric_gate(injected, slarr, 1e-08)
    tests.append({'name': 'actual_voxel_line_control_and_1mm_injection', 'reference_gate': 'PASS', 'injected_gate': 'FAIL', 'max_injected_error_mm': float(np.max(abs(injected - slarr)))})
    model = json.load(open(HERE / 'MOTION_MODELS.json'))[0]
    angles = np.array([-15.0, 0.0, 5.0, 15.0])
    (lo, hi) = interval(model, angles)
    (cl, ch) = matrix_control(model, angles)
    assert numeric_gate(np.r_[lo, hi], np.r_[cl, ch], 1e-09)
    inj = lo.copy()
    inj[2] += 1.0
    assert not numeric_gate(np.r_[inj, hi], np.r_[cl, ch], 1e-09)
    tests.append({'name': 'actual_matrix_control_and_1mm_injection', 'reference_gate': 'PASS', 'injected_gate': 'FAIL'})
    rows = [json.loads(l) for l in open(HERE / 'RAW_R2.jsonl')]
    dis = [r[n][e]['HM_gradient_discrepancy_mm'] for r in rows for n in ('buccal', 'lingual') for e in ('root_image_edge', 'bone_image_edge') if r[n].get(e, {}).get('status') == 'MEASURED_IMAGE_EDGE_SURROGATE']
    assert edge_gate(dis)
    assert not edge_gate(np.array(dis) + 1.0)
    tests.append({'name': 'all_identified_edge_control_and_1mm_injection', 'reference_gate': 'PASS', 'injected_gate': 'FAIL'})
    s2 = json.load(open(HERE / 'SUMMARY_R2.json'))
    assert coverage_gate(s2['bilateral_teeth_two_heights'], s2['covered_cases'], 20, 4)
    assert not coverage_gate(0, s2['covered_cases'], 20, 4)
    tests.append({'name': 'coverage_and_zero_teeth_injection', 'reference_gate': 'PASS', 'injected_gate': 'FAIL'})
    s4 = json.load(open(HERE / 'SUMMARY_R4.json'))
    assert not survival_gate(s4['admissions_surviving_sampled_R4'], s4['central_R3_admitted_queries'])
    assert survival_gate(s4['central_R3_admitted_queries'], s4['central_R3_admitted_queries'])
    assert not survival_gate(0, s4['central_R3_admitted_queries'])
    tests.append({'name': 'central_to_surface_survival_veto', 'actual_gate': 'FAIL', 'all_survive_positive_control': 'PASS', 'zero_survive_injection': 'FAIL'})
    src = json.load(open(HERE / 'SOURCES.json'))['external_referents'][0]
    assert not biology_gate(src['values']['dehiscence_sensitivity'])
    assert biology_gate(1.0)
    assert not biology_gate(0.1)
    tests.append({'name': 'external_static_model_biology_veto', 'actual_external_gate': 'FAIL', 'positive_control': 'PASS', 'wrong_low_sensitivity_injection': 'FAIL', 'locator': src['locator']})
    q = query(model, 0, 0, 100)
    assert q['conditional_geometric_status'] != 'WITHIN_SECTION_CONSTRAINTS'
    assert q['anatomical_status'] == 'UNDETERMINED'

    def admissible_port(q):
        return q['anatomical_status'] == 'UNDETERMINED' and q['biological_status'] == 'UNDETERMINED'
    assert admissible_port(q)
    bad = dict(q, anatomical_status='SAFE')
    assert not admissible_port(bad)
    tests.append({'name': 'unknown_physical_port_and_safe_status_injection', 'reference_gate': 'PASS', 'injected_gate': 'FAIL'})
    r5 = json.load(open(HERE / 'RAW_R5.json'))
    crown_actual = np.array([q['target_crown_advancement_mm'] + q['crown_target_error_mm'] for r in r5 for q in r['targets']])
    crown_target = np.array([q['target_crown_advancement_mm'] for r in r5 for q in r['targets']])
    assert numeric_gate(crown_actual, crown_target, 1e-09)
    assert not numeric_gate(crown_actual, crown_target + 1.0, 1e-09)
    tests.append({'name': 'inverse_crown_target_and_1mm_target_injection', 'reference_gate': 'PASS', 'injected_gate': 'FAIL'})
    for row in r5:
        for q5 in row['targets']:
            assert numeric_gate(q5['chosen_central_margin_mm'], q5['matrix_margin_mm'], 1e-09)
            assert not numeric_gate(q5['chosen_central_margin_mm'] + 1.0, q5['matrix_margin_mm'], 1e-09)
    tests.append({'name': 'inverse_root_margin_and_1mm_injection', 'reference_gate': 'PASS', 'injected_gate': 'FAIL'})
    manifest = json.load(open(HERE / 'DATA_MANIFEST_R1.json'))
    for f in manifest['files']:
        assert sha(f['path']) == f['sha256'], f['path']
    with zipfile.ZipFile(tf2_io.ZIP) as z:
        for p in manifest['provenance']:
            for m in p['members']:
                zi = z.getinfo(m['member'])
                assert zi.CRC == m['crc32'] and zi.file_size == m['bytes']
    tests.append({'name': 'raw_profiles_and_dataset_member_identity', 'status': 'PASS', 'files': len(manifest['files'])})
    probe = json.load(open(HERE / 'PROBE_EXPORT.json'))
    assert probe['non_two_incident_edges'] == 0
    assert abs(probe['mesh_signed_volume_mm3'] - probe['voxel_union_volume_mm3']) <= 1e-07
    for f in probe['files']:
        assert sha(f['path']) == f['sha256']
    tests.append({'name': 'closed_voxel_mesh_volume_and_hash', 'status': 'PASS'})
    reference = HERE / 'REFERENCE_NUMERICS.json'
    if reference.exists():
        saved = json.load(open(reference))
        for (file, h) in saved['raw_sha256'].items():
            assert sha(HERE / file) == h, file
        tests.append({'name': 'full_replay_raw_hash_agreement', 'status': 'PASS', 'files': len(saved['raw_sha256'])})
    report = {'status': 'PASS_REPRODUCIBILITY_AND_FALSIFIABLE_GATES', 'scientific_acceptance': 'PENDING_INDEPENDENT_REVIEW', 'clinical_or_anatomical_accuracy': 'UNKNOWN', 'tests': tests}
    dump(HERE / 'VERIFICATION.json', report)
    print(json.dumps(report, indent=2))
if __name__ == '__main__':
    main()
