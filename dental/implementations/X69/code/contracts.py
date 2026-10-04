"""Reuse Patient360 validator and expose provenance that its schema cannot verify."""
from dental_release.paths import expand as _release_expand
import copy, sys, time, types, json
from registration import *
LEGACY = pathlib.Path(_release_expand('@DENTAL_IMPLEMENTATIONS@/PATIENT360_R2/code'))

def strict_bridge(bridge, a, b):
    if bridge.get('metric_unit_evidence') != 'SOURCE_MM_CALIBRATED':
        raise ValueError('UNKNOWN_STL_METRIC_SCALE')
    if bridge.get('heldout_landmarks', {}).get('provenance_kind') != 'INDEPENDENT_ANNOTATED_LANDMARKS':
        raise ValueError('INDEPENDENT_LANDMARK_PROVENANCE_MISSING')
    if bridge.get('heldout_landmarks', {}).get('reference_locator') == bridge.get('fit_reference_locator'):
        raise ValueError('LANDMARK_REFERENCE_REUSED_FOR_FIT')
    return legacy_validate(bridge, a, b)

def legacy_validate(bridge, a, b):
    if str(LEGACY) not in sys.path:
        sys.path.insert(0, str(LEGACY))
    from region_field import validate_bridge
    return validate_bridge(bridge, a, b)

def make_bridge(jaw, r, sha_ios, sha_cbct):
    ident = {'patient_id': _release_expand('@DENTAL_CASE_ID@'), 'source_sha256': sha_ios, 'target_sha256': sha_cbct, 'pairing_status': 'SOURCE_DECLARED_SAME_SUBJECT', 'locator': 'https://doi.org/10.5281/zenodo.8027553', 'scope': 'Same subject per publisher Demo_1 folder; reconstruction provenance from primary paper. No independent pose or anatomical trueness.'}
    path = ROOT / 'raw' / f'IDENTITY_{jaw}.json'
    dump(path, ident)
    w = np.load(DATA / 'r2_witnesses.npz')
    p = w[jaw + '_source_points']
    d = w[jaw + '_residual_mm']
    q = w[jaw + '_closest_cbct_surface']
    fdi = w[jaw + '_fdi']
    idx = np.flatnonzero(np.isin(fdi, [11, 21, 31, 41]))[::100]
    b = {'schema': 'Patient360-IOS-CBCT-registration-v2', 'patient_id': ident['patient_id'], 'unit': 'mm', 'source_frame_id': 'IOS_NATIVE_' + jaw, 'target_frame_id': 'CBCT_DEMO1', 'source_source_sha256': sha_ios, 'target_source_sha256': sha_cbct, 'source_to_target': r['source_to_target'], 'identity_evidence': {'locator': ident['locator'], 'path': str(path), 'sha256': sha(path)}, 'fit_landmark_ids': ['FIT_CROWN_OUTER_TEETH'], 'heldout_landmarks': {'ids': ['NEAREST_SURFACE_' + str(int(i)) for i in idx], 'source_xyz_mm': p[idx].tolist(), 'target_xyz_mm': q[idx].tolist(), 'provenance_kind': 'NEAREST_SURFACE_WITNESS_NOT_INDEPENDENT_LANDMARK', 'reference_locator': 'Same CBCT-only mesh used by ICP'}, 'heldout_acceptance_mm': 0.5, 'uniform_registration_error_bound_mm': None, 'time_scale': 'SIMULTANEOUS', 'acquisition_state': 'Same-subject acquisitions; exact timing and bite-state relation UNKNOWN; independent jaw coordinate link only', 'evidence_kind': 'published_dataset', 'metric_unit_evidence': 'ASSUMED_MM_FROM_PUBLISHER_STL_PIPELINE', 'fit_reference_locator': 'Same CBCT-only mesh used by ICP', 'warning': 'This is a candidate bridge; nearest-surface witnesses are deliberately not called independent landmarks.'}
    a = types.SimpleNamespace(meta={'patient_id': b['patient_id'], 'unit': 'mm', 'frame_id': b['source_frame_id'], 'source_sha256': sha_ios})
    target = types.SimpleNamespace(meta={'patient_id': b['patient_id'], 'unit': 'mm', 'frame_id': b['target_frame_id'], 'source_sha256': sha_cbct})
    dump(ROOT / 'raw' / f'PATIENT360_CANDIDATE_BRIDGE_{jaw}.json', b)
    return (b, a, target)

def run_contracts():
    t = time.monotonic()
    r = json.loads((ROOT / 'raw/R2_RESULTS.json').read_text())
    files = json.loads((ROOT / 'raw/HAO_ACQUISITION_MANIFEST.json').read_text())['files']
    cbct = next((x['sha256'] for x in files if x['member'].endswith('out_smoothed.stl')))
    rows = []
    for jawresult in r['jaws']:
        jaw = jawresult['jaw']
        bundle = {'files': jawresult['source_labels']}
        path = ROOT / 'raw' / f'IOS_SOURCE_BUNDLE_{jaw}.json'
        dump(path, bundle)
        (b, a, target) = make_bridge(jaw, jawresult, sha(path), cbct)
        for (name, fn) in [('original_surface_only_bridge', legacy_validate), ('strict_provenance_bridge', strict_bridge)]:
            try:
                value = fn(b, a, target)
                rows.append({'jaw': jaw, 'check': name, 'status': 'ACCEPTED', 'value': value})
            except ValueError as e:
                rows.append({'jaw': jaw, 'check': name, 'status': 'ABSTAIN', 'reason': str(e)})
        scale_declared = copy.deepcopy(b)
        scale_declared['metric_unit_evidence'] = 'SOURCE_MM_CALIBRATED'
        try:
            strict_bridge(scale_declared, a, target)
            rows.append({'jaw': jaw, 'check': 'strict_landmark_gate_with_scale_declared', 'status': 'FAULT_ACCEPTED'})
        except ValueError as e:
            rows.append({'jaw': jaw, 'check': 'strict_landmark_gate_with_scale_declared', 'status': 'ABSTAIN', 'reason': str(e)})
        good = copy.deepcopy(b)
        hold = good['heldout_landmarks']
        hold['target_xyz_mm'] = transform(np.array(hold['source_xyz_mm']), np.array(b['source_to_target'])).tolist()
        good['heldout_landmarks']['provenance_kind'] = 'INDEPENDENT_ANNOTATED_LANDMARKS'
        original = legacy_validate(good, a, target)
        rows.append({'jaw': jaw, 'check': 'injected_fit_derived_zero_residual_false_landmarks', 'legacy_status': 'ACCEPTED', 'producer_verdict': 'CONTROL_FAIL_INDEPENDENCE_UNCHECKED', 'evidence_kind': 'our_own_fixture', 'scope': 'Original validator cannot establish independence from ids alone', 'value': original})
        for (name, mut) in [('wrong_patient', lambda x: x.update(patient_id='different')), ('wrong_unit', lambda x: x.update(unit='cm')), ('wrong_hash', lambda x: x.update(source_source_sha256='0' * 64)), ('wrong_frame', lambda x: x.update(source_frame_id='other')), ('translation_error', lambda x: x['source_to_target'][0].__setitem__(3, x['source_to_target'][0][3] + 10)), ('scale_error', lambda x: x['source_to_target'][0].__setitem__(0, x['source_to_target'][0][0] * 1.1 + 0.01)), ('fit_holdout_overlap', lambda x: x.update(fit_landmark_ids=x['heldout_landmarks']['ids'])), ('unsupported_uniform_bound', lambda x: x.update(uniform_registration_error_bound_mm=0.01))]:
            bad = copy.deepcopy(good)
            mut(bad)
            try:
                legacy_validate(bad, a, target)
                ok = False
                reason = 'FAULT_ACCEPTED'
            except ValueError as e:
                ok = True
                reason = str(e)
            rows.append({'jaw': jaw, 'check': name, 'fault_rejected': ok, 'reason': reason, 'evidence_kind': 'our_own_fixture'})
    support = np.array([[-1.0, -1.0, 0.0], [-1.0, 1.0, 0.0], [1.0, -1.0, 0.0], [1.0, 1.0, 0.0]])
    A = np.eye(4)
    B = np.eye(4)
    B[:3, :3] = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
    dA = cKDTree(support).query(transform(support, A))[0]
    dB = cKDTree(support).query(transform(support, B))[0]
    root = np.array([[1.0, 0.0, -10.0]])
    s = {'summary': 'sorted residual vector, RMS, p95, coverage', 'summary_identity_error_mm': float(np.max(np.abs(np.sort(dA) - np.sort(dB)))), 'machine_exact_identical': bool(np.array_equal(dA, dB)), 'summary_A_mm': dA.tolist(), 'summary_B_mm': dB.tolist(), 'downstream_root_A_mm': transform(root, A).tolist(), 'downstream_root_B_mm': transform(root, B).tolist(), 'downstream_difference_mm': float(np.linalg.norm(transform(root, A) - transform(root, B))), 'resolution': 'PER_POINT', 'minimal_extension': 'Directed anatomically named asymmetric correspondences plus separate jaw state; root correspondence is not identified by a surface residual histogram.', 'external_referent': {'kind': 'our_own_fixture', 'locator': str(ROOT / 'code/contracts.py'), 'compared_quantity': 'exact residual identity and distinct transformed root coordinate; not external measurement truth', 'refutes_us': True}, 'evidence_kind': 'our_own_fixture', 'physical_scope': 'NONE; mathematical insufficiency witness'}
    dump(ROOT / 'raw/SUFFICIENCY.json', s)
    dump(ROOT / 'raw/CONTRACT_RESULTS.json', {'rows': rows, 'sufficiency': s, 'seconds': time.monotonic() - t, 'original_validator_path': str(LEGACY / 'region_field.py'), 'original_validator_sha256': sha(LEGACY / 'region_field.py'), 'clinical_release': 'ABSTAIN_INDEPENDENT_LANDMARKS_AND_METRIC_CALIBRATION_MISSING', 'source_metadata_validation': 'Same-subject source declaration is accepted separately from registration trueness'})
    print('Contract fault checks', len(rows), 'root displacement', s['downstream_difference_mm'])
if __name__ == '__main__':
    run_contracts()
