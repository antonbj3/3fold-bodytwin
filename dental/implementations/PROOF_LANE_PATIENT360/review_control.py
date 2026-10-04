"""Independent retrospective review of delivered artifacts; no new model or experiment."""
from dental_release.paths import expand as _release_expand
import argparse, copy, hashlib, json, struct
from pathlib import Path
import numpy as np

def load(p):
    return json.loads(Path(p).read_text())

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def refs(x):
    if isinstance(x, dict):
        if {'path', 'sha256'} <= x.keys():
            yield x
        for v in x.values():
            yield from refs(v)
    elif isinstance(x, list):
        for v in x:
            yield from refs(v)

def check_identity(frame, contract, frame_path=None):
    assert frame['patient_id'] == contract['patient_id'], 'FRAME_PATIENT_MISMATCH'
    assert frame['unit'] == 'mm', 'FRAME_UNIT_MISMATCH'
    assert sha(frame_path or contract['frame']['path']) == contract['frame']['sha256'], 'FRAME_HASH_MISMATCH'
    if contract['track'] == 'CBCT':
        assert frame['array_axes'] == 'zyx' and frame['canonical_axes'] == 'xyz'

def semantic_correct(rows, frame):
    for r in rows:
        expected = r['fdi'] if frame['semantic_map'] == 'identity' else r['fdi'] + 10 if r['fdi'] < 40 else r['fdi'] - 10
        if r['semantic_map'] != frame['semantic_map'] or r['source_pulpy_fdi'] != expected:
            return False
    return True

def audit(root, run=None):
    root = Path(root).resolve()
    x = load(root / 'results.json')
    run = Path(run or x['run_directory'])
    checks = []
    errors = []
    metrics = {}
    suff = []

    def record(name, passed, scope):
        checks.append({'check': name, 'pass': bool(passed), 'scope': scope})
        if not passed:
            errors.append(name)
    contracts = {}
    for p in sorted((root / 'patients').glob('*/DECISION_CONTRACT.json')):
        c = load(p)
        pid = c['patient_id']
        contracts[pid] = c
        f = load(run / pid / 'FRAME.json')
        try:
            check_identity(f, c, run / pid / 'FRAME.json')
            passed = True
        except AssertionError:
            passed = False
        record(pid + '_actual_frame_binding', passed, 'PER_POINT identity and declared units; scanner accuracy UNKNOWN')
        record(pid + '_release', c['release'] == 'ABSTAIN' and c['cross_patient_edges'] == [], 'No physical or cross-patient release')
        for ref in refs(c):
            pth = Path(ref['path'])
            record(pid + '_artifact_' + pth.name, pth.is_file() and sha(pth) == ref['sha256'], 'Delivered evidence bytes')
    ios = run / _release_expand('@DENTAL_CASE_C@')
    pairs = load(ios / 'X21_REPLAY.json')['pairs']
    report = load(ios / 'REPORT_COMPARISON.json')
    comparisons = [a for r in report['reports'] for a in r['comparison'].values()]
    metrics['tooth_pairs'] = len(pairs)
    metrics['report_covered'] = sum((v['covered'] for v in comparisons))
    metrics['report_total'] = len(comparisons)
    record('IOS_table_counts', metrics['tooth_pairs'] == x['IOS']['tooth_pairs'] and metrics['report_covered'] == x['IOS']['report_fields_covered'] and (metrics['report_total'] == x['IOS']['report_fields_total']), 'PER_SURFACE_REGION / PER_ARCH diagnostic only')
    cr = load(ios / 'CROWN_RESULTS.json')
    frozen = load(ios / 'FROZEN_CROWN_PREDICTIONS.json')
    record('FE_frozen_result_equality', cr['FE'] == frozen['FE'], 'Conditional numerical predictions')
    record('X34_actual_status', load(ios / 'DESIGN_GATE.json')['verdict'] == x['IOS']['X34'] == 'FAIL', 'Design failure retained')
    r3 = load(run / 'R3/RESULTS_R3.json')
    exp = load(root / 'data' / run.name / 'R3/export/FROZEN_PREDICTIONS.json')
    rf = load(ios / 'FRAME.json')
    record('actual_export_patient_frame', exp['patient_id'] == contracts[_release_expand('@DENTAL_CASE_C@')]['patient_id'] and exp['frame'] == rf, 'Exported patient geometry frame')
    record('actual_export_theta', exp['theta_rows'] == r3['theta_rows'], 'Full downstream scenario rows')
    record('FE_force_theta', all((r['FE']['status'] == 'UNKNOWN_NO_CONTACT' or abs(r['FE']['total_load_N'] - r['shared_theta']['force_N']) < 1e-10 for r in r3['theta_rows'])), 'PHENOMENOLOGICAL force')
    delta = cr['FE']['100']['maximum_tensile_MPa'] - cr['FE']['50']['maximum_tensile_MPa']
    roof = np.load(root / 'data' / run.name / _release_expand('@DENTAL_CASE_C@/roof.npz'))
    area = cr['generated_contact_area_mm2']
    summary_a = np.array([area], dtype=np.float64)
    summary_b = np.array([area], dtype=np.float64)
    suff.append({'patient_id': 'synthetic', 'summary': 'generated contact area and same roof/contact mask', 'state_locators': [str(ios / 'FROZEN_CROWN_PREDICTIONS.json') + '#/FE/50', str(ios / 'FROZEN_CROWN_PREDICTIONS.json') + '#/FE/100'], 'same_geometry_sha256': sha(root / 'data' / run.name / _release_expand('@DENTAL_CASE_C@/roof.npz')), 'summary_identity_error': float(abs(summary_a - summary_b).max()), 'bitwise_identical': bool(np.array_equal(summary_a.view(np.uint64), summary_b.view(np.uint64))), 'downstream': 'maximum tensile MPa', 'downstream_difference': delta, 'decision': 'INSUFFICIENT', 'smallest_extension_in_this_fixed_normalized_load_family': 'total force N; full spatial load field still required beyond this fixed family', 'scope': 'Existing 50/100N scenarios, no new FE solve; not a claim that arbitrary geometry compression is sufficient'})
    bone = []
    pulp = []
    for b in x['CBCT']:
        pid = b['patient_id']
        f = load(run / pid / 'FRAME.json')
        m = load(run / pid / 'X12_MEASUREMENTS.json')
        cb = load(run / pid / 'CBCT_RESULTS.json')
        record(pid + '_semantic_lineage', semantic_correct(m['rows'], f), 'Per-tooth source identity')
        record(pid + '_pulp_count', len(m['rows']) == b['pulp_teeth'], 'PER_TOOTH')
        br = load(run / 'R4' / (pid + '_EXACT_BONE.json'))
        n = sum((v[k]['primary_exact']['status'] == 'MEASURED_ANNOTATION_SURROGATE' for v in br['rows'] for k in ['buccal', 'lingual']))
        record(pid + '_bone_count', n == b['bone_sides'], 'PER_SURFACE_REGION')
        record(pid + '_sinus_abstention', len(cb['sinus']['eligible_sites']) == 0 and len(cb['sinus']['rejected']) == 8 and all((r['reason'] == 'PRESENT_TOOTH_NOT_EDENTULOUS' for r in cb['sinus']['rejected'])), 'Matched missing-site condition')
        q = np.load(root / 'data' / run.name / pid / '36_paired.npz')
        whole = np.load(root / 'data' / run.name / pid / 'labels.npz')['labels']
        coords = np.argwhere(q['tooth'])
        sgn = f['crops']['36']['crown_axis_sign']
        if sgn < 0:
            coords[:, 0] = q['tooth'].shape[0] - 1 - coords[:, 0]
        coords += q['crop_origin_zyx']
        record(pid + '_crop_source_identity', np.all(whole[tuple(coords.T)] == 36), 'All retained crop tooth voxels, source label identity')
        xyz = coords[:, ::-1] * q['spacing'][::-1] @ np.array(f['direction']).T + f['offset_xyz_mm']
        back = (xyz - f['offset_xyz_mm']) @ np.linalg.inv(np.array(f['direction'])).T / q['spacing'][::-1]
        err = float(abs(back[:, ::-1] - coords).max() * q['spacing'].max())
        record(pid + '_all_crop_world_roundtrip', err < 1e-08, 'PER_POINT transform, no registration accuracy claim')
        th = load(run / 'R2' / pid / 'THERMAL_RESULTS.json')
        hist = np.loadtxt(run / 'R2' / pid / 'raw/history_h0.csv', delimiter=',', skiprows=1)
        peak = float(hist[:, 1].max())
        record(pid + '_thermal_raw_peak', abs(peak - th['scenarios'][0]['peak_pulp_delta_C']) < 1e-12, 'Saved raw time samples, PHENOMENOLOGICAL')
        states = th['scenarios']
        a = float(states[0]['input_J'])
        z = float(states[-1]['input_J'])
        diff = states[0]['peak_pulp_delta_C'] - states[-1]['peak_pulp_delta_C']
        suff.append({'patient_id': pid, 'summary': 'input heat J, same crop, same source path and retained fraction', 'state_locators': [str(run / 'R2' / pid / 'THERMAL_RESULTS.json') + '#/scenarios/0', str(run / 'R2' / pid / 'THERMAL_RESULTS.json') + '#/scenarios/2'], 'same_geometry_sha256': th['crop']['sha256'], 'summary_values': [a, z], 'summary_identity_error': abs(a - z), 'bitwise_identical': struct.pack('d', a) == struct.pack('d', z), 'downstream': 'peak pulp temperature rise degC', 'downstream_difference': diff, 'decision': 'INSUFFICIENT', 'smallest_extension_in_this_fixed_cooling_family': 'cooling coefficient h with inherited boundary temperature and boundary face mask; outside this family retain full cooling boundary field', 'scope': 'Existing h=0 and h=2000 states only; geometric-summary sufficiency with h held fixed remains UNTESTED'})
        bone.append(n)
        pulp.append(len(m['rows']))
        metrics[pid] = {'pulp_count': len(m['rows']), 'bone_sides': n, 'crop_roundtrip_mm': err, 'raw_thermal_peak_C': peak, 'canal_bracket_mm': cb['canal_clearance_mm']}
    record('producer_controls', x['controls']['controls'] == 27 and all((c['valid_pass'] and c['fault_rejected'] for c in x['controls']['checks'])), 'Specified injected faults only')
    return {'status': 'PASS' if not errors else 'FAIL', 'checks': checks, 'failed_checks': errors, 'metrics': metrics, 'sufficiency': suff, 'rigorous_FE_enclosure': 'MISSING; two scalar load replays show numerical linearity only; fixed linear model scaling is not a certified floating-point enclosure', 'claim_type': 'capability'}
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('root', type=Path)
    ap.add_argument('--run')
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    r = audit(a.root, a.run)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(r, indent=2, ensure_ascii=False) + '\n')
    print(r['status'], r['failed_checks'])
    raise SystemExit(0 if r['status'] == 'PASS' else 1)
