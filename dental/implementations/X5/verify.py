import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '1'
from pathlib import Path
import argparse, hashlib, json
import numpy as np
from run_r1 import sha, put
from voxel_decisions import verify_patch_selection, witness_subset_report, verify_witness_subset_claim, local_cube_field, direct_cube_field, paired_evaluate
from decision_error import fixed_design_record
P = Path(__file__).resolve().parent
D = Path(__import__('os').environ.get('DENTAL_PROJECT_ROOT', str(P.parent.parent))).resolve()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--inputs-only', action='store_true')
    a = ap.parse_args()
    checks = {}
    for name in ['PREREG_R1', 'PREREG_R2', 'PREREG_R3', 'FROZEN_PREDICTIONS', 'FROZEN_PREDICTIONS_R3']:
        checks[name + '_hash'] = sha(P / (name + '.json')) == (P / (name + '.sha256')).read_text().strip()
    source = json.load(open(P / 'SOURCE_DESIGN_MANIFEST.json'))
    checks['all_369_archived_source_hashes'] = all((sha(r['path']) == r['sha256'] for r in source))
    if (P / 'SOURCE_SMALL_MANIFEST.json').exists():
        source = json.load(open(P / 'SOURCE_SMALL_MANIFEST.json'))
        checks['all_small_source_hashes'] = all((sha(r['path']) == r['sha256'] for r in source))
    for name in ['DATA_MANIFEST.json', 'DATA_MANIFEST_R3.json']:
        checks[name + '_hashes'] = all((sha(r['path']) == r['sha256'] for r in json.load(open(P / name))))
    if a.inputs_only:
        assert all(checks.values()), checks
        print(json.dumps(checks))
        return
    s1 = json.load(open(P / 'SUMMARY_R1.json'))
    s2 = json.load(open(P / 'SUMMARY_R2.json'))
    s3 = json.load(open(P / 'SUMMARY_R3.json'))
    required_computational_gates = [['inverse_replay', 'fault_injection'], ['source_replay', 'direct_allpairs', 'equal_metrics_decision_disagreement', 'patch_compression', 'direct_patch_replay', 'direction_crossing', 'wrong_direction_rejected', 'corrupted_gap_rejected', 'corrupted_field_rejected', 'omitted_dangerous_cell_rejected', 'forged_full_surface_witness_rejected', 'prediction_control'], ['minimum_connected_extension', 'global_metric_disagreement', 'postprocessing_survival', 'fault_injection']]
    checks['all_computational_gates_pass'] = all((all((s.get('gates', {}).get(k) is True for k in keys)) for (s, keys) in zip([s1, s2, s3], required_computational_gates)))
    rr = np.array([2.1, 3.4])
    good = np.array([True, False])
    bad = np.array([False, False])
    checks['omitted_dangerous_patch_rejected'] = verify_patch_selection(rr, 2.0, 0.3, good) and (not verify_patch_selection(rr, 2.0, 0.3, bad))
    report = witness_subset_report([3.2, 4.0], 2.0)
    checks['forged_far_claim_rejected'] = verify_witness_subset_claim([3.2, 4.0], 2.0, report) and (not verify_witness_subset_claim([3.2, 4.0], 2.0, {**report, 'can_certify_full_mask_above_threshold': True}))
    sentinel = fixed_design_record({'canal_body': None, 'lateral_bone': 2.0}, {'canal_body': 2.0, 'lateral_bone': 1.0})
    checks['missing_hazard_cannot_pass'] = sentinel['nominal_pass'] is None and sentinel['adverse_radius_mm'] is None
    X = np.array([[0, 0, 0], [1, 2, 1]])
    Y = np.array([[4, 1, 3], [2, 4, 6]])
    sp = np.array([0.15, 0.3, 0.6])
    f = local_cube_field(X, Y, sp)['distance_mm']
    q = direct_cube_field(X, Y, sp)
    checks['anisotropic_cube_control'] = bool(np.max(abs(f - q)) < 1e-12)
    try:
        local_cube_field(X, np.empty((0, 3), int), sp)
        empty_rejected = False
    except ValueError:
        empty_rejected = True
    checks['empty_hazard_rejected'] = empty_rejected
    big = np.stack(np.meshgrid(np.arange(-5, 6), np.arange(-5, 6), np.arange(-5, 6), indexing='ij'), -1).reshape(-1, 3)
    shell = big[np.max(abs(big), axis=1) == 5]
    contained = paired_evaluate(np.array([[0, 0, 0]]), big, shell, [0.3, 0.3, 0.3], 1.0)
    checks['solid_containment_rejected'] = not contained['reference']['nominal_pass'] and contained['prediction']['nominal_pass'] and contained['decision_flip']
    rows = list(map(json.loads, open(P / 'RAW_R1_SITE_GUARDS.jsonl')))
    ids = [(r['case'], r['site'], r.get('selection', 'unknown')) for r in rows]
    checks['unique_source_selection_rows'] = len(ids) == len(set(ids))
    checks['R1_physical_unknown_all'] = all((r.get('physical_status', 'UNKNOWN').startswith('UNKNOWN') for r in rows))
    checks['R2_cube_replay_max_below_gate'] = max((r['controls']['direct_allpairs_max_error_mm'] for r in json.load(open(P / 'RAW_R2_FIELDS.json')))) <= 1e-09
    checks['R3_true_minimality'] = all((r['controls']['minimality_exact'] and r['controls']['k_minus_one_claim_rejected'] for r in json.load(open(P / 'R3_CONNECTED_COUNTEREXAMPLES.json'))))
    out = {'passed': all(checks.values()), 'checks': checks, 'scope': 'Producer numerical/integrity verification; independent scientific review and matched physical accuracy not provided'}
    put('VERIFICATION.json', out)
    assert out['passed'], out
    print(json.dumps(out))
if __name__ == '__main__':
    main()
