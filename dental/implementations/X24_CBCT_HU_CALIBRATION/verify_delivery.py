from pathlib import Path
import json, hashlib, sys
import numpy as np
from run_r2 import weighted_median
from calibration import bone_branch_b
P = Path(__file__).resolve().parent

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def verify(outdir=None):
    outdir = Path(outdir or P)
    checks = {}
    for name in ['PREREG_R1', 'PREREG_R2', 'PREREG_R3', 'FROZEN_PREDICTIONS', 'FROZEN_RESPONSE_R2', 'FROZEN_PREDICTIONS_R3']:
        checks[name + '_hash'] = sha(P / (name + '.json')) == (P / (name + '.sha256')).read_text().strip()
    if (P / 'INPUT_MANIFEST.json').exists():
        for row in json.loads((P / 'INPUT_MANIFEST.json').read_text()):
            checks['input:' + row['path']] = sha(row['path']) == row['sha256']
    bas = json.loads((P / 'PHANTOM_BASELINE.json').read_text())
    new = json.loads((P / 'PHANTOM_NEW_R2.json').read_text())
    err = 0.0
    n = 0
    for d in bas + new:
        rr = d['raw_samples']
        checks['raw:' + Path(rr['path']).name] = sha(rr['path']) == rr['sha256']
        with np.load(rr['path']) as a:
            for row in d['inserts']:
                m = row['material']
                pred = weighted_median(a[m + '_gray'], a[m + '_legacy_multiplicity'])
                err = max(err, abs(pred - row['gray_median']))
                n += 1
    checks['raw_median_replay_132_ROIs'] = err <= 1e-12 and n == 132
    with np.load(bas[0]['raw_samples']['path']) as a:
        v = a['PMP_gray']
        w = a['PMP_legacy_multiplicity'].copy()
        cleanmedian = weighted_median(v, w)
        w[np.argmax(v)] += 100 * w.sum()
        badmedian = weighted_median(v, w)
    checks['corrupt_multiplicity_rejected'] = abs(badmedian - cleanmedian) > 1
    r1 = json.loads((outdir / 'results_R1.json').read_text())
    r2 = json.loads((outdir / 'results_R2.json').read_text())
    r3 = json.loads((outdir / 'results_R3.json').read_text())
    checks['R1_negative_preserved'] = r1['outcome'] == 'NEGATIVE_TWO_ANCHOR_SUFFICIENCY' and (not r1['gates']['holdout_40HU'])
    checks['R2_negative_preserved'] = r2['outcome'] == 'NEGATIVE_THREE_ANCHOR_SUFFICIENCY' and (not r2['gates']['heldout40HU'])
    checks['R2_100HU_injection_failure_preserved'] = not r2['gates']['error_injection']
    checks['R3_controls_pass'] = all(r3['gates'].values())
    checks['all_variant_1000HU_errors_rejected'] = all((abs(r['error_HU'] + 1000) > 40 for r in r2['rows']))
    checks['K42_no_false_acceptance'] = all((r['status'] == 'ABSTAIN_UNKNOWN' for r in r3['K42_certificates']))
    checks['joint_domain_conservation'] = r3['dropout']['source_sites'] == r3['dropout']['retained_sites'] + r3['dropout']['rejected_sites']
    budgetchecks = []
    for r in r3['answers']:
        H = r['H_reference_closure']
        dh = r['HU_symmetric_budget_5percent_E']
        e = float(bone_branch_b(H)['E_MPa'])
        lo = float(bone_branch_b(H - dh)['E_MPa'])
        hi = float(bone_branch_b(H + dh)['E_MPa'])
        budgetchecks.append(max(abs(lo / e - 1), abs(hi / e - 1)) <= 0.05 + 1e-12)
    checks['five_percent_E_budget_endpoints'] = all(budgetchecks)
    checks['external_archive_md5_identity'] = json.loads((P / 'PHANTOM_ARCHIVE_IDENTITY.json').read_text())['identity_match']
    links = json.loads((P / 'EXTERNAL_ARCHIVE_MEMBER_LINKS.json').read_text())
    checks['external_archive_baseline_and_ROI_links'] = links['status'] == 'PASS' and len(links['members']) == 28
    for row in links['members']:
        checks['archive_link_current:' + row['path']] = row['matches_original_archive'] and sha(row['path']) == row['sha256']
    if (outdir / 'CONTROL_FAULT_AUDIT.json').exists():
        audit = json.loads((outdir / 'CONTROL_FAULT_AUDIT.json').read_text())
        checks['control_fault_audit'] = audit['status'] == 'PASS' and all((r['fault_rejected'] for r in audit['controls']))
    if outdir != P:
        for name in ['results_R1', 'results_R2', 'results_R3']:
            ref = json.loads((P / (name + '.json')).read_text())
            got = json.loads((outdir / (name + '.json')).read_text())
            ref.pop('cost', None)
            got.pop('cost', None)
            checks['rerun_science_identical:' + name] = ref == got
    status = {'ok': all(checks.values()), 'checks': checks, 'raw_ROIs_replayed': n, 'raw_median_max_error_gray': err, 'multiplicity_fault_shift_gray': badmedian - cleanmedian}
    (outdir / 'DELIVERY_CHECKS.json').write_text(json.dumps(status, indent=2) + '\n')
    print('PASS' if status['ok'] else 'FAIL', len(checks), 'checks')
    if not status['ok']:
        print({k: v for (k, v) in checks.items() if not v})
        raise SystemExit(1)
    return status
if __name__ == '__main__':
    verify(sys.argv[1] if len(sys.argv) > 1 else None)
