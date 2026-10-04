import json
from pathlib import Path
from common import ROOT, OUTPUT, sha, write

def check():
    if OUTPUT == ROOT:
        write('REPRODUCTION_RECEIPT.json', {'mode': 'original_compute', 'status': 'COMPUTED', 'replay': 'Run ./run_all.sh to compare scientific artifacts independently of timing'})
        return
    exact = ['raw/OBSERVATIONS_R2.json', 'raw/PATIENT_MODELS_R2.json', 'raw/DIAGNOSIS_CONCORDANCE_R2.json', 'raw/CONTROLS_R2.json', 'OBSERVATION_EDGES.json']
    rows = []
    for name in exact:
        same = sha((ROOT / name).read_bytes()) == sha((OUTPUT / name).read_bytes())
        rows.append({'path': name, 'source_sha256': sha((ROOT / name).read_bytes()), 'replay_sha256': sha((OUTPUT / name).read_bytes()), 'same': same})
    for (name, keys) in [('raw/RESULTS_R2.json', ['by_chain', 'by_quantity', 'paired_cases_by_chain', 'joint_localized_snapshot_rows', 'same_information_control_disagreements']), ('raw/RESULTS_R3.json', ['chains', 'certificate_failures', 'certificate_gate']), ('raw/RESULTS_R4.json', ['worlds', 'counterexample', 'real_shared_case', 'gate']), ('raw/VERIFICATION.json', ['status', 'sufficiency', 'source_scope_audit', 'injections', 'frozen_certificate_checks'])]:
        a = json.loads((ROOT / name).read_text())
        b = json.loads((OUTPUT / name).read_text())
        rows.append({'path': name, 'scientific_fields': keys, 'same': all((a[k] == b[k] for k in keys))})
    passed = all((x['same'] for x in rows))
    write('REPRODUCTION_RECEIPT.json', {'status': 'PASS' if passed else 'FAIL', 'comparisons': rows, 'timing_and_freeze_timestamps_expected_to_differ': True, 'originals_not_overwritten': True})
    print('REPRODUCTION', 'PASS' if passed else 'FAIL')
    if not passed:
        raise SystemExit(1)
if __name__ == '__main__':
    check()
