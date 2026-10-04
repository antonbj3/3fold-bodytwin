"""Verify delivered hashes after a reporting-only change; no numerical rerun."""
import subprocess, sys
from shared import *
latest = load(R / 'LATEST_RUN.json')
run = Path(latest['path'])
cmd = [sys.executable, str(R / 'code/package.py'), str(run)]
with (run / 'package_scope_audit.log').open('a') as f:
    subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=True)
revision = run / 'REPORTING_REVISION.json'
rec = load(revision)
rec.update(finalized_utc=now(), script=artifact(R / 'code/package.py'), audit_script=artifact(Path(__file__)), current_results=artifact(R / 'results.json'), freeze_index=artifact(R / 'FROZEN_PREDICTIONS_INDEX.json'))
dump(revision, rec)
latest['reporting_revision'] = artifact(revision)
latest['results'] = artifact(R / 'results.json')
dump(R / 'LATEST_RUN.json', latest)
for p in sorted(R.glob('PREREG_*.json')):
    h = p.with_suffix('.sha256')
    if not h.exists():
        h = Path(str(p) + '.sha256')
    assert h.read_text().strip() == sha(p), 'Prereg hash mismatch'
assert load(R / 'GRAPH_FEEDBACK.json')['result_sha256'] == sha(R / 'results.json')
for p in R.glob('patients/*/DECISION_CONTRACT.json'):
    x = load(p)
    assert x['release'] == 'ABSTAIN' and x['cross_patient_edges'] == []
    print(x['patient_id'], 'words', len((p.parent / 'PATIENT360.md').read_text().split()))
print('Frozen entries', len(load(R / 'FROZEN_PREDICTIONS_INDEX.json')['files']))
print('Frozen and feedback hashes verified; reports and releases consistent')
state('FOUR_ROUNDS_COMPLETE', 'R4 exact bone PASS / sampled bone FAIL; 27 numerical controls and injections PASS; crown design FAIL; full scope PARTIAL; clinical/manufacturing ABSTAIN', 'Independent review; same-patient IOS+CBCT registration and measured preparation, loaded bite/force or matched pulp thermometry')
current = load(R / 'CURRENT_WORK_STATE.json')
current.update(latest_run=str(run), reporting_revision=artifact(revision))
dump(R / 'CURRENT_WORK_STATE.json', current)
