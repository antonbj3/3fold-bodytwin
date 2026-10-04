from pathlib import Path
import hashlib, json, datetime
ROOT = Path(__file__).resolve().parents[1]

def hash_file(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda : f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()
files = [{'file': str(p.relative_to(ROOT)), 'bytes': p.stat().st_size, 'sha256': hash_file(p)} for p in sorted(ROOT.rglob('*')) if p.is_file() and '__pycache__' not in str(p) and (p.name not in ['ARTIFACT_MANIFEST.json', 'CURRENT_WORK_STATE.json'])]
(ROOT / 'ARTIFACT_MANIFEST.json').write_text(json.dumps({'sealed_at': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'files': files}, indent=2) + '\n')
s = json.loads((ROOT / 'CURRENT_WORK_STATE.json').read_text())
s.update(replay_launcher='PASS ./run_all.sh; raw/ONE_COMMAND_RUN.log; exit0', status='ROUND_COMPLETE_PENDING_REVIEW', updated_at=datetime.datetime.now(datetime.timezone.utc).isoformat(), ongoing='Five rounds complete and one-command replay verified', latest_gate='Relative tracking available; absolute movement and physical force-response UNKNOWN;97/108 force-stage numerical cases', next_operation='Independent annotated repeat-reference scan and matched seated wrench/time-history measurement;then freeze heldout prediction and validate', handoff='HANDOFF.md', review_state='PENDING_INDEPENDENT_REVIEW', artifact_manifest_sha256=hash_file(ROOT / 'ARTIFACT_MANIFEST.json'))
(ROOT / 'CURRENT_WORK_STATE.json').write_text(json.dumps(s, indent=2) + '\n')
