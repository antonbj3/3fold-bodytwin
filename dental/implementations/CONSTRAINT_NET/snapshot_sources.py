"""Freeze small contracts only; never copy source trees or arrays."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
import datetime, hashlib, json, subprocess
HERE = Path(__file__).resolve().parent
DENT = HERE.parents[1]
R = DENT / 'results'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def save(p, obj):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n')

def main():
    if (HERE / 'SOURCE_MANIFEST.json').exists():
        raise SystemExit('Already frozen')
    files = [DENT / 'notes/NORTHSTAR_SPEC.md', DENT / 'notes/OLD_DENTAL_CELL_INDEX.md', R / 'GRAPH_WORKING_VIEW_20260923/WORKING_VIEW.json', R / 'DEMO48_PACKAGE/demos.json', R / 'DEMO48_PACKAGE/FACIT.json']
    reviews = []
    for p in sorted(R.glob('LANE_XREVIEW*/REVIEW_*.json')):
        d = json.loads(p.read_text())
        files.append(p)
        decision = d.get('graph_decision', d.get('decision', ''))
        q = d.get('result_file')
        if q:
            q = Path(q)
            q = q if q.is_absolute() else DENT / q
            if q.exists() and q.stat().st_size < 5000000:
                files.append(q)
        reviews.append({'review_file': str(p), 'decision': decision, 'result_file': str(q) if q else None, 'declared_sha256': d.get('result_sha256'), 'actual_sha256': sha(q) if q and q.exists() else None})
    for p in (R / 'DEMO48_PACKAGE').glob('demos/*/results.json'):
        files.append(p)
    for p in (R / 'DEMO48_PACKAGE').glob('batch*/demos/*/results.json'):
        if p.stat().st_size < 5000000:
            files.append(p)
    for name in ['CURRENT_WORK_STATE.json', 'INPUT_LOCK.json', 'inputs/DECISIONS.json', 'inputs/MEASUREMENTS.json', 'inputs/LAB_ASSUMPTIONS.json', 'inputs/EXTERNAL_FACIT.json']:
        p = R / _release_expand('X71') / name
        if p.exists():
            files.append(p)
    manifest = []
    for p in sorted(set(files)):
        h = sha(p)
        dest = HERE / 'sources' / f'{h}{p.suffix}'
        dest.parent.mkdir(exist_ok=True)
        if not dest.exists():
            dest.write_bytes(p.read_bytes())
        manifest.append({'source_file': str(p), 'sha256': h, 'bytes': p.stat().st_size, 'snapshot': str(dest.relative_to(HERE))})
    save(HERE / 'SOURCE_MANIFEST.json', {'frozen_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'files': manifest})
    save(HERE / 'REVIEW_CENSUS.json', reviews)
    packet = subprocess.run(['./graph', 'working', 'packet', '--id', 'DENT-UNC-PROPAGATION'], cwd=DENT, capture_output=True, text=True)
    (HERE / 'GRAPH_PACKET.json').write_text(packet.stdout)
    (HERE / 'GRAPH_PACKET.stderr').write_text(packet.stderr)
    print(json.dumps({'snapshot_files': len(manifest), 'bytes': sum((x['bytes'] for x in manifest)), 'packet_exit': packet.returncode}))
if __name__ == '__main__':
    main()
