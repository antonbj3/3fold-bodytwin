"""Hash lane-owned scientific data without copying any source tree."""
from dental_release.paths import expand as _release_expand
import hashlib, json, datetime
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X42_gen_beat_standard'))

def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda : f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def run():
    files = []
    for p in sorted(DATA.rglob('*')):
        if not p.is_file() or 'replay_evaluation' in p.relative_to(DATA).parts:
            continue
        files.append({'path': str(p.relative_to(DATA)), 'bytes': p.stat().st_size, 'sha256': sha(p)})
    out = {'scope': 'Lane-owned scientific arrays, models, frozen predictions and first trusted evaluation; replays separately recorded', 'data_root': str(DATA), 'files': files, 'total_bytes': sum((f['bytes'] for f in files)), 'max_intermediate_bytes': 3000000000, 'within_budget': sum((f['bytes'] for f in files)) < 3000000000}
    dest = ROOT / 'DATA_MANIFEST.json'
    dest.write_text(json.dumps(out, indent=2) + '\n')
    print({'file_count': len(files), 'total_bytes': out['total_bytes'], 'within_budget': out['within_budget'], 'manifest_sha256': sha(dest)})
if __name__ == '__main__':
    run()
