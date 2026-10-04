import sys, json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gencad_bench.io import sha, DATA
for prereg in ROOT.glob('PREREG_*.json'):
    if sha(prereg) != prereg.with_suffix('.sha256').read_text().strip():
        raise ValueError('Preregistration drift: ' + str(prereg))
for row in json.loads((ROOT / 'data/DEPENDENCIES_FROZEN.json').read_text())['files']:
    path = ROOT / row['path'][6:] if row['path'].startswith('@lane/') else Path(row['path'])
    if sha(path) != row['sha256']:
        raise ValueError('Inherited dependency drift: ' + str(path))
for row in json.loads((ROOT / '_vendor/MANIFEST.json').read_text())['files']:
    if sha(ROOT / row['path']) != row['sha256']:
        raise ValueError('Local vendor drift: ' + row['path'])
size = sum((p.stat().st_size for base in [ROOT, DATA] for p in base.rglob('*') if p.is_file() and (not p.is_symlink())))
if size >= 2000000000:
    raise RuntimeError('Lane intermediate size limit reached')
print('preflight: inherited hashes match; lane bytes', size, flush=True)
