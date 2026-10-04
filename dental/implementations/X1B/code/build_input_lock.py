"""Once all local inputs and exports exist, freeze their portable manifest."""
from common import R, dump, sha
inputs = []
for folder in ['inputs', 'exports', 'vendor', 'history']:
    for p in sorted((R / folder).rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts:
            inputs.append(dict(path=str(p.relative_to(R)), sha256=sha(p), bytes=p.stat().st_size))
for p in sorted(R.glob('PREREG_*.json')):
    inputs.append(dict(path=p.name, sha256=sha(p), bytes=p.stat().st_size))
inputs.append(dict(path='RUNTIME_LOCK.json', sha256=sha(R / 'RUNTIME_LOCK.json'), bytes=(R / 'RUNTIME_LOCK.json').stat().st_size))
dump('INPUT_LOCK.json', dict(files=inputs, total_bytes=sum((p['bytes'] for p in inputs))))
(R / 'INPUT_LOCK.sha256').write_text(sha(R / 'INPUT_LOCK.json') + '\n')
print('Pinned', len(inputs), 'files', sum((p['bytes'] for p in inputs)) / 1000000.0, 'MB')
