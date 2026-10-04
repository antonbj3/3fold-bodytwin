"""Acquire the frozen small test split only. No archive copies or training data."""
from dental_release.paths import expand as _release_expand
import json, pathlib, urllib.request, urllib.parse, time, hashlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = pathlib.Path(_release_expand('@DENTAL_WORK_ROOT@/LANE_XBREAK_HUNT_5'))
DATA.mkdir(parents=True, exist_ok=True)
tree = json.loads((ROOT / 'sources/PREP_TREE.json').read_text())
rows = []
for entry in tree['tree']:
    name = entry['path']
    if not (name.startswith('test Data/') and name.endswith('.vtp')):
        continue
    target = DATA / name.split('/')[-1]
    url = 'https://raw.githubusercontent.com/intellident-ai/teethPreparationData/' + tree['sha'] + '/' + urllib.parse.quote(name)
    start = time.perf_counter()
    if not target.exists():
        content = urllib.request.urlopen(url, timeout=40).read()
        if len(content) != entry['size']:
            raise ValueError('source byte mismatch ' + name)
        if not content.startswith(b'<?xml'):
            raise ValueError('not VTP ' + name)
        target.write_bytes(content)
    content = target.read_bytes()
    if len(content) != entry['size']:
        raise ValueError('existing byte mismatch ' + name)
    rows.append({'name': name, 'file': str(target), 'bytes': len(content), 'sha256': hashlib.sha256(content).hexdigest(), 'git_blob': entry['sha'], 'elapsed_seconds': time.perf_counter() - start, 'url': url})
(ROOT / 'sources/INPUTS.json').write_text(json.dumps({'commit': tree['sha'], 'test_files': rows, 'claimed_test_n': 13, 'actual_test_n': len(rows), 'missing_as_reported': 13 - len(rows), 'unit': 'mm (source dental geometry convention; scanner calibration unavailable)'}, indent=2))
print(json.dumps({'files': len(rows), 'bytes': sum((r['bytes'] for r in rows)), 'acquisition_s': sum((r['elapsed_seconds'] for r in rows))}))
