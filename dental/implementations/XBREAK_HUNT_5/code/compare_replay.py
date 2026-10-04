from common import *
import sys
replay = pathlib.Path(sys.argv[1])
checks = []

def compare(a, b, path):
    if isinstance(a, dict):
        for k in a:
            if k in {'wall_seconds', 'elapsed_seconds', 'stl', 'stl_sha256'}:
                continue
            compare(a[k], b[k], path + '/' + k)
    elif isinstance(a, list):
        if len(a) != len(b):
            raise ValueError('list length mismatch ' + path)
        for (i, (x, y)) in enumerate(zip(a, b)):
            compare(x, y, path + '/' + str(i))
    elif isinstance(a, float):
        if abs(a - b) > 1e-08:
            raise ValueError('numeric replay mismatch ' + path)
    elif a != b:
        raise ValueError('replay mismatch ' + path)
for name in ['R1', 'R2', 'R3', 'R4']:
    original = json.loads((ROOT / f'raw/{name}.json').read_text())
    fresh = json.loads((replay / f'raw/{name}.json').read_text())
    compare(original, fresh, name)
    checks.append(name)
out = {'gate': 'PASS', 'stages': checks, 'numeric_tolerance': 1e-08, 'ignored_fields': ['wall_seconds', 'elapsed_seconds', 'stl output path', 'STL hash (separately check closed export/witnesses)'], 'replay': str(replay), 'physical_gate': 'UNKNOWN'}
(replay / 'REPLAY_COMPARISON.json').write_text(json.dumps(out, indent=2))
(ROOT / 'raw/LATEST_REPLAY.json').write_text(json.dumps(out, indent=2))
print(json.dumps(out))
