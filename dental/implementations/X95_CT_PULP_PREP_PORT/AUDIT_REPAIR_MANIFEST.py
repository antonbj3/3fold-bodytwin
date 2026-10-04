"""Repair the one stale final-note hash after verifying every other binding."""
from dental_release.paths import expand as _release_expand
import hashlib
import json
from pathlib import Path
import sys
ORIGINAL = Path(_release_expand('@DENTAL_IMPLEMENTATIONS@/X95_CT_PULP_PREP_PORT'))
FINAL_NOTE = 'SOL_NIGHT_FINAL_X95-ct-pulp-prep-port.md'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    package = Path(sys.argv[1]).resolve()
    if package == ORIGINAL.resolve():
        raise SystemExit('Refuse to modify original delivery')
    manifest = package / 'DELIVERY_MANIFEST.json'
    data = json.loads(manifest.read_text())
    bad = []
    for row in data['files']:
        path = Path(row['path'])
        path = path if path.is_absolute() else package / path
        if not path.is_file() or sha(path) != row['sha256']:
            bad.append(row['path'])
    if bad != [FINAL_NOTE]:
        raise SystemExit('Expected only the stale final note, found ' + repr(bad))
    row = next((r for r in data['files'] if r['path'] == FINAL_NOTE))
    before = row.copy()
    row.update(sha256=sha(package / FINAL_NOTE), bytes=(package / FINAL_NOTE).stat().st_size)
    manifest.write_text(json.dumps(data, indent=2) + '\n')
    print(json.dumps({'only_changed_binding': FINAL_NOTE, 'before': before, 'after': row, 'scientific_result_sha256_unchanged': sha(package / 'results.json')}))
if __name__ == '__main__':
    main()
