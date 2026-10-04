"""Version operational bindings in our copied package; preserve all scientific hashes."""
import json, shutil
from pathlib import Path
from prepare_replay import HERE, PACKAGE, DATA, sha
OPERATIONAL = {'PACKAGE_ACTIVE_INDEX.json', 'demos.json', 'package_entry.py', 'run_demos.sh'}

def contract():
    old = json.loads((PACKAGE / 'PACKAGE_CODE_LOCK.json').read_text())
    current = json.loads((HERE / 'PACKAGE_CODE_LOCK_RELEASE_R7.json').read_text())
    assert sha(HERE / 'PACKAGE_CODE_LOCK_RELEASE_R7.json') == (HERE / 'PACKAGE_CODE_LOCK_RELEASE_R7.json.sha256').read_text().split()[0]
    a = {r['path']: r['sha256'] for r in old['files']}
    b = {r['path']: r['sha256'] for r in current['files']}
    assert a.keys() == b.keys()
    changed = {name for name in a if a[name] != b[name]}
    assert changed == OPERATIONAL and len(a) - len(changed) == 69
    return (current, b)

def install(clone):
    (current, expected) = contract()
    for name in sorted(OPERATIONAL):
        src = DATA / 'metadata_snapshot' / name
        assert sha(src) == expected[name], 'Current frozen metadata snapshot drift: ' + name
        dst = clone / name
        assert not dst.is_symlink()
        if dst.exists():
            dst.unlink()
        shutil.copy2(src, dst)
    dst = clone / 'PACKAGE_CODE_LOCK.json'
    assert not dst.is_symlink()
    if dst.exists():
        dst.unlink()
    shutil.copy2(HERE / 'PACKAGE_CODE_LOCK_RELEASE_R7.json', dst)
    return {'version': 'R7_FROZEN_BATCH25_OPERATIONAL_BINDING', 'new_lock_sha256': sha(dst), 'changed_operational_entries': 4, 'preserved_original_entries': 69, 'scientific_input_or_threshold_changes': 0}
