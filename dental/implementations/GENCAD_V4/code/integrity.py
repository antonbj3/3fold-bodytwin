"""Release integrity against an explicit external trust anchor (TRUSTED_RELEASE.sha256)."""
from pathlib import Path
import hashlib, json, os
ROOT = Path(__file__).resolve().parents[1]

class IntegrityError(RuntimeError):
    pass

def digest_file(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda : f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def lock():
    try:
        expected = (ROOT / 'TRUSTED_RELEASE.sha256').read_text().split()[0]
        blob = (ROOT / 'RELEASE_LOCK.json').read_bytes()
    except OSError as e:
        raise IntegrityError('MISSING_RELEASE_ANCHOR') from e
    if hashlib.sha256(blob).hexdigest() != expected:
        raise IntegrityError('RELEASE_MANIFEST_HASH_MISMATCH')
    return json.loads(blob)

def verify():
    m = lock()
    for (rel, item) in m['files'].items():
        p = ROOT / rel
        if not p.is_file() or p.stat().st_size != item['bytes'] or digest_file(p) != item['sha256']:
            raise IntegrityError('INPUT_HASH_MISMATCH ' + rel)
    for directory in m['closed_directories']:
        base = ROOT / directory
        actual = set()
        for p in base.rglob('*'):
            if p.is_symlink():
                raise IntegrityError('INTERNAL_SYMLINK ' + str(p))
            if p.is_file():
                actual.add(str(p.relative_to(ROOT)))
        expected = {r for r in m['files'] if r.startswith(directory + '/')}
        if actual != expected:
            raise IntegrityError('CLOSED_DIRECTORY_CONTENTS ' + directory)
    return dict(files=len(m['files']), sha256=digest_file(ROOT / 'RELEASE_LOCK.json'), bytes=sum((x['bytes'] for x in m['files'].values())))

def checked_json_bytes(path):
    p = Path(path)
    blob = p.read_bytes()
    rel = str(p.absolute().relative_to(ROOT))
    m = lock()
    item = m['files'].get(rel)
    if item and hashlib.sha256(blob).hexdigest() != item['sha256']:
        raise IntegrityError('READ_AFTER_VERIFY_MISMATCH ' + rel)
    return blob
if __name__ == '__main__':
    print(json.dumps(verify()))
