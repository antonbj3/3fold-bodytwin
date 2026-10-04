"""Always hash current bytes; never trust a previous PASS receipt."""
from common import *
import time, sys

def verify_sources():
    t = __import__('time').monotonic()
    checks = []
    for r in read(ROOT / 'inputs/PARENT_LOCK.json'):
        assert sha(r['path']) == r['sha256'], f"Parent drift: {r['path']}"
    for r in read(ROOT / 'inputs/AUXILIARY_LOCK.json'):
        assert sha(r['path']) == r['sha256'], f"Auxiliary source drift: {r['path']}"
    for n in [1, 2, 3, 4]:
        p = ROOT / f'PREREG_R{n}.json'
        if p.exists():
            assert sha(p) == read(ROOT / f'PREREG_R{n}.sha256.json')['sha256'], f'Prereg drift R{n}'
    expected = {line.split()[1]: line.split()[0] for line in (DATA / 'SHA256SUMS').read_text().splitlines()}
    for (name, digest) in expected.items():
        actual = sha(DATA / name)
        assert actual == digest, f'ZIP bytes mismatch {name}'
        checks.append(dict(path=str(DATA / name), sha256=actual, bytes=(DATA / name).stat().st_size, match=True))
    names = {Path(r['cache']).name for r in read(X37 / 'raw/SOURCE_SURFACE_HASHES.json')}
    used = []
    for r in read(X37 / 'DATA_MANIFEST.json'):
        if Path(r['path']).name in names:
            assert sha(r['path']) == r['sha256'], f"Cloud drift {r['path']}"
            used.append(r)
    assert len(used) == 1040
    write(ROOT / 'raw/USED_CLOUD_MANIFEST.json', used)
    out = dict(current_bytes_verified=True, zips=checks, clouds_verified=len(used), cloud_bytes=sum((r['bytes'] for r in used)), wall_s=time.monotonic() - t, stale_receipt_accepted=False)
    write(ROOT / 'raw/INTEGRITY.json', out)
    print('Current source/cache bytes verified', len(used), round(out['wall_s'], 2), flush=True)
    return out
if __name__ == '__main__':
    verify_sources()
