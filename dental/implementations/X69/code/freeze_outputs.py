"""Freeze every numerical run's own outputs, without copying acquired source trees."""
import datetime, shutil
from registration import *

def main():
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    dest = DATA / 'runs' / stamp
    local = ROOT / 'runs' / stamp
    own = ['r1_geometry.npz', 'r1_witnesses.npz', 'r2_witnesses.npz', 'figure_r1.npz', 'figure_r2.npz']
    paths = [DATA / name for name in own] + list((DATA / 'exports').glob('*'))
    total = sum((p.stat().st_size for p in DATA.rglob('*') if p.is_file()))
    copybytes = sum((p.stat().st_size for p in paths))
    if total + copybytes > 3000000000:
        raise RuntimeError('lane immutable snapshot cap')
    if shutil.disk_usage(DATA).free - copybytes < 5000000000:
        raise RuntimeError('disk reserve')
    dest.mkdir(parents=True)
    local.mkdir(parents=True)
    mapping = {}
    for p in paths:
        q = dest / ('exports' if p.parent.name == 'exports' else '') / p.name
        q.parent.mkdir(exist_ok=True)
        shutil.copy2(p, q)
        if sha(p) != sha(q):
            raise ValueError('snapshot hash mismatch')
        mapping[str(p)] = str(q)
    for file in ['raw/R1_RESULTS.json', 'raw/R2_RESULTS.json', 'EXPORT_SIDECAR.json']:
        p = ROOT / file
        s = p.read_text()
        for (a, b) in mapping.items():
            s = s.replace(a, b)
        p.write_text(s)
    (local / 'raw').mkdir()
    for p in (ROOT / 'raw').glob('*.json'):
        shutil.copy2(p, local / 'raw' / p.name)
    shutil.copy2(ROOT / 'EXPORT_SIDECAR.json', local / 'EXPORT_SIDECAR.json')
    for p in ROOT.glob('PREREG*.json'):
        shutil.copy2(p, local / p.name)
    (local / 'code').mkdir()
    for p in (ROOT / 'code').glob('*.py'):
        shutil.copy2(p, local / 'code' / p.name)
    dump(ROOT / 'raw/LATEST_IMMUTABLE_RUN.json', {'run': stamp, 'arrays': mapping, 'local_results': str(local), 'data_root': str(dest), 'copied_own_bytes': copybytes, 'immutable_source_data': 'Acquired source members kept on original read-only-by-convention paths and verified by SHA; no source trees copied'})
    print('Frozen own outputs', stamp, copybytes, 'bytes')
if __name__ == '__main__':
    main()
