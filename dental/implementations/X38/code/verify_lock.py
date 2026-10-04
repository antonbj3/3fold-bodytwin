"""Refuse reruns with changed frozen numerical implementation/input artifacts."""
from pathlib import Path
import hashlib, json, sys
if __name__ == '__main__':
    R = Path(__file__).resolve().parents[1]
    p = R / 'CODE_LOCK.json'
    if not p.exists():
        raise RuntimeError('CODE_LOCK.json missing')
    for q in json.loads(p.read_text())['files']:
        if hashlib.sha256((R / q['path']).read_bytes()).hexdigest() != q['sha256']:
            raise RuntimeError('code hash changed: ' + q['path'])
    import numpy, scipy, trimesh, statsmodels
    actual = {'numpy': numpy.__version__, 'scipy': scipy.__version__, 'trimesh': trimesh.__version__, 'statsmodels': statsmodels.__version__}
    if actual != json.loads((R / 'RUNTIME_LOCK.json').read_text())['packages']:
        raise RuntimeError('numerical runtime version mismatch; review before unpinned replay')
    if sys.version != json.loads((R / 'RUNTIME_LOCK.json').read_text())['python']:
        raise RuntimeError('Python interpreter version mismatch')
    print('Frozen numerical code and Python runtime verified.')
