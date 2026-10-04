"""Execute unchanged full-kernel participant over verified public bytes."""
import sys, json, time, io, hashlib, runpy
from pathlib import Path
import numpy as np
lock = json.loads(Path('/public_lock.json').read_text())
read_text = Path.read_text
load = np.load

def safe_blob(p):
    rel = str(p)[len('/inputs/'):]
    entry = lock[rel]
    for i in range(4):
        try:
            b = Path(p).read_bytes()
            if len(b) != entry['bytes'] or hashlib.sha256(b).hexdigest() != entry['sha256']:
                raise ValueError('immutable public byte drift ' + rel)
            return b
        except FileNotFoundError:
            if i == 3:
                raise
            time.sleep(0.5)

def safe_text(self, *args, **kw):
    if str(self).startswith('/inputs/'):
        return safe_blob(self).decode(kw.get('encoding') or 'utf-8')
    return read_text(self, *args, **kw)

def safe_load(file, *args, **kw):
    if isinstance(file, (str, Path)) and str(file).startswith('/inputs/'):
        file = io.BytesIO(safe_blob(file))
    return load(file, *args, **kw)
Path.read_text = safe_text
np.load = safe_load
sys.argv[0] = '/candidate/participant.py'
runpy.run_path('/candidate/participant.py', run_name='__main__')
