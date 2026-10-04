from dental_release.paths import expand as _release_expand
from pathlib import Path
import json, hashlib, datetime, time, resource, os
L = Path(__file__).resolve().parents[1]
D = L / 'local_data'
ROOT = Path(_release_expand('@DENTAL_INPUT_ROOT@/workspace'))

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda : f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def write(p, o):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(o, indent=2, ensure_ascii=False, allow_nan=False, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')

def verify(n):
    p = L / n
    expected = (L / (n + '.sha256')).read_text().strip()
    assert sha(p) == expected, (n, 'frozen hash drift')
    return json.loads(p.read_text())

def state(stage, gate, nxt):
    write(L / 'CURRENT_WORK_STATE.json', dict(lane='X44-missing-chains', status=stage, latest_gate=gate, next_operation=nxt, updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))

def footprint():
    return sum((x.stat().st_size for x in D.rglob('*') if x.is_file()))

def cost(st):
    return dict(wall_seconds=time.perf_counter() - st, max_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, threads=4, gpu_seconds=0, intermediate_bytes=footprint())
