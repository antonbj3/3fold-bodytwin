"""Public-only plugin process; cannot import scorer or read withheld references."""
import os, sys, json, time, resource, importlib.util
from pathlib import Path
import numpy as np
from .tasks import load_task
from .common import dump, clean
from .generators import PARTICIPANTS

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--participant', required=True)
    ap.add_argument('--inputs', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--plugin')
    a = ap.parse_args()
    fn = PARTICIPANTS.get(a.participant)
    if a.plugin:
        (path, name) = a.plugin.rsplit(':', 1)
        spec = importlib.util.spec_from_file_location('external_generator', path)
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        fn = getattr(m, name)
    if fn is None:
        raise ValueError('unknown participant')
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    start = time.perf_counter()
    for p in sorted(Path(a.inputs).glob('*.json')):
        t = load_task(p)
        tic = time.perf_counter()
        try:
            if t['status'] != 'READY':
                d = dict(task_id=t['task_id'], status='ABSTAIN', reason='source site unavailable')
            else:
                d = fn(t)
            dest = out / (t['task_id'] + '.json')
            dump(dest, d)
            if dest.stat().st_size > 2000000:
                raise ValueError('oversized output')
            rows.append(dict(task_id=t['task_id'], status=d.get('status'), seconds=time.perf_counter() - tic))
        except Exception as e:
            dump(out / (t['task_id'] + '.json'), dict(task_id=t['task_id'], status='ERROR', error=type(e).__name__ + ': ' + str(e)))
            rows.append(dict(task_id=t['task_id'], status='ERROR', reason=str(e), seconds=time.perf_counter() - tic))
        if len(rows) % 36 == 0:
            print(a.participant, len(rows), flush=True)
    dump(out / '_cost.json', dict(participant=a.participant, seconds=time.perf_counter() - start, max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, per_task=rows))
if __name__ == '__main__':
    main()
