"""Execution wrapper around unchanged frozen participant, verified one-case snapshots.

Only public bytes are included, not original source trees or referee manifests.
Temporary snapshots are deleted after each case and never cache the corpus.
"""
import sys, os, time, json, tempfile, hashlib, resource, collections, concurrent.futures, argparse
from pathlib import Path
import numpy as np
import participant as p
LOCK = {}

def blob(path, entry):
    for attempt in range(4):
        try:
            b = path.read_bytes()
            if len(b) != entry['bytes'] or hashlib.sha256(b).hexdigest() != entry['sha256']:
                raise ValueError('public input drift: ' + str(path))
            return b
        except FileNotFoundError:
            if attempt == 3:
                raise
            time.sleep(0.5)

def guarded_case(path):
    original = p.INPUTS
    with tempfile.TemporaryDirectory(prefix='verified_case_', dir='/tmp') as tmp:
        tmp = Path(tmp)
        rel = 'tasks/' + path.name
        task_blob = blob(path, LOCK[rel])
        ts = json.loads(task_blob)
        (tmp / 'tasks').mkdir()
        (tmp / 'scenes').mkdir()
        (tmp / rel).write_bytes(task_blob)
        for name in {t['geometry_file'] for t in ts if t['status'] == 'READY'}:
            (tmp / name).write_bytes(blob(original / name, LOCK[name]))
        p.INPUTS = tmp
        try:
            return p.casework(tmp / rel)
        finally:
            p.INPUTS = original

def run(inputs, models, output, config, lock):
    global LOCK
    LOCK = json.loads(Path(lock).read_text())
    p.INPUTS = inputs
    p.CFG = json.loads(Path(config).read_text())
    start = time.perf_counter()
    cost = collections.Counter()
    calls = collections.Counter()
    rss = 0.0
    n = 0
    for f in models.glob('*/*.npz'):
        with np.load(f, allow_pickle=False) as a:
            p.MODELS[f.parent.name, f.stem] = dict(a)
    paths = sorted((inputs / 'tasks').glob('*.json'))
    pending = [f for f in paths if not all(((output / name / (f.stem + '.npz')).exists() for name in p.CFG['participants']))]
    print('resuming pending', len(pending), 'completed', len(paths) - len(pending), flush=True)
    with concurrent.futures.ProcessPoolExecutor(max_workers=3) as pool:
        for (key, arrays, cc, nn, rr, seconds) in pool.map(guarded_case, pending, chunksize=1):
            for (name, a) in arrays.items():
                np.savez_compressed(output / name / (key + '.npz'), **a)
            n += 1
            cost.update(cc)
            calls.update(nn)
            rss = max(rss, rr)
            if n % 5 == 0:
                print('resumed', n, 'total', len(paths) - len(pending) + n, 'seconds', round(time.perf_counter() - start), flush=True)
    receipt = dict(seconds=time.perf_counter() - start, cost_seconds=dict(cost), calls=dict(calls), cases_resumed=n, original_completed=len(paths) - len(pending), peak_worker_rss_MiB=rss, snapshot_peak='one case public arrays per worker only; deleted after return', candidate_sha_preserved=True)
    (output / 'COST.json').write_text(json.dumps(receipt, indent=2) + '\n')
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--inputs')
    ap.add_argument('--models')
    ap.add_argument('--output')
    ap.add_argument('--config')
    ap.add_argument('--lock')
    a = ap.parse_args()
    run(Path(a.inputs), Path(a.models), Path(a.output), a.config, a.lock)
