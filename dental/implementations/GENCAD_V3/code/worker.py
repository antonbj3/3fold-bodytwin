"""Participant process: only public inputs, local runtime and output are mounted."""
import sys, os, time, resource, argparse
from pathlib import Path
import numpy as np
from util import read, dump, FAMILIES, LEVELS
from task_io import attach
from legacy.generators import population, parametric
from fast_geometry import optimize
METHODS = {'population': population, 'parametric': parametric, 'constraint_optimizer': optimize}

def run(public, output):
    start = time.perf_counter()
    calls = 0
    for path in sorted((public / 'tasks').glob('*.json')):
        ts = read(path)
        arrays = {name: {'task_ids': np.array([t['task_id'] for t in ts])} for name in METHODS}
        statuses = {name: [] for name in METHODS}
        scenes = {}
        for (i, t) in enumerate(ts):
            loaded = t
            if t['status'] == 'READY':
                fn = t['geometry_file']
                if fn not in scenes:
                    with np.load(public / fn, allow_pickle=False) as a:
                        scenes[fn] = dict(a)
                loaded = attach(t, scenes[fn])
            for (name, fn) in METHODS.items():
                if t['status'] != 'READY':
                    d = {'status': 'UNKNOWN_SITE'}
                else:
                    d = fn(loaded)
                    calls += 1
                statuses[name].append(d['status'])
                if d['status'] == 'DESIGN':
                    arrays[name]['outer_' + str(i)] = np.asarray(d['outer_vertices'])[:, 2]
                    arrays[name]['inner_' + str(i)] = np.asarray(d['inner_vertices'])[:, 2]
        for name in METHODS:
            dest = output / name
            dest.mkdir(parents=True, exist_ok=True)
            arrays[name]['status'] = np.array(statuses[name])
            np.savez_compressed(dest / (path.stem + '.npz'), **arrays[name])
        print('generated', path.stem, flush=True)
    dump(output / 'COST.json', dict(seconds=time.perf_counter() - start, calls=calls, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--inputs', required=True)
    p.add_argument('--output', required=True)
    a = p.parse_args()
    run(Path(a.inputs), Path(a.output))
