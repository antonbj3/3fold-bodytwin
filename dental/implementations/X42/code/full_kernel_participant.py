"""Additional frozen standard control; public-only inputs and models."""
import sys, json, time, resource, argparse
from pathlib import Path
import numpy as np
sys.path.insert(0, '/baseline')
from task_io import attach
from generator import branch

def run(inputs, models, contacts, output):
    shape = {}
    contact = {}
    for p in models.glob('*.npz'):
        with np.load(p, allow_pickle=False) as a:
            shape[p.stem] = dict(a)
    for p in contacts.glob('*.npz'):
        with np.load(p, allow_pickle=False) as a:
            contact[p.stem] = dict(a)
    start = time.perf_counter()
    calls = 0
    for (ii, p) in enumerate(sorted((inputs / 'tasks').glob('*.json'))):
        ts = json.loads(p.read_text())
        arrays = {'task_ids': np.array([t['task_id'] for t in ts])}
        statuses = []
        scenes = {}
        for (i, t) in enumerate(ts):
            if t['status'] != 'READY':
                d = {'status': 'UNKNOWN_SITE'}
            else:
                f = t['geometry_file']
                if f not in scenes:
                    with np.load(inputs / f, allow_pickle=False) as a:
                        scenes[f] = dict(a)
                d = branch(attach(t, scenes[f]), shape[t['family']], contact[t['family']], 1.0)
                calls += 1
            statuses.append(d['status'])
            if d['status'] == 'DESIGN':
                arrays['outer_' + str(i)] = np.asarray(d['outer_vertices'])[:, 2]
                arrays['inner_' + str(i)] = np.asarray(d['inner_vertices'])[:, 2]
        arrays['status'] = np.array(statuses)
        np.savez_compressed(output / (p.stem + '.npz'), **arrays)
        if (ii + 1) % 25 == 0:
            print('full kernel', ii + 1, 'seconds', round(time.perf_counter() - start), flush=True)
    return dict(seconds=time.perf_counter() - start, calls=calls, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--inputs', required=True)
    ap.add_argument('--models', required=True)
    ap.add_argument('--contacts', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--cost-output', required=True)
    a = ap.parse_args()
    c = run(Path(a.inputs), Path(a.models), Path(a.contacts), Path(a.output))
    Path(a.cost_output).write_text(json.dumps(c, indent=2) + '\n')
