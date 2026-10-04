"""Unchanged Field participant in a process with its own HiGHS scheduler."""
import sys, json, time, resource
from pathlib import Path
import numpy as np
sys.path.insert(0, '/participants')
sys.path.insert(0, '/participants/field')
from task_io import attach
import field_generator as field

def run():
    start = time.perf_counter()
    output = Path('/output/field_generator')
    output.mkdir(exist_ok=True, parents=True)
    reasons = []
    diags = []
    for (idx, path) in enumerate(sorted(Path('/inputs/tasks').glob('*.json'))):
        ts = json.loads(path.read_text())
        a = {'task_ids': np.array([t['task_id'] for t in ts])}
        status = []
        scenes = {}
        for (i, t) in enumerate(ts):
            if t['status'] == 'READY':
                fn = t['geometry_file']
                if fn not in scenes:
                    with np.load(Path('/inputs', fn), allow_pickle=False) as sc:
                        scenes[fn] = dict(sc)
                t = attach(t, scenes[fn])
            d = field.generate(t)
            status.append(d['status'])
            if d['status'] == 'DESIGN':
                a['outer_' + str(i)] = np.asarray(d['outer_vertices'])[:, 2]
                a['inner_' + str(i)] = np.asarray(d['inner_vertices'])[:, 2]
                diags.append(dict(task=t['task_id'], **d['diagnostics']))
            else:
                reasons.append(dict(task=t['task_id'], status=d['status'], reason=d.get('reason')))
        a['status'] = np.array(status)
        np.savez_compressed(output / (path.stem + '.npz'), **a)
        field._cache.clear()
        if idx % 8 == 0:
            print('isolated field', idx + 1, 'seconds', round(time.perf_counter() - start), flush=True)
    Path('/output/FIELD_DIAGNOSTICS.json').write_text(json.dumps(diags, indent=2) + '\n')
    Path('/output/FIELD_REASONS.json').write_text(json.dumps(reasons, indent=2) + '\n')
    Path('/output/FIELD_COST.json').write_text(json.dumps(dict(seconds=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024), indent=2) + '\n')
if __name__ == '__main__':
    run()
