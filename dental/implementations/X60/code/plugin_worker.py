"""Process-isolated runner, using the existing module.py:callable plugin contract."""
import sys, os, json, time, resource, importlib.util
from pathlib import Path
import numpy as np
sys.path.insert(0, '/participants')
sys.path.insert(0, '/plugin')
from task_io import attach
spec = importlib.util.spec_from_file_location('external_generator', '/plugin/plugin.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
fn = getattr(m, 'generate')

def run():
    cfg = json.loads(Path('/plugin/PANEL_SELECTION.json').read_text())
    out = Path('/output')
    cost = []
    errors = []
    diagnostics = []
    start = time.perf_counter()
    try:
        Path('/private/references/secret.npz').read_bytes()
        private_denied = False
    except FileNotFoundError:
        private_denied = True
    for case in cfg['selected']:
        key = case['case_key']
        tasks = json.loads(Path('/inputs/tasks', key + '.json').read_text())
        arrays = {}
        scenes = {}
        ids = []
        statuses = []
        for t in tasks:
            if t['family'] not in cfg['families']:
                continue
            i = len(ids)
            ids.append(t['task_id'])
            tic = time.perf_counter()
            try:
                if t['status'] != 'READY':
                    d = dict(status='UNKNOWN_SITE', reason=t['site_error'])
                else:
                    name = t['geometry_file']
                    if name not in scenes:
                        with np.load('/inputs/' + name, allow_pickle=False) as a:
                            scenes[name] = dict(a)
                    t = attach(t, scenes[name])
                    d = fn(t)
                statuses.append(d['status'])
                if d['status'] == 'DESIGN':
                    arrays['outer_' + str(i)] = np.asarray(d['outer_vertices'])[:, 2]
                    arrays['inner_' + str(i)] = np.asarray(d['inner_vertices'])[:, 2]
                    diagnostics.append(dict(task_id=t['task_id'], **d['diagnostics']))
                else:
                    diagnostics.append(dict(task_id=t['task_id'], status=d['status'], reason=d.get('reason')))
            except Exception as e:
                statuses.append('ERROR')
                errors.append(dict(task_id=t['task_id'], error=repr(e)))
            cost.append(dict(task_id=t['task_id'], seconds=time.perf_counter() - tic))
        arrays.update(task_ids=np.array(ids), status=np.array(statuses))
        np.savez_compressed(out / (key + '.npz'), **arrays)
        Path(out / 'DIAGNOSTICS.json').write_text(json.dumps(diagnostics, indent=2) + '\n')
        Path(out / 'COST.json').write_text(json.dumps(dict(wall_seconds=time.perf_counter() - start, max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, per_task=cost, errors=errors, private_reference_unmounted=private_denied), indent=2) + '\n')
        print('case', key, 'tasks', len(cost), 'errors', len(errors), 'seconds', round(time.perf_counter() - start, 1), flush=True)
if __name__ == '__main__':
    run()
