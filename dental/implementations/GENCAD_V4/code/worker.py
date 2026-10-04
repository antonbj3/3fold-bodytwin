"""Public-only participant runner; no reference or scoring code is mounted."""
import os, sys, json, time, resource, argparse, subprocess
from pathlib import Path
import numpy as np
sys.path.insert(0, '/participants')
sys.path.insert(0, '/participants/field')
from task_io import attach
from legacy import generators as g
from fast_geometry import optimize
import field_generator as field
import x42_generator as x42

def run(inputs, output):
    begin = time.perf_counter()
    cfg = json.loads((inputs / 'PARTICIPANTS.json').read_text())
    names = [n for n in cfg['participants'] if n != 'field_generator']
    models = {}
    for kind in ['kernel', 'kernel_contact']:
        for p in Path('/models', kind).glob('*.npz'):
            with np.load(p, allow_pickle=False) as a:
                models[kind, p.stem] = dict(a)
    costs = {n: 0.0 for n in names}
    errors = []
    diagnostics = []
    for (idx, path) in enumerate(sorted((inputs / 'tasks').glob('*.json'))):
        tasks = json.loads(path.read_text())
        scenes = {}
        cached18 = {}
        arrays = {n: {'task_ids': np.array([t['task_id'] for t in tasks])} for n in names}
        statuses = {n: [] for n in names}
        for (i, t0) in enumerate(tasks):
            t = t0
            if t['status'] == 'READY':
                fn = t['geometry_file']
                if fn not in scenes:
                    with np.load(inputs / fn, allow_pickle=False) as a:
                        scenes[fn] = dict(a)
                t = attach(t, scenes[fn])
            for name in names:
                st = time.perf_counter()
                try:
                    if t['status'] != 'READY':
                        d = {'status': 'UNKNOWN_SITE'}
                    elif name == 'population':
                        d = g.population(t)
                    elif name == 'parametric':
                        d = g.parametric(t)
                    elif name == 'constraint_optimizer':
                        d = optimize(t)
                    elif name == 'field_generator':
                        d = field.generate(t)
                    elif name == 'X1B_morphology':
                        tt = dict(t, requirements=dict(t['requirements'], tool_radius_mm=cfg['x1b_tool_radius_mm'][t['level']]))
                        d = g.crown_loop(tt)
                    elif name == 'X18_antagonist':
                        if fn not in cached18:
                            cached18[fn] = np.asarray(g.antagonist(t)['outer_vertices'])[:, 2]
                        d = g.submit(t, cached18[fn])
                    elif name == 'x42_shape':
                        d = x42.conditional(t, models['kernel', t['family']])
                    elif name == 'x42_contact_branch':
                        d = x42.branch(t, models['kernel', t['family']], models['kernel_contact', t['family']], 1.0)
                    else:
                        raise ValueError(name)
                except Exception as e:
                    d = {'status': 'ERROR'}
                    errors.append(dict(task=t['task_id'], participant=name, error=repr(e)))
                costs[name] += time.perf_counter() - st
                statuses[name].append(d['status'])
                if d['status'] == 'DESIGN':
                    arrays[name]['outer_' + str(i)] = np.asarray(d['outer_vertices'])[:, 2]
                    arrays[name]['inner_' + str(i)] = np.asarray(d['inner_vertices'])[:, 2]
                if name == 'field_generator' and d.get('diagnostics'):
                    diagnostics.append(dict(task=t['task_id'], **d['diagnostics']))
        for name in names:
            arrays[name]['status'] = np.array(statuses[name])
            (output / name).mkdir(exist_ok=True, parents=True)
            np.savez_compressed(output / name / (path.stem + '.npz'), **arrays[name])
        g._x1_cache.clear()
        field._cache.clear()
        receipt = dict(cases_completed=idx + 1, wall_seconds=time.perf_counter() - begin, participant_seconds=costs, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, errors=errors)
        (output / 'COST.json').write_text(json.dumps(receipt, indent=2) + '\n')
        print('generated', idx + 1, path.stem, 'seconds', round(receipt['wall_seconds'], 1), 'errors', len(errors), flush=True)
    subprocess.run([sys.executable, '-s', '-E', '-B', '/runner/field_worker.py'], check=True)
if __name__ == '__main__':
    run(Path('/inputs'), Path('/output'))
