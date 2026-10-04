"""Isolated generator: only public scenes and frozen numerical models are mounted."""
import os, sys, time, resource, json, argparse, concurrent.futures, collections
from pathlib import Path
import numpy as np
sys.path.insert(0, '/baseline')
from task_io import attach
from legacy.generators import population, parametric, crown_loop, antagonist, submit
from fast_geometry import optimize
from generator import conditional, branch
INPUTS = None
MODELS = {}
CFG = {}

def casework(path):
    start = time.perf_counter()
    ts = json.loads(path.read_text())
    names = CFG['participants']
    arrays = {name: {'task_ids': np.array([t['task_id'] for t in ts])} for name in names}
    status = {name: [] for name in names}
    scenes = {}
    outer_x18 = {}
    cost = collections.Counter()
    calls = collections.Counter()
    for (i, t) in enumerate(ts):
        loaded = t
        if t['status'] == 'READY':
            fn = t['geometry_file']
            if fn not in scenes:
                with np.load(INPUTS / fn, allow_pickle=False) as a:
                    scenes[fn] = dict(a)
            loaded = attach(t, scenes[fn])
        for name in names:
            begin = time.perf_counter()
            if t['status'] != 'READY':
                d = {'status': 'UNKNOWN_SITE'}
            elif name in CFG['model_methods']:
                d = conditional(loaded, MODELS[CFG['model_methods'][name], t['family']])
            elif name in CFG['branch_methods']:
                b = CFG['branch_methods'][name]
                d = branch(loaded, MODELS[b['shape'], t['family']], MODELS[b['contact'], t['family']], b['beta'])
            elif name == 'population':
                d = population(loaded)
            elif name == 'parametric':
                d = parametric(loaded)
            elif name == 'constraint_optimizer':
                d = optimize(loaded)
            elif name == 'x18_antagonist':
                fn = t['geometry_file']
                if fn not in outer_x18:
                    q = antagonist(loaded)
                    outer_x18[fn] = np.asarray(q['outer_vertices'])[:, 2]
                d = submit(loaded, outer_x18[fn])
            elif name == 'x1b_v2_radius':
                inp = dict(loaded)
                inp['requirements'] = dict(inp['requirements'], tool_radius_mm=CFG['x1b_radius_mm'][t['level']])
                d = crown_loop(inp)
            else:
                raise ValueError(name)
            cost[name] += time.perf_counter() - begin
            calls[name] += 1
            status[name].append(d['status'])
            if d['status'] == 'DESIGN':
                arrays[name]['outer_' + str(i)] = np.asarray(d['outer_vertices'])[:, 2]
                arrays[name]['inner_' + str(i)] = np.asarray(d['inner_vertices'])[:, 2]
    for name in names:
        arrays[name]['status'] = np.array(status[name])
    return (path.stem, arrays, dict(cost), dict(calls), resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, time.perf_counter() - start)

def run(inputs, models, output, config):
    global INPUTS, MODELS, CFG
    INPUTS = inputs
    CFG = json.loads(Path(config).read_text())
    start = time.perf_counter()
    cost = collections.Counter()
    calls = collections.Counter()
    rss = 0.0
    cpu_worker = 0.0
    for p in sorted(models.glob('*/*.npz')):
        with np.load(p, allow_pickle=False) as a:
            MODELS[p.parent.name, p.stem] = dict(a)
    output.mkdir(exist_ok=True, parents=True)
    for name in CFG['participants']:
        (output / name).mkdir(exist_ok=True)
    with concurrent.futures.ProcessPoolExecutor(max_workers=CFG['max_cpu_workers']) as pool:
        for (idx, (key, arrays, cc, nn, rr, ss)) in enumerate(pool.map(casework, sorted((inputs / 'tasks').glob('*.json')), chunksize=1)):
            for (name, a) in arrays.items():
                np.savez_compressed(output / name / (key + '.npz'), **a)
            cost.update(cc)
            calls.update(nn)
            rss = max(rss, rr)
            cpu_worker += ss
            receipt = dict(seconds=time.perf_counter() - start, sum_worker_wall_seconds=cpu_worker, cost_seconds=dict(cost), calls=dict(calls), peak_worker_rss_MiB=rss, peak_parent_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, cases=idx + 1, max_workers=CFG['max_cpu_workers'], numeric_threads_each=1, scope='Generation only; initial labels/templates UNKNOWN cost, all development fit/discovery separately recorded')
            (output / 'COST.json').write_text(json.dumps(receipt, indent=2) + '\n')
            if (idx + 1) % 5 == 0:
                print('generated', idx + 1, 'case', key, 'seconds', round(receipt['seconds']), flush=True)
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--inputs', required=True)
    ap.add_argument('--models', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--config', required=True)
    a = ap.parse_args()
    run(Path(a.inputs), Path(a.models), Path(a.output), a.config)
