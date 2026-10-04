import os
for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    os.environ[key] = '4'
os.environ['MPLBACKEND'] = 'Agg'
import argparse, copy, importlib, json, time, resource, unittest, io, subprocess, sys
from pathlib import Path
import numpy as np
from gencad_bench.io import ROOT, DATA, dump, sha, digest, now, state
from gencad_bench.schema import validate_task, validate_split
from gencad_bench.generators import BASELINES
from gencad_bench.evaluate import evaluate
from gencad_bench.geometry import shell
from gencad_bench.data_local import prepare
from gencad_bench.literature import extract
from gencad_bench.calibration import predict_loso, score_loso

def tests(out):
    stream = io.StringIO()
    suite = unittest.defaultTestLoader.discover(str(ROOT / 'tests'))
    t = time.perf_counter()
    r = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    (out / 'tests.log').write_text(stream.getvalue())
    result = dict(tests=r.testsRun, failures=len(r.failures), errors=len(r.errors), wall_s=time.perf_counter() - t)
    dump(out / 'tests.json', result)
    if not r.wasSuccessful():
        raise RuntimeError('Core tests failed; see ' + str(out / 'tests.log'))
    return result

def material_shape_summary(rows, name):
    out = {}
    for mat in sorted({r['material'] for r in rows}):
        shapes = [r['shape'] for r in rows if r['generator'] == name and r['material'] == mat and r['shape']]
        out[mat] = dict(RMS_mm=float(np.mean([s['symmetric_sampled_RMS_mm'] for s in shapes])) if shapes else None, p95_mm=float(np.mean([s['symmetric_sampled_p95_mm'] for s in shapes])) if shapes else None, n_scored=len(shapes), n_attempted=sum((r['generator'] == name and r['material'] == mat for r in rows)))
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--plugin', help='module:generate, trusted local Python code')
    ap.add_argument('--prepare', action='store_true', help='Rebuild local inputs; frozen manifest must match')
    ap.add_argument('--skip-context', action='store_true', help='Skip standalone registered-bite import demonstration')
    args = ap.parse_args()
    tic = time.perf_counter()
    stamp = now().replace(':', '').replace('+', '_')
    out = ROOT / 'raw/runs' / stamp
    out.mkdir(parents=True)
    state('R1_RUNNING', 'PREREG hash checked next', 'Run exact tests and all 72 baseline-task pairs')
    if sha(ROOT / 'PREREG_R1.json') != (ROOT / 'PREREG_R1.sha256').read_text().strip():
        raise ValueError('PREREG drift')
    testresult = tests(out)
    print('tests', testresult, flush=True)
    if args.prepare or not (ROOT / 'data/tasks.json').exists():
        prepare()
    frozen = json.loads((ROOT / 'data/TASKSET_FROZEN.json').read_text())['payload']
    if sha(ROOT / 'data/tasks.json') != frozen['tasks_sha256']:
        raise ValueError('Task-set drift')
    tasks = json.loads((ROOT / 'data/tasks.json').read_text())
    validate_split(tasks)
    private = json.loads((ROOT / 'data/evaluator_private.json').read_text())
    if digest(private) != digest(frozen['references']):
        raise ValueError('Reference manifest drift')
    generators = dict(BASELINES)
    if args.plugin:
        (mod, func) = args.plugin.split(':', 1)
        generators[args.plugin] = getattr(importlib.import_module(mod), func)
    rows = []
    cache = {}
    exported = []
    for (index, task) in enumerate(tasks):
        validate_task(task)
        ref = private[task['geometry_id']]
        if sha(ref['path']) != ref['sha256']:
            raise ValueError('Reference geometry drift')
        with np.load(ref['path']) as f:
            reference = {k: f[k] for k in ['vertices', 'faces']}
        for (name, generate) in generators.items():
            given = copy.deepcopy(task)
            before = digest(given)
            t = time.perf_counter()
            try:
                design = generate(given)
                generation = time.perf_counter() - t
                if digest(given) != before:
                    raise ValueError('Generator mutated task')
                result = evaluate(task, design, reference, cache)
                result['design_sha256'] = digest(design)
                dump(out / (task['task_id'] + '__' + name.replace(':', '_') + '_design.json'), design)
                if index == 0 and design['kind'] == 'radial_crown_v1':
                    path = out / (name.replace(':', '_') + '.stl')
                    m = shell(design['outer'], design['intaglio'])
                    m.export(path)
                    error = float(np.linalg.norm(m.vertices - m.vertices.astype(np.float32), axis=1).max())
                    exported.append(dict(generator=name, path=str(path), sha256=sha(path), max_vertex_quantization_mm=error, status='digital candidate; no clinical or CAM release; analytic/exact score is on design JSON'))
            except Exception as e:
                generation = time.perf_counter() - t
                result = dict(level1=dict(status='INVALID', score=0, eligible=False), shape=None, checks={}, error=type(e).__name__ + ': ' + str(e), timings_s={})
            row = dict(task_id=task['task_id'], geometry_id=task['geometry_id'], split=task['split'], tooth_FDI=task['tooth_FDI'], material=task['material']['product'], generator=name, generation_s=generation, **result)
            rows.append(row)
            dump(out / 'baseline_rows_partial.json', rows)
        if (index + 1) % 3 == 0:
            print('tasks', index + 1, '/', len(tasks), 'last', task['geometry_id'], flush=True)
            state('R1_BASELINES_RUNNING', f'{index + 1}/{len(tasks)} tasks recorded', 'Continue complete frozen task set')
    extract()
    predict_loso()
    calibration = score_loso()
    context = None
    if not args.skip_context:
        from gencad_bench.registered_bite import run as context_run
        context = context_run()
    summary = {}
    for name in generators:
        rr = [r for r in rows if r['generator'] == name]
        shapes = [r['shape'] for r in rr if r['shape'] is not None]
        unique = {r['geometry_id']: r['shape'] for r in rr if r['shape'] is not None}
        summary[name] = dict(n_tasks=len(rr), **{s: sum((r['level1']['status'] == s for r in rr)) for s in ['PASS', 'FAIL', 'UNKNOWN', 'INVALID']}, n_anatomies=len(unique), check_counts={k: {s: sum((r.get('checks', {}).get(k, {}).get('status', 'UNKNOWN') == s for r in rr)) for s in ['PASS', 'FAIL', 'UNKNOWN', 'INVALID']} for k in task['requirements']['required_checks']})
    for name in summary:
        summary[name]['shape_aggregation'] = 'Eight anatomical references per material; no material-replica independence assumed'
        summary[name]['by_material'] = material_shape_summary(rows, name)
    result = dict(schema='DentalGenCAD-Bench-0.1', run_utc=now(), review_state='PENDING_INDEPENDENT_REVIEW', tests=testresult, n_tasks=len(tasks), n_designs=len(rows), subject_count=len({t['subject_group'] for t in tasks}), taskset_sha256=sha(ROOT / 'data/tasks.json'), summary=summary, rows=rows, exports=exported, calibration_summary=calibration['summary'], calibration_admission=calibration['calibration_admission'], registered_bite=context, cost=dict(run_wall_s=time.perf_counter() - tic, peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, data_preparation=json.loads((ROOT / 'raw/data_preparation.json').read_text())['preparation_wall_s'], source_reading_and_reasoning='UNKNOWN; no total discovery speedup claimed', threads=4), external_referent=dict(kind='published_dataset', locator='https://osf.io/xctdy/; data/evaluator_private.json with native file hashes', compared_quantity='withheld native labelled tooth surface, mm; not optimal crown design', refutes_us=True), external_referents=[calibration['external_referent']], conclusion='Scoped digital feasibility and anatomical reconstruction; no clinical-quality ranking or 10x claim', run_directory=str(out))
    dump(out / 'results.json', result)
    dump(ROOT / 'results.json', result)
    subprocess.run([sys.executable, '-s', str(ROOT / 'scripts/render_report.py')], check=True)
    state('R1_COMPLETE', 'All frozen baseline tasks evaluated; see results.json for failures and UNKNOWN', 'Preserve round-1 handoff; attack first structural failure with round-2 preregistration')
    print(json.dumps(dict(summary=summary, cost=result['cost']), indent=2), flush=True)
if __name__ == '__main__':
    main()
