from pathlib import Path
import datetime, hashlib, json, os, resource, shutil, subprocess, sys, time
package = Path(__file__).resolve().parents[1]
record = '--record' in sys.argv[1:]
if record:
    if (package / 'results.json').exists():
        raise SystemExit('Archived results already exist; run without --record for a replay')
    out = package
else:
    out = package / 'replays' / datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    out.mkdir(parents=True)
    for pat in ['PREREG*.json', 'PREREG*.sha256', 'DECOMPOSITION*.json', 'FROZEN*.json', 'FROZEN*.sha256']:
        for p in package.glob(pat):
            shutil.copy2(p, out / p.name)
    for name in ['RESULTS_R3.json', 'CURRENT_WORK_STATE.json', 'GOAL_TO_STRESS_POINT.md']:
        shutil.copy2(package / name, out / name)
env = dict(os.environ, X10_OUTPUT_DIR=str(out), OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
(out / 'raw').mkdir(exist_ok=True)
start = time.perf_counter()
steps = []
for name in ['prepare_sources.py', 'run_r1.py', 'run_r2.py', 'run_r4.py', 'export_fixture.py', 'validate_exports.py', 'plot_demo.py']:
    t = time.perf_counter()
    se = dict(env)
    if name == 'plot_demo.py':
        se['PYTHONNOUSERSITE'] = '1'
    cmd = [sys.executable, str(package / 'code' / name)]
    p = subprocess.run(cmd, cwd=package, env=se, text=True, capture_output=True)
    (out / 'raw' / (name + '.stdout.txt')).write_text(p.stdout)
    (out / 'raw' / (name + '.stderr.txt')).write_text(p.stderr)
    steps.append({'command': cmd, 'elapsed_seconds': time.perf_counter() - t, 'returncode': p.returncode, 'environment_overrides': {k: se[k] for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'X10_OUTPUT_DIR']}})
    print(p.stdout.strip())
    if p.returncode:
        (out / 'RUN_FAILURE.json').write_text(json.dumps({'steps': steps, 'stderr': p.stderr}, indent=2) + '\n')
        print(p.stderr, file=sys.stderr)
        raise SystemExit(p.returncode)
ledger = {'steps': steps, 'pre_assembly_wall_seconds': time.perf_counter() - start, 'max_child_rss_bytes': resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * 1024, 'gpu': False, 'threads': 1, 'external_physical_measurement': 'NOT_RUN'}
(out / 'RUN_LEDGER.json').write_text(json.dumps(ledger, indent=2) + '\n')
assembly = subprocess.run([sys.executable, str(package / 'code/assemble.py')], cwd=package, env=env)
if assembly.returncode:
    raise SystemExit(assembly.returncode)
(out / 'RUN_COMPLETION.json').write_text(json.dumps({'full_command_wall_seconds': time.perf_counter() - start, 'max_child_rss_bytes': resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * 1024, 'results_sha256': hashlib.sha256((out / 'results.json').read_bytes()).hexdigest(), 'completed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat()}, indent=2) + '\n')
print('Output:', out)
