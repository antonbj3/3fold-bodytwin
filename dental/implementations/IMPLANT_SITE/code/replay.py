import os, pathlib, subprocess, sys, time, json, resource, datetime
P = pathlib.Path(__file__).resolve().parents[1]
for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    os.environ[k] = '1'
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
start = time.perf_counter()
steps = []
for (name, cmd) in [('R1', [sys.executable, 'code/run_r1.py']), ('R2', [sys.executable, 'code/run_r2.py']), ('R3', [sys.executable, 'code/run_r3.py']), ('clinical', [sys.executable, 'code/clinical_reference.py']), ('verification', [sys.executable, 'code/verify.py']), ('report', [sys.executable, 'code/report.py']), ('figure', ['/usr/bin/python3', '-s', 'code/plot.py'])]:
    t = time.perf_counter()
    log = P / 'raw' / f'replay_{name}.log'
    with log.open('w') as f:
        r = subprocess.run(cmd, cwd=P, stdout=f, stderr=subprocess.STDOUT)
    steps.append({'step': name, 'command': cmd, 'returncode': r.returncode, 'seconds': time.perf_counter() - t, 'log': str(log)})
    print(name, 'PASS' if r.returncode == 0 else 'FAIL', round(steps[-1]['seconds'], 3), 's', flush=True)
    if r.returncode:
        print(log.read_text()[-5000:])
        sys.exit(r.returncode)
cost = {'timestamp_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'steps': steps, 'wall_seconds': time.perf_counter() - start, 'child_peak_rss_kib': resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss, 'threads': 1, 'gpu': False}
(P / 'raw/run_all_cost.json').write_text(json.dumps(cost, indent=2) + '\n')
subprocess.run([sys.executable, 'code/finalize.py'], cwd=P, check=True)
print('All steps finished. Physical clinical plan remains UNKNOWN. See README_DEMO.md.')
