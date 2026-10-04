"""One-command replay in a small isolated own-artifact directory."""
from pathlib import Path
import datetime, hashlib, json, os, resource, shutil, subprocess, sys, time
HERE = Path(__file__).resolve().parent

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def read(p):
    return json.loads(p.read_text())
start = time.perf_counter()
preregs = {f'PREREG_R{i}.json': sha(HERE / f'PREREG_R{i}.json') for i in range(1, 6)}
for i in range(1, 6):
    assert read(HERE / f'PREREG_R{i}.sha256')['sha256'] == preregs[f'PREREG_R{i}.json']
lock = read(HERE / 'INPUT_LOCK.json')
for x in lock['inputs']:
    assert sha(HERE / x['path']) == x['sha256'], x['path']
original_plan = next((x for x in read(HERE / 'SOURCE_MANIFEST.json') if x['path'].endswith('DEMO48_PACKAGE/MEASUREMENT_PLAN.md')))
assert sha(Path(original_plan['path'])) == original_plan['sha256'], 'X40 original hash drift'
replay = HERE / 'replay' / datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
replay.mkdir(parents=True)
names = ['planning.py', 'run_r1.py', 'run_r2.py', 'run_r3.py', 'run_r4.py', 'run_r5.py', 'measurement_port.py', 'report.py', 'plot_demo.py', 'INPUT_LOCK.json', 'SOURCE_MANIFEST.json', 'FROZEN_PREDICTIONS.json'] + list(preregs)
for n in names:
    shutil.copy2(HERE / n, replay / n)
shutil.copytree(HERE / 'inputs', replay / 'inputs')
(replay / 'raw').mkdir()
shutil.copy2(HERE / 'raw/INVENTORY_FINAL_SCOPE.json', replay / 'raw/INVENTORY_FINAL_SCOPE.json')
shutil.copy2(HERE / 'raw/INVENTORY_V2_CHANGE.json', replay / 'raw/INVENTORY_V2_CHANGE.json')
env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
logs = []
for script in ['run_r1.py', 'run_r2.py', 'run_r3.py', 'run_r4.py', 'run_r5.py', 'report.py']:
    run = subprocess.run([sys.executable, script], cwd=replay, env=env, text=True, capture_output=True, timeout=60)
    logs.append(dict(script=script, returncode=run.returncode, stdout=run.stdout, stderr=run.stderr))
    if run.returncode:
        raise RuntimeError(run.stdout + run.stderr)
run = subprocess.run(['/usr/bin/python3', '-s', 'plot_demo.py'], cwd=replay, env=env, text=True, capture_output=True, timeout=30)
logs.append(dict(script='plot_demo.py', returncode=run.returncode, stdout=run.stdout, stderr=run.stderr))
if run.returncode:
    raise RuntimeError(run.stderr)
base = read(HERE / 'results.json')
actual = read(replay / 'results.json')
fields = ['conditional_first_week_evaluable_queries', 'physically_new_decisions', 'full_chains_closed', 'first_week_crown72_only', 'sufficiency', 'inventory']
checks = {k: actual[k] == base[k] for k in fields}
checks['preregs_unchanged'] = all((sha(HERE / n) == h for (n, h) in preregs.items()))
checks['original_X40_unchanged'] = sha(Path(original_plan['path'])) == original_plan['sha256']
checks['figure_exists'] = (replay / 'MEASUREMENT_PRIORITY.png').stat().st_size > 1000
checks['source_small_intermediates'] = sum((p.stat().st_size for p in HERE.rglob('*') if p.is_file())) < 3000000000
receipt = dict(status='PASS' if all(checks.values()) else 'FAIL', checks=checks, replay_dir=str(replay), runtime_s=time.perf_counter() - start, max_child_RSS_MiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024, threads=1, network=False, physical_measurement=False, logs=logs)
(HERE / 'DELIVERY_VERIFICATION.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(receipt['status'], 'crown72 pilot5 potential queries; physical decisions0; exact summary difference24')
print('Figure:', replay / 'MEASUREMENT_PRIORITY.png')
if not all(checks.values()):
    raise SystemExit(1)
