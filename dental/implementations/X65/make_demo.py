"""One-command offline replay. All required scientific inputs are packaged."""
import datetime, hashlib, json, os, resource, subprocess, sys, time
from pathlib import Path
ROOT = Path(__file__).resolve().parent
os.chdir(ROOT)

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    start = time.perf_counter()
    os.environ.update(OMP_NUM_THREADS='4', OPENBLAS_NUM_THREADS='4', MKL_NUM_THREADS='4', NUMEXPR_NUM_THREADS='4', MPLBACKEND='Agg')
    (ROOT / 'raw').mkdir(exist_ok=True)
    manifest = json.loads((ROOT / 'inputs/SOURCE_MANIFEST.json').read_text())
    checked = []
    for r in manifest['files']:
        if 'path' not in r:
            continue
        p = ROOT / r['path']
        match = sha(p) == r['sha256']
        checked.append({'path': r['path'], 'sha256': sha(p), 'match': match})
    if not all((r['match'] for r in checked)):
        raise ValueError('packaged scientific source hash changed')
    (ROOT / 'raw/PACKAGE_VERIFY.json').write_text(json.dumps(checked, indent=2) + '\n')
    costs = []
    for script in ['extract_primary.py', 'run_r1.py', 'run_r2.py', 'run_r3.py', 'run_r4.py']:
        t = time.perf_counter()
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        with (ROOT / 'raw' / ('run_' + script.replace('.py', '.log'))).open('w') as log:
            proc = subprocess.run([sys.executable, str(ROOT / 'code' / script)], stdout=log, stderr=subprocess.STDOUT)
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        costs.append({'script': script, 'wall_s': time.perf_counter() - t, 'cpu_s': after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime, 'largest_child_RSS_MiB': after.ru_maxrss / 1024, 'exit_code': proc.returncode})
        print(f"{script}: exit {proc.returncode}, {costs[-1]['wall_s']:.3f}s", flush=True)
        if proc.returncode:
            raise RuntimeError(f"{script} failed; see raw/run_{script.replace('.py', '.log')}")
    cost = {'created_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'stages': costs, 'numerical_pipeline_wall_s': time.perf_counter() - start, 'cpu_s': sum((c['cpu_s'] for c in costs)), 'peak_child_RSS_MiB': resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024, 'threads_per_process': 4, 'gpu': False, 'packaged_input_bytes': manifest['total_bytes'], 'one_time_selective_source_preparation_wall_s': manifest['wall_s'], 'agent_discovery_cost': 'UNKNOWN; file searches and code construction documented', 'human_questions': 0, 'fallback': 'abstain; no further laboratory measurement performed', 'future_acquisition_cost': 'UNKNOWN until same-montage protocol and instrument plan are fixed'}
    (ROOT / 'raw/COST.json').write_text(json.dumps(cost, indent=2) + '\n')
    proc = subprocess.run([sys.executable, str(ROOT / 'code/report.py')])
    if proc.returncode:
        raise RuntimeError('report generation failed')
    print('Demo finished: README_DEMO.md, RESULTS.md, results.json, figure.png', flush=True)
if __name__ == '__main__':
    main()
