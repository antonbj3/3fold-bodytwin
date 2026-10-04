import json, os, resource, subprocess, sys, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
env = os.environ.copy()
for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    env[k] = '4'
env['MPLCONFIGDIR'] = str(ROOT / 'raw/matplotlib_cache')
start = time.perf_counter()
steps = []
for name in ['run_r1.py', 'run_r2.py', 'run_r3.py', 'validate.py']:
    stamp = time.perf_counter()
    p = subprocess.run([sys.executable, str(ROOT / 'code' / name)], cwd=ROOT, env=env)
    steps.append({'command': [sys.executable, str(ROOT / 'code' / name)], 'exit_code': p.returncode, 'wall_seconds': time.perf_counter() - stamp})
    if p.returncode:
        raise SystemExit(p.returncode)
cost = {'steps': steps, 'wall_seconds_before_render': time.perf_counter() - start, 'max_child_rss_MiB': resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024, 'threads_max': 4, 'gpu': False, 'new_physical_observations': 0, 'physical_cost': 'UNKNOWN', 'lab_cost_ranges': 'PHENOMENOLOGICAL_X71_ASSUMPTIONS'}
(ROOT / 'COST_RUN.json').write_text(json.dumps(cost, indent=2) + '\n')
plotenv = env.copy()
plotenv['PYTHONNOUSERSITE'] = '1'
stamp = time.perf_counter()
p = subprocess.run([sys.executable, str(ROOT / 'code/report.py')], cwd=ROOT, env=plotenv)
cost['steps'].append({'command': [sys.executable, str(ROOT / 'code/report.py')], 'exit_code': p.returncode, 'wall_seconds': time.perf_counter() - stamp, 'PYTHONNOUSERSITE': '1'})
cost['wall_seconds_total'] = time.perf_counter() - start
cost['max_child_rss_MiB'] = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024
(ROOT / 'COST_RUN.json').write_text(json.dumps(cost, indent=2) + '\n')
raise SystemExit(p.returncode)
