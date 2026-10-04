import subprocess, time, resource, shlex, sys, os
from common import *
cmd = sys.argv[1:]
run = ROOT / 'runs' / ('run_%03d' % (len(list((ROOT / 'runs').glob('run_*'))) + 1))
run.mkdir()
env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1', NUMBA_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1', MPLCONFIGDIR=str(ROOT / 'runs/mplconfig'))
utc0 = now()
t = time.perf_counter()
cpu0 = resource.getrusage(resource.RUSAGE_CHILDREN)
result = subprocess.run(cmd, cwd=ROOT, env=env, stdout=(run / 'stdout').open('w'), stderr=(run / 'stderr').open('w'))
cpu1 = resource.getrusage(resource.RUSAGE_CHILDREN)
record = {'command': shlex.join(cmd), 'cwd': str(ROOT), 'exitcode': result.returncode, 'wall_seconds': time.perf_counter() - t, 'cpu_seconds': cpu1.ru_utime + cpu1.ru_stime - cpu0.ru_utime - cpu0.ru_stime, 'peak_child_RSS_KiB': cpu1.ru_maxrss, 'started_utc': utc0, 'completed_utc': now(), 'thread_limits': 1}
(run / 'record.json').write_text(json.dumps(record, indent=2) + '\n')
with (ROOT / 'COMMANDS.md').open('a') as f:
    f.write('\n```bash\n' + shlex.join(cmd) + '\n```\nexit ' + str(result.returncode) + '; logs `' + str(run.relative_to(ROOT)) + '/{stdout,stderr,record.json}`.\n')
print(json.dumps(record))
print((run / 'stdout').read_text()[-1400:])
print((run / 'stderr').read_text()[-2400:])
sys.exit(result.returncode)
