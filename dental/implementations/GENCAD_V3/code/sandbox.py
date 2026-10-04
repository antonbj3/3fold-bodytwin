import subprocess, shutil, os, sys, time
from pathlib import Path
from util import *

def public_code():
    dest = PAYLOAD / 'participant_code'
    dest.mkdir(exist_ok=True)
    for name in ['util.py', 'task_io.py', 'worker.py', 'fast_geometry.py']:
        shutil.copy2(ROOT / 'code' / name, dest / name)
    if (dest / 'legacy').exists():
        shutil.rmtree(dest / 'legacy')
    shutil.copytree(ROOT / 'code/legacy', dest / 'legacy', ignore=shutil.ignore_patterns('__pycache__'))
    return dest

def command(output, argv=None, inputs=None):
    cmd = [str(PAYLOAD / 'runtime/bin/bwrap'), '--die-with-parent', '--unshare-all', '--new-session', '--clearenv', '--ro-bind', '/usr', '/usr', '--ro-bind', '/lib', '/lib', '--ro-bind', '/lib64', '/lib64', '--proc', '/proc', '--dev', '/dev', '--tmpfs', '/tmp', '--dir', '/work', '--chdir', '/work', '--ro-bind', str(PAYLOAD / 'runtime'), '/runtime', '--ro-bind', str(PAYLOAD / 'participant_code'), '/code', '--ro-bind', str(inputs or PAYLOAD / 'public'), '/inputs', '--bind', str(output), '/output', '--setenv', 'HOME', '/tmp', '--setenv', 'PYTHONDONTWRITEBYTECODE', '1', '--setenv', 'PYTHONWARNINGS', 'ignore']
    for name in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'THREADS']:
        cmd += ['--setenv', name, '2']
    return cmd + ['/runtime/bin/python3', '-s', '-E'] + (argv or ['/code/worker.py', '--inputs', '/inputs', '--output', '/output'])

def generate(output, receipt_name='GENERATION_PROCESS.json', log_name='generation.log'):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    cmd = command(output)
    start = time.perf_counter()
    with (ROOT / 'raw' / log_name).open('w') as f:
        p = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, timeout=14400)
    dump(ROOT / 'raw' / receipt_name, dict(argv=cmd, exit_code=p.returncode, seconds=time.perf_counter() - start))
    if p.returncode:
        raise RuntimeError('participant failed, see raw/generation.log')
