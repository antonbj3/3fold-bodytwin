"""A local namespace boundary: generators see no archive, private reference, or scorer."""
import os, subprocess, time, sys, site, shutil
from pathlib import Path
from .common import *

def command(inputs, output, argv, plugin=None):
    cmd = ['bwrap', '--die-with-parent', '--unshare-all', '--new-session', '--ro-bind', '/usr', '/usr', '--ro-bind', '/lib', '/lib', '--ro-bind', '/lib64', '/lib64', '--proc', '/proc', '--dev', '/dev', '--tmpfs', '/tmp', '--dir', '/home', '--dir', '/work', '--chdir', '/work']
    bundle = ROOT / 'runtime_public'
    cmd += ['--ro-bind', str(bundle), '/code', '--ro-bind', str(inputs), '/inputs', '--bind', str(output), '/output']
    user_site = site.getusersitepackages()
    cmd += ['--ro-bind', user_site, '/pydeps', '--setenv', 'PYTHONPATH', '/code:/pydeps', '--setenv', 'HOME', '/tmp', '--setenv', 'PYTHONDONTWRITEBYTECODE', '1']
    for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'THREADS']:
        cmd += ['--setenv', k, '4']
    if plugin:
        cmd += ['--ro-bind', str(plugin), '/external_plugin.py']
    return cmd + ['/usr/bin/python3', '-s'] + argv

def bundle():
    dest = ROOT / 'runtime_public/gencad_bench_v2'
    dest.mkdir(parents=True, exist_ok=True)
    for name in ['__init__.py', 'common.py', 'geometry.py', 'tasks.py', 'generators.py', 'worker.py']:
        shutil.copy2(ROOT / 'gencad_bench_v2' / name, dest / name)
    (dest / 'vendor').mkdir(exist_ok=True)
    for name in ['continuous.py', 'obstacle.py', 'crown_fit_geometry.py']:
        shutil.copy2(ROOT / 'gencad_bench_v2/vendor' / name, dest / 'vendor' / name)

def run(participant, inputs=None, output=None, plugin=None):
    bundle()
    inputs = Path(inputs or DATA / 'public')
    output = Path(output or DATA / 'predictions' / participant)
    output.mkdir(parents=True, exist_ok=True)
    argv = ['-m', 'gencad_bench_v2.worker', '--participant', participant, '--inputs', '/inputs', '--output', '/output']
    if plugin:
        argv += ['--plugin', '/external_plugin.py:' + plugin[1]]
    cmd = command(inputs, output, argv, Path(plugin[0]) if plugin else None)
    tic = time.perf_counter()
    log = ROOT / 'raw' / ('participant_' + participant + '.log')
    with log.open('w') as f:
        p = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, timeout=1800)
    record = dict(argv=cmd, returncode=p.returncode, seconds=time.perf_counter() - tic, log=str(log))
    with (ROOT / 'raw/commands.jsonl').open('a') as f:
        f.write(json.dumps(record) + '\n')
    if p.returncode:
        raise RuntimeError('Participant process failed; see ' + str(log))
    for p in output.iterdir():
        if p.is_symlink() or not p.is_file() or p.stat().st_size > 2000000:
            raise ValueError('unsafe output file ' + str(p))
    return record
