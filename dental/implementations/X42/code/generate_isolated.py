from common import *
import subprocess

def run():
    out = DATA / 'PREDICTIONS'
    out.mkdir(exist_ok=True)
    b = BENCH / 'payload'
    cmd = [str(b / 'runtime/bin/bwrap'), '--die-with-parent', '--unshare-all', '--new-session', '--clearenv', '--ro-bind', '/usr', '/usr', '--ro-bind', '/lib', '/lib', '--ro-bind', '/lib64', '/lib64', '--proc', '/proc', '--dev', '/dev', '--tmpfs', '/tmp', '--dir', '/work', '--chdir', '/work', '--ro-bind', str(b / 'runtime'), '/runtime', '--ro-bind', str(b / 'participant_code'), '/baseline', '--dir', '/candidate', '--ro-bind', str(ROOT / 'code/participant.py'), '/candidate/participant.py', '--ro-bind', str(ROOT / 'code/generator.py'), '/candidate/generator.py', '--ro-bind', str(b / 'public'), '/inputs', '--ro-bind', str(DATA / 'FINAL_MODELS'), '/models', '--ro-bind', str(ROOT / 'FINAL_CONFIG.json'), '/config.json', '--bind', str(out), '/output', '--setenv', 'HOME', '/tmp', '--setenv', 'PYTHONDONTWRITEBYTECODE', '1', '--setenv', 'PYTHONWARNINGS', 'ignore']
    for var in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
        cmd += ['--setenv', var, '1']
    cmd += ['/runtime/bin/python3', '-s', '-E', '-B', '/candidate/participant.py', '--inputs', '/inputs', '--models', '/models', '--output', '/output', '--config', '/config.json']
    tic = time.perf_counter()
    with (ROOT / 'raw/generation.log').open('w') as log:
        p = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    dump(ROOT / 'raw/ISOLATION_RECEIPT.json', dict(argv=cmd, exit_code=p.returncode, seconds=time.perf_counter() - tic, reference_mount=False, scorer_mount=False, network='unshare-all', original_checkout_writes=False))
    if p.returncode:
        raise RuntimeError('isolated generator failed; preserve outputs/log before any repair')
if __name__ == '__main__':
    run()
