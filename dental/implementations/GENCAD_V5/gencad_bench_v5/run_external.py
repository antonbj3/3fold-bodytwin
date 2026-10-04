from common import *
import subprocess, shutil

def run():
    out = DATA / 'external_predictions'
    out.mkdir(exist_ok=True)
    cmd = ['bwrap', '--die-with-parent', '--unshare-all', '--new-session', '--clearenv', '--ro-bind', '/usr', '/usr', '--ro-bind', '/lib', '/lib', '--ro-bind', '/lib64', '/lib64', '--ro-bind', '/etc/ld.so.cache', '/etc/ld.so.cache', '--proc', '/proc', '--dev-bind', '/dev', '/dev', '--tmpfs', '/tmp', '--dir', '/work', '--chdir', '/work']
    venv = str(Path(PYTHON).parent.parent)
    for (a, b) in [(venv, venv), (str(DATA / 'deps'), '/deps'), (str(DATA / 'ToothCraft'), '/model'), (str(DATA / 'weights'), '/weights'), (str(DATA / 'external_inputs'), '/inputs'), (str(ROOT / 'gencad_bench_v5'), '/runner')]:
        cmd += ['--ro-bind', a, b]
    cmd += ['--bind', str(out), '/output', '--setenv', 'PYTHONPATH', '/deps:/model', '--setenv', 'PYTHONDONTWRITEBYTECODE', '1']
    for n in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
        cmd += ['--setenv', n, '4']
    cmd += [PYTHON, '-s', '-B', '/runner/external_worker.py']
    lock = {'prereg_sha256': sha(ROOT / 'PREREG_V5.json'), 'inputs_sha256': sha(ROOT / 'EXTERNAL_INPUTS.json'), 'model_lock_sha256': sha(ROOT / 'raw/EXTERNAL_MODEL_LOCK.json'), 'plugin_sha256': sha(ROOT / 'gencad_bench_v5/toothcraft_plugin.py'), 'worker_sha256': sha(ROOT / 'gencad_bench_v5/external_worker.py'), 'argv': cmd}
    if not (ROOT / 'EXTERNAL_GENERATOR_R3.json').exists():
        freeze(ROOT / 'EXTERNAL_GENERATOR_R3.json', lock)
    tic = time.perf_counter()
    with (ROOT / 'raw/external_generation.log').open('w') as log:
        p = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    dump(ROOT / 'raw/EXTERNAL_PROCESS.json', dict(argv=cmd, seconds=time.perf_counter() - tic, returncode=p.returncode))
    if p.returncode:
        raise RuntimeError('External generation failed: see raw/external_generation.log')
    freeze(ROOT / 'FROZEN_PREDICTIONS.json', dict(model='ToothCraft normal author checkpoints', generator_sha256=sha(ROOT / 'EXTERNAL_GENERATOR_R3.json'), files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir()}, physical_measurement='NOT_RUN', reference_scoring_before_freeze=False))
    state('EXTERNAL_FROZEN', 'Actual author network outputs frozen', 'Functional scoring and physical support tests')
if __name__ == '__main__':
    run()
