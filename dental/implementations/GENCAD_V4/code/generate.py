from common import *
import subprocess, shutil

def command(output, worker='/runner/worker.py'):
    cmd = [str(DATA / 'runtime/bin/bwrap'), '--die-with-parent', '--unshare-all', '--new-session', '--clearenv', '--ro-bind', '/usr', '/usr', '--ro-bind', '/lib', '/lib', '--ro-bind', '/lib64', '/lib64', '--proc', '/proc', '--dev', '/dev', '--tmpfs', '/tmp', '--dir', '/work', '--chdir', '/work', '--ro-bind', str(DATA / 'runtime'), '/runtime', '--ro-bind', str(DATA / 'participant_code'), '/participants', '--ro-bind', str(DATA / 'public'), '/inputs', '--ro-bind', str(DATA / 'models'), '/models', '--ro-bind', str(DATA / 'runner'), '/runner', '--bind', str(output), '/output', '--setenv', 'HOME', '/tmp', '--setenv', 'PYTHONDONTWRITEBYTECODE', '1', '--setenv', 'PYTHONWARNINGS', 'ignore', '--setenv', 'FIELD_MESH_SDF_LIBRARY', '/participants/field/native_mesh_sdf.so']
    for n in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'THREADS', 'FIELD_ENGINE_NATIVE_THREADS']:
        cmd += ['--setenv', n, '4' if n == 'FIELD_ENGINE_NATIVE_THREADS' else '1']
    return cmd + ['/runtime/bin/python3', '-s', '-E', '-B', worker]

def run():
    start = time.perf_counter()
    out = DATA / 'predictions'
    out.mkdir(exist_ok=True)
    (DATA / 'runner').mkdir(exist_ok=True)
    shutil.copy2(ROOT / 'code/worker.py', DATA / 'runner/worker.py')
    shutil.copy2(ROOT / 'code/field_worker.py', DATA / 'runner/field_worker.py')
    inputs = {}
    for d in ['public', 'participant_code', 'models', 'runner', 'runtime']:
        inputs.update({d + '/' + k: v for (k, v) in inventory(DATA / d).items()})
    freeze(ROOT / 'FROZEN_GENERATOR.json', dict(prereg_sha256=sha(ROOT / 'PREREG_R2.json'), files=inputs, scope='Public inputs and executable bytes bound before generation; private references not mounted'))
    cmd = command(out)
    with (ROOT / 'raw/generation.log').open('w') as f:
        r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    dump(ROOT / 'raw/GENERATION_PROCESS.json', dict(argv=cmd, exit_code=r.returncode, seconds=time.perf_counter() - start))
    if r.returncode:
        raise RuntimeError('generation failed; preserve log')
    freeze(ROOT / 'FROZEN_PREDICTIONS.json', dict(generator_sha256=sha(ROOT / 'FROZEN_GENERATOR.json'), files=inventory(out), physical_measurement='NOT_RUN', quality_reference_queries_before_freeze=0))
    state('PREDICTIONS_FROZEN', 'All participant outputs frozen before quality scores', 'Evaluate reference-owned spatial metrics')
if __name__ == '__main__':
    run()
