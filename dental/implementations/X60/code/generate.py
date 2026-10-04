from common import *
import shutil, subprocess, time

def run():
    for (p, item) in read(ROOT / 'INPUT_LOCK.json')['files'].items():
        if sha(p) != item['sha256']:
            raise ValueError('input hash drift ' + p)
    public = DATA / 'plugin'
    public.mkdir(exist_ok=True)
    for name in ['plugin.py', 'plugin_worker.py', 'mechanics.py']:
        shutil.copy2(ROOT / 'code' / name, public / name)
    shutil.copy2(ROOT / 'PANEL_SELECTION.json', public / 'PANEL_SELECTION.json')
    freeze(ROOT / 'FROZEN_GENERATOR.json', dict(files=inventory(public), prereg_sha256=sha(ROOT / 'PREREG_R2.json'), source_checks_sha256=sha(V4 / 'payload/participant_code/legacy/checks.py')))
    for factor in [1.0, 1.1]:
        tag = 'volume_' + str(factor)
        out = DATA / 'predictions' / tag
        out.mkdir(parents=True, exist_ok=True)
        runtime = V4 / 'payload/runtime'
        cmd = [str(runtime / 'bin/bwrap'), '--die-with-parent', '--unshare-all', '--new-session', '--clearenv', '--ro-bind', '/usr', '/usr', '--ro-bind', '/lib', '/lib', '--ro-bind', '/lib64', '/lib64', '--proc', '/proc', '--dev', '/dev', '--tmpfs', '/tmp', '--dir', '/work', '--chdir', '/work', '--ro-bind', str(runtime), '/runtime', '--ro-bind', str(V4 / 'payload/participant_code'), '/participants', '--ro-bind', str(V4 / 'payload/public'), '/inputs', '--ro-bind', str(public), '/plugin', '--bind', str(out), '/output', '--setenv', 'HOME', '/tmp', '--setenv', 'X60_VOLUME_BUDGET', str(factor), '--setenv', 'PYTHONDONTWRITEBYTECODE', '1']
        for n in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
            cmd += ['--setenv', n, '1']
        cmd += ['/runtime/bin/python3', '-s', '-E', '-B', '/plugin/plugin_worker.py']
        tic = time.perf_counter()
        with (ROOT / 'raw' / f'generate_{tag}.log').open('w') as f:
            r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
        dump(ROOT / 'raw' / f'PROCESS_{tag}.json', dict(argv=cmd, exit_code=r.returncode, seconds=time.perf_counter() - tic))
        if r.returncode:
            raise RuntimeError('preserved generator error ' + tag)
        freeze(ROOT / f'FROZEN_PREDICTIONS_{tag}.json', dict(files=inventory(out), generator_sha256=sha(ROOT / 'FROZEN_GENERATOR.json'), private_quality_queries_before_freeze=0, physical_measurements='NOT_RUN'))
        state('R2_PREDICTIONS_FROZEN_' + tag, 'All generated coordinates frozen before source comparison', 'Score unmodified GenCAD predicates and quality')
    budget()
if __name__ == '__main__':
    run()
