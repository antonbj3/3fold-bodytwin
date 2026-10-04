from common import *
import shutil, subprocess
from generate import command

def run():
    start = time.perf_counter()
    out = DATA / 'whole_predictions'
    out.mkdir(exist_ok=True)
    runner = DATA / 'whole_runner'
    runner.mkdir(exist_ok=True)
    shutil.copy2(ROOT / 'code/whole_worker.py', runner / 'whole_worker.py')
    shutil.copy2(ROOT / 'PREREG_R3.json', runner / 'WHOLE_CONFIG.json')
    cmd = command(out, worker='/runner/whole_worker.py')
    idx = cmd.index(str(DATA / 'runner'))
    cmd[idx] = str(runner)
    insertion = cmd.index('/runtime/bin/python3')
    cmd[insertion:insertion] = ['--ro-bind', str(DATA / 'whole_inputs'), '/whole_inputs', '--ro-bind', str(DATA / 'predictions'), '/roofs']
    freeze(ROOT / 'FROZEN_WHOLE_GENERATOR.json', dict(runner_files=inventory(runner), input_freeze_sha256=sha(ROOT / 'FROZEN_WHOLE_INPUTS.json'), roof_freeze_sha256=sha(ROOT / 'FROZEN_PREDICTIONS.json')))
    with (ROOT / 'raw/whole_generation.log').open('w') as f:
        r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    dump(ROOT / 'raw/WHOLE_GENERATION_PROCESS.json', dict(argv=cmd, exit_code=r.returncode, seconds=time.perf_counter() - start))
    if r.returncode:
        raise RuntimeError('Whole generation failed')
    freeze(ROOT / 'FROZEN_WHOLE_PREDICTIONS.json', dict(files=inventory(out), generator_sha256=sha(ROOT / 'FROZEN_WHOLE_GENERATOR.json'), full_surface_reference_queries_before_freeze=0, physical_measurement='NOT_RUN'))
    state('WHOLE_CROWNS_FROZEN', 'Full shells generated without reference mount', 'Score axial and occlusal exterior against untouched target and neighbors')
if __name__ == '__main__':
    run()
