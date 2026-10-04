from common import *
from generate import command
import shutil, subprocess

def run():
    out = DATA / 'prep_predictions'
    out.mkdir(exist_ok=True)
    runner = DATA / 'prep_runner'
    runner.mkdir(exist_ok=True)
    shutil.copy2(ROOT / 'code/prep_worker.py', runner / 'prep_worker.py')
    shutil.copy2(ROOT / 'PREREG_R5.json', runner / 'PREP_CONFIG.json')
    freeze(ROOT / 'FROZEN_PREP_GENERATOR.json', dict(runner_files=inventory(runner), input_freeze_sha256=sha(ROOT / 'FROZEN_WHOLE_INPUTS.json')))
    cmd = command(out, '/runner/prep_worker.py')
    cmd[cmd.index(str(DATA / 'runner'))] = str(runner)
    i = cmd.index('/runtime/bin/python3')
    cmd[i:i] = ['--ro-bind', str(DATA / 'whole_inputs'), '/whole_inputs']
    start = time.perf_counter()
    with (ROOT / 'raw/prep_generation.log').open('w') as f:
        r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    dump(ROOT / 'raw/PREP_GENERATION_PROCESS.json', dict(argv=cmd, exit_code=r.returncode, seconds=time.perf_counter() - start))
    if r.returncode:
        raise RuntimeError('Preparation-conditioned generation failed')
    freeze(ROOT / 'FROZEN_PREP_PREDICTIONS.json', dict(files=inventory(out), generator_sha256=sha(ROOT / 'FROZEN_PREP_GENERATOR.json'), new_track_reference_queries_before_freeze=0, retrospective_source_panel=True))
if __name__ == '__main__':
    run()
