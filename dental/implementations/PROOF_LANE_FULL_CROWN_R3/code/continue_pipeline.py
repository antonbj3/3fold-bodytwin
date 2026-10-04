from common_r3 import *
import subprocess, time
while not (ROOT / 'FROZEN_PREDICTIONS_A.json').exists():
    p = subprocess.run(['pgrep', '-f', 'python code/resume_a.py'], capture_output=True, text=True)
    if p.returncode:
        raise RuntimeError('A resume launcher ended without prediction freeze; inspect raw/resume_A.log')
    time.sleep(5)
from score_round import run as score
score('A')
state('A_DECIDED_B_READY', 'A scored; all failures and handoff preserved', 'Execute frozen B explicit-boundary construction')
from run_b import run as generate_b
generate_b()
score('B')
from run_c import run as generate_c
generate_c()
score('C')
subprocess.run(['/usr/bin/time', '-v', '-o', str(ROOT / 'raw/SCORER_CONTROLS_MEMORY.txt'), PYTHON, str(ROOT / 'code/scorer_controls.py')], check=True)
from lab_export_r3 import run as export
export()
from report_r3 import run as report
report()
state('A_B_C_REPORTED', 'Three frozen constructions scored and exports selected', 'Run full one-command replay and bind review snapshot')
