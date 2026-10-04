import json, time
import numpy as np
from contact_model import P, write, sha, state
from run_r3 import continuous

def main():
    rows = []
    tick = time.perf_counter()
    for case in [1, 2, 3]:
        g = json.loads((P / 'raw' / f'geometry_{case:03d}_h02.json').read_text())
        (new, info) = continuous(g)
        old = json.loads((P / 'raw' / f'geometry_{case:03d}_continuous.json').read_text())
        error = float(np.max(np.abs(np.array([p['gap_mm'] for p in new['patches']]) - np.array([p['gap_mm'] for p in old['patches']]))))
        rows.append(dict(case=case, controls=info['controls'], prediction_identity_error_mm=error, cost=info['wall_s'], all_pass=bool(all((x['pass_gate'] for x in info['controls'])) and error == 0.0)))
        print('R3B', case, 'controls', len(info['controls']), 'identical gap', error, 'pass', rows[-1]['all_pass'], flush=True)
    write(P / 'rounds/R3B/results.json', dict(claim_type='capability', rows=rows, all_pass=all((x['all_pass'] for x in rows)), wall_s=time.perf_counter() - tick, prereg_sha256=sha(P / 'PREREG_R3B_CONTROL_REPAIR.json'), original_failure='rounds/R3/results.json all_LP_controls_passFalse; NaN control selection retained'))
if __name__ == '__main__':
    main()
