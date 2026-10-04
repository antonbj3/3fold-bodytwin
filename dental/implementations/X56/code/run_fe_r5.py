"""Separate initial-state intervention; reused own fixed-contact K27 module."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
import json, sys, hashlib, importlib.util, datetime
LANE = Path(__file__).resolve().parents[1]
SRC = Path(_release_expand('@DENTAL_INPUT_ROOT@/workspace/cells/physics'))
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X56-micromotion/FE_R5'))

def main():
    sys.path.insert(0, str(SRC))
    own = LANE / 'code/k27_fixed_contact.py'
    spec = importlib.util.spec_from_file_location('k27fixed', own)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    cfgs = json.loads((LANE / 'raw/R4_FE_specs.json').read_text())
    base = cfgs[-1]
    output = LANE / 'raw/R5_FE_runs.jsonl'
    done = set()
    if output.exists():
        done = {json.loads(x)['id'] for x in output.read_text().splitlines()}
    cases = []
    for fac in [0.1, 0]:
        d = {**base, 'id': 'R5_prestrain' + str(fac), 'T_c': base['T_c'] * fac, 'T_s': base['T_s'] * fac}
        cases.append(d)
    (LANE / 'raw/R5_FE_specs.json').write_text(json.dumps(cases, indent=2) + '\n')
    for d in cases:
        if d['id'] in done:
            continue
        print('START', d['id'], flush=True)
        try:
            r = mod.solve(d, str(DATA / d['id']), tag='mm')
        except Exception:
            import traceback
            r = {'ok': False, 'fail': traceback.format_exc(), 'spec': d}
        r['id'] = d['id']
        r['code_sha256'] = hashlib.sha256(own.read_bytes()).hexdigest()
        with output.open('a') as f:
            f.write(json.dumps(r, default=float) + '\n')
        print('RESULT', d['id'], r.get('ok'), r.get('u_load_along_F_um'), r.get('slip_max_um'), r.get('t_total'), flush=True)
        state = json.loads((LANE / 'CURRENT_WORK_STATE.json').read_text())
        state['updated_at'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        state['last_FE_run'] = {'id': d['id'], 'ok': r.get('ok'), 'u_load_along_F_um': r.get('u_load_along_F_um')}
        (LANE / 'CURRENT_WORK_STATE.json').write_text(json.dumps(state, indent=2) + '\n')
if __name__ == '__main__':
    main()
