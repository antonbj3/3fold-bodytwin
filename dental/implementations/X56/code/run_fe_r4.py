"""Reuse K27 read-only, with own fixed per-layer contact closure inputs."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
import json, importlib.util, sys, hashlib
LANE = Path(__file__).resolve().parents[1]
SOURCE = Path(_release_expand('@DENTAL_INPUT_ROOT@/workspace/cells/physics/micromotion_contact_fe.py'))
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X56-micromotion/FE_R4'))

def main():
    DATA.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(SOURCE.parent))
    src = SOURCE.read_text()
    old = 'Kd = {"CORT": p["kfac"] * p["E_c"] / p["h_if"], "CAN": p["kfac"] * p["E_s"] / p["h_if"]}'
    new = 'Kd = dict(p["K_pen"]) if isinstance(p["K_pen"], dict) else {"CORT": p["kfac"] * p["E_c"] / p["h_if"], "CAN": p["kfac"] * p["E_s"] / p["h_if"]}'
    old2 = 'lamd = {g: p["lamfac"] * K / 10.0 for g, K in Kd.items()}'
    new2 = 'lamd = dict(p["lam"]) if isinstance(p["lam"], dict) else {g: p["lamfac"] * K / 10.0 for g, K in Kd.items()}'
    assert src.count(old) == 1 and src.count(old2) == 1, 'Source assignments changed; inspect'
    own = LANE / 'code/k27_fixed_contact.py'
    own.write_text(src.replace(old, new).replace(old2, new2))
    modspec = importlib.util.spec_from_file_location('k27fixed', own)
    mod = importlib.util.module_from_spec(modspec)
    modspec.loader.exec_module(mod)
    default = json.loads((LANE / 'raw/historical_runs.json').read_text())[1]['spec']
    cfg = json.loads((LANE / 'PREREG_R4.json').read_text())['parameters']
    record = LANE / 'raw/R4_FE_runs.jsonl'
    done = set()
    if record.exists():
        done = {json.loads(x)['id'] for x in record.read_text().splitlines()}
    names = ['default', 'foam_only', 'cortex_only', 'both', 'both_mu06']
    specs = []
    for (name, (ec, es), mu) in zip(names, cfg['material_scenarios'], cfg['frictions']):
        p = {**mod.DEFAULT, **default, 'id': 'R4_' + name, 'E_c': ec, 'E_s': es, 'mu': mu, 'h_if': 0.4, 'K_pen': dict(zip(['CORT', 'CAN'], cfg['penalty_fixed_CORT_CAN_MPa_per_mm'])), 'lam': dict(zip(['CORT', 'CAN'], cfg['stick_fixed_CORT_CAN'])), 'T_c': cfg['expansion_temperatures_c_s'][0], 'T_s': cfg['expansion_temperatures_c_s'][1], 'calib': False}
        specs.append(p)
    (LANE / 'raw/R4_FE_specs.json').write_text(json.dumps(specs, indent=2) + '\n')
    for p in specs:
        if p['id'] in done:
            continue
        print('START', p['id'], flush=True)
        try:
            r = mod.solve(p, str(DATA / p['id']), tag='mm')
        except Exception:
            import traceback
            r = {'ok': False, 'fail': traceback.format_exc(), 'spec': p}
        r['id'] = p['id']
        r['source_sha256'] = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
        with record.open('a') as f:
            f.write(json.dumps(r, default=float) + '\n')
        print('RESULT', p['id'], r.get('ok'), r.get('u_load_along_F_um'), r.get('slip_max_um'), r.get('t_total'), flush=True)
        s = json.loads((LANE / 'CURRENT_WORK_STATE.json').read_text())
        s['last_FE_run'] = {'id': p['id'], 'ok': r.get('ok'), 'u_load_along_F_um': r.get('u_load_along_F_um')}
        (LANE / 'CURRENT_WORK_STATE.json').write_text(json.dumps(s, indent=2) + '\n')
if __name__ == '__main__':
    main()
