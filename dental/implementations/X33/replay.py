"""One-command replay without rewriting the frozen first-run evidence."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
import subprocess, os, datetime, json, shutil, time
import numpy as np
ROOT = Path(__file__).resolve().parent

def main():
    start = time.time()
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out = ROOT / 'replay' / stamp
    out.mkdir(parents=True)
    for d in ['raw', 'figures', 'inputs']:
        (out / d).mkdir()
    for p in ROOT.glob('PREREG_R*'):
        if p.is_file():
            shutil.copyfile(p, out / p.name)
    for p in (ROOT / 'inputs').glob('*'):
        if p.is_file():
            shutil.copyfile(p, out / 'inputs' / p.name)
    env = os.environ.copy()
    env['X33_RUN_ROOT'] = str(out)
    env['X33_MATRIX_DIR'] = _release_expand('@DENTAL_WORK_ROOT@/X33/replay/') + stamp
    for script in ['laser_heat.py', 'inverse_measurement.py', 'heat_bin_budget.py']:
        with (out / 'raw' / (script + '.log')).open('w') as log:
            subprocess.run(['python3', str(ROOT / script)], env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    checks = {}
    for (name, keys) in [('R1', ['last_refinement_change_C']), ('R2', ['uniform_peak_C', 'minimum_possible_maximum_C', 'maximum_possible_maximum_C']), ('R3', ['coarse_acceptances', 'coarse_acceptances_refuted_by_finer_feasible_history'])]:
        reference = json.loads((ROOT / f'raw/{name}_RESULTS.json').read_text())
        got = json.loads((out / f'raw/{name}_RESULTS.json').read_text())
        checks[name] = {key: bool(np.isclose(reference[key], got[key], rtol=1e-10, atol=1e-07)) for key in keys}
    assert all((all(v.values()) for v in checks.values()))
    subprocess.run(['/usr/bin/python3', '-s', str(ROOT / 'report.py')], env=env, check=True)
    result = {'status': 'PASS_REPLAY', 'scientific_values_match': checks, 'wall_s': time.time() - start, 'output': str(out), 'original_results_unchanged': True}
    (out / 'REPLAY_VALIDATION.json').write_text(json.dumps(result, indent=2) + '\n')
    (ROOT / 'LATEST_REPLAY.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)
if __name__ == '__main__':
    main()
