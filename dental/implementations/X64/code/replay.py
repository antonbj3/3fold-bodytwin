"""One-command rerun, preserving immutable first-run outcomes and predictions."""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '4'
import datetime, json, subprocess, sys, time
from pathlib import Path
from freeze import ROOT, sha, dump

def main():
    run = ROOT / 'runs' / datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    run.mkdir(parents=True)
    cmds = []
    for round_ in ('R1', 'R2'):
        out = run / round_
        out.mkdir()
        env = dict(os.environ, X64_REPLAY_DIR=str(out))
        cmd = [sys.executable, str(ROOT / 'code' / f'run_{round_.lower()}.py')]
        tic = time.perf_counter()
        r = subprocess.run(cmd, env=env, cwd=ROOT, capture_output=True, text=True)
        (out / 'stdout.log').write_text(r.stdout)
        (out / 'stderr.log').write_text(r.stderr)
        cmds.append({'argv': cmd, 'cwd': str(ROOT), 'env_overrides': {'X64_REPLAY_DIR': str(out), 'OMP_NUM_THREADS': '4', 'OPENBLAS_NUM_THREADS': '4', 'MKL_NUM_THREADS': '4'}, 'exit_code': r.returncode, 'seconds': time.perf_counter() - tic})
        dump(run / 'commands.json', cmds)
        assert r.returncode == 0, f'{round_} failed: {out}/stderr.log'
        original = json.loads((ROOT / 'raw' / f'RESULTS_{round_}.json').read_text())
        new = json.loads((out / 'raw' / f'RESULTS_{round_}.json').read_text())
        original.pop('cost')
        new.pop('cost')
        assert original == new, f'{round_} numerical/scientific outcome drift'
    cmd = [sys.executable, '-s', str(ROOT / 'code/make_figure.py')]
    tic = time.perf_counter()
    figure_env = dict(os.environ, X64_REPLAY_DIR=str(run / 'figure'))
    r = subprocess.run(cmd, cwd=ROOT, env=figure_env, capture_output=True, text=True)
    (run / 'figure_stdout.log').write_text(r.stdout)
    (run / 'figure_stderr.log').write_text(r.stderr)
    cmds.append({'argv': cmd, 'cwd': str(ROOT), 'exit_code': r.returncode, 'seconds': time.perf_counter() - tic})
    dump(run / 'commands.json', cmds)
    assert r.returncode == 0, r.stderr
    original_prediction_sha = (ROOT / 'FROZEN_PREDICTIONS.sha256').read_text().strip()
    assert sha(ROOT / 'FROZEN_PREDICTIONS.json') == original_prediction_sha
    assert sha(run / 'R2/FROZEN_PREDICTIONS.json') == original_prediction_sha
    verdict = {'pass': True, 'replay_directory': str(run), 'compared': 'Exact scientific outputs; timing/RSS excluded from equality, recorded separately', 'frozen_predictions_unchanged': True, 'commands': cmds}
    dump(run / 'VERIFICATION.json', verdict)
    dump(ROOT / 'VERIFICATION.json', verdict)
    print('PASS: both chains rerun, source/numeric controls and adverse injections pass; frozen predictions unchanged.')
    print('Physical D1 optical/fracture validation: UNKNOWN. Figure: manufacturing_demo.png')
if __name__ == '__main__':
    main()
