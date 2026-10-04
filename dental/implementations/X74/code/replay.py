"""One-command offline regeneration and scientific replay, with immutable first receipt."""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '4'
import json, subprocess, sys, time, shutil, datetime
from common import ROOT, dump, sha, stamp, check_frozen, state

def scientific(x):
    x = dict(x)
    x.pop('cost', None)
    return x

def main():
    tic = time.perf_counter()
    first = ROOT / 'first_execution'
    first.mkdir(exist_ok=True)
    for r in ['R1', 'R2', 'R3', 'R4']:
        dest = first / (r + '.json')
        if not dest.exists():
            shutil.copyfile(ROOT / 'rounds' / (r + '.json'), dest)
    for f in ['PREREG_R1.json', 'PREREG_R2.json', 'PREREG_R3.json', 'PREREG_R4.json', 'SOURCE_MANIFEST.json', 'DECOMPOSITION.json', 'FROZEN_PREDICTIONS_R1.json', 'FROZEN_PREDICTIONS.json', 'FROZEN_LAB_PREDICTIONS.json']:
        check_frozen(f)
    if (ROOT / 'CODE_MANIFEST.json').exists():
        for (path, digest) in check_frozen('CODE_MANIFEST.json').items():
            assert sha(ROOT / path) == digest, 'Code drift: ' + path
    run = ROOT / 'runs' / datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    run.mkdir(parents=True)
    commands = []
    for script in ['run_r1.py', 'run_r2.py', 'run_r3_source_gate.py', 'run_r4.py', 'make_figure.py']:
        cmd = [sys.executable] + (['-s'] if script == 'make_figure.py' else []) + [str(ROOT / 'code' / script)]
        start = time.perf_counter()
        r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
        (run / (script + '.stdout.log')).write_text(r.stdout)
        (run / (script + '.stderr.log')).write_text(r.stderr)
        commands.append({'argv': cmd, 'cwd': str(ROOT), 'threads_env': {k: os.environ[k] for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']}, 'exit_code': r.returncode, 'seconds': time.perf_counter() - start})
        dump(run / 'commands.json', commands)
        assert r.returncode == 0, f'{script} failed; see {run}'
    for r in ['R1', 'R2', 'R3', 'R4']:
        old = json.loads((first / (r + '.json')).read_text())
        new = json.loads((ROOT / 'rounds' / (r + '.json')).read_text())
        assert scientific(old) == scientific(new), 'Scientific replay drift: ' + r
        shutil.copyfile(ROOT / 'rounds' / (r + '.json'), run / (r + '.json'))
    check_frozen('FROZEN_PREDICTIONS.json')
    check_frozen('FROZEN_LAB_PREDICTIONS.json')
    size = sum((p.stat().st_size for p in ROOT.rglob('*') if p.is_file()))
    assert size < 3000000000, 'Lane intermediate budget exceeded'
    receipt = {'pass': True, 'scientific_rounds_exactly_reproduced': ['R1', 'R2', 'R3', 'R4'], 'frozen_predictions_unchanged': True, 'commands': commands, 'elapsed_s': time.perf_counter() - tic, 'package_bytes': size, 'physical_validation': 'NOT_RUN', 'patient_transfer': 'UNKNOWN', 'review_state': 'PENDING_INDEPENDENT_REVIEW'}
    dump(run / 'VERIFICATION.json', receipt)
    dump(ROOT / 'VERIFICATION.json', receipt)
    state(status='COMPLETE_DELIVERY_PENDING_INDEPENDENT_REVIEW', latest_gate={'offline_replay': 'PASS', 'R1_transfer_rules': 'FAIL', 'R3_source_geometry': 'FAIL', 'patient_transfer': 'UNKNOWN'}, next_operation='Manufacture/measure source-matched coupons from LAB_REQUESTS.json, then same-batch curved crown with local force-displacement/colour maps; graph needs coordinator working-view refresh before receipts')
    print('PASS: source tables, local geometry/traction, nonlinear enclosures, all fault injections and retained failures replay exactly.')
    print('Patient physical transfer: UNKNOWN; future measurements: NOT_RUN. Figure: figures/restorative_demo.png')
if __name__ == '__main__':
    main()
