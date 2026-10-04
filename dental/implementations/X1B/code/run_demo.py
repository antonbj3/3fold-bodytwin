"""One-command deterministic rebuild of analysis, exports checks and fresh FE replay."""
from dental_release.paths import expand as _release_expand
import argparse, sys, os, platform, importlib, time, resource
from pathlib import Path
from common import R, read, dump, sha, state

def validate():
    lock = read('RUNTIME_LOCK.json')
    if platform.python_version() != lock['python_version'] or sha(Path(sys.executable).resolve()) != lock['interpreter_sha256']:
        raise ValueError('Wrong Python runtime; use the pinned interpreter or reproduce its lock explicitly')
    for (package, version) in lock['packages'].items():
        if importlib.import_module(package).__version__ != version:
            raise ValueError('Wrong package version: ' + package)
    if sha(R / 'INPUT_LOCK.json') != (R / 'INPUT_LOCK.sha256').read_text().strip():
        raise ValueError('Input manifest changed')
    for f in read('INPUT_LOCK.json')['files']:
        if sha(R / f['path']) != f['sha256']:
            raise ValueError('Input hash drift: ' + f['path'])
    for p in R.glob('PREREG_*.json'):
        if sha(p) != p.with_suffix('.sha256').read_text().strip():
            raise ValueError('Frozen prereg drift')
    return len(read('INPUT_LOCK.json')['files'])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--skip-fresh-solve', action='store_true', help='analysis-only replay; explicitly labelled in output')
    ap.add_argument('--work-dir', default=os.environ.get('CROWN_WORKDATA', _release_expand('@DENTAL_WORK_ROOT@/LANE_X1B_CROWN_LOOP/fresh_solver')))
    a = ap.parse_args()
    start = time.monotonic()
    n = validate()
    state('PORTABLE_INPUTS_VERIFIED', dict(input_files=n), 'Recompute all 90 raw responses and physical export surfaces')
    from extract_literature import extract
    from calibrate import calibrate
    from matched_contrasts import run as matched
    from freeze_predictions import run as freeze
    from verify_demo import verify
    extract()
    calibrate()
    matched()
    freeze()
    check = verify()
    if not a.skip_fresh_solve:
        state('GATES_VERIFIED', dict(checks=len(check['checks']), passed=True), 'Fresh vendored CPU solve of two D1 load cases')
        from fresh_solve import solve
        fresh = solve(a.work_dir)
    else:
        fresh = dict(fresh_solve=False, explicit_analysis_only=True)
    runtime = dict(seconds=time.monotonic() - start, python_max_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024, threads_max=4, input_files=n, fresh_solve=fresh)
    from report import write
    write(check, runtime)
    print('X1b computation and fault checks PASS; absolute fracture transfer UNKNOWN. Read README_DEMO.md and crown_loop.png.')
if __name__ == '__main__':
    main()
