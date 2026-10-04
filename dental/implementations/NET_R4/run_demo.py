import time, resource, copy, subprocess, sys
from common import *
from extract_sources import extract
from gate_sources import run as source_gate
from build_net import build
from validate import validate, mutations
from range_control import run as range_control
from sufficiency import run as sufficiency
from fatigue_countercheck import run as fatigue_countercheck
from plans import build as plans
from make_figure import make
from report import build as report
from export_bindings import build as bindings

def main():
    start = time.perf_counter()
    verify()
    state('REPLAY_RUNNING', 'Frozen source hashes verified', 'Replay primary source and range gates')
    extract()
    a = source_gate()
    n = build()
    v = validate(n)
    v['fault_controls'] = mutations(n)
    write('VALIDATION.json', v)
    r = range_control(n)
    s = sufficiency()
    f = fatigue_countercheck()
    plans(n)
    make()
    (H / 'examples').mkdir(exist_ok=True)
    for (name, i, val) in [('valid_preload', 0, None), ('invalid_preload', 0, 80100), ('invalid_internal_gap', 3, 801000), ('valid_internal_gap', 3, None)]:
        o = copy.deepcopy(n['r4_observations'][i])
        if val is not None:
            o['value'] = val
        write('examples/' + name + '.json', o)
    cli = []
    for (name, expected) in [('valid_preload', 0), ('invalid_preload', 1), ('invalid_internal_gap', 1), ('valid_internal_gap', 0)]:
        cmd = [sys.executable, str(H / 'validate.py'), '--observation', 'examples/' + name + '.json']
        p = subprocess.run(cmd, cwd=H, text=True, capture_output=True)
        cli.append({'command': cmd, 'exit_code': p.returncode, 'expected_exit_code': expected, 'stdout': p.stdout, 'stderr': p.stderr, 'pass': p.returncode == expected and (expected == 0 or 'OUT_OF_PLAUSIBLE_RANGE' in p.stdout)})
    p = subprocess.run([sys.executable, str(H / 'validate.py'), '--strict'], cwd=H, text=True, capture_output=True)
    write('STRICT_CLI_RECEIPT.json', {'exit_code': p.returncode, 'expected_exit_code': 1, 'stdout': p.stdout, 'stderr': p.stderr, 'meaning': 'Full chain remains unvalidated; expected negative gate.'})
    write('CLI_TESTS.json', cli)
    passed = a['pass'] and v['integrity_pass'] and v['fault_controls']['all_rejected'] and r['all_pass'] and s['all_pass'] and f['all_pass'] and all((x['pass'] for x in cli)) and (p.returncode == 1)
    runtime = {'wall_seconds_before_report': time.perf_counter() - start, 'peak_RSS_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'threads': 1, 'gpu': False, 'input_snapshot_bytes': sum((e['bytes'] for e in entries().values())), 'preparation_manual_hours': 'UNKNOWN', 'fit_seconds': 0, 'discovery_total_manual_and_agent_cost': 'UNKNOWN; not reduced to replay runtime', 'validation': 'Source replay + interval membership + exact witness + actual mutations + public CLI', 'physical_measurements': 0, 'fallback': 'ABSTAIN for unsupported scopes/ranges; no imputation', 'run_utc': utc()}
    write('COST_RUN.json', runtime)
    result = report()
    bindings()
    manifest = {str(p.relative_to(H)): sha(p) for p in H.rglob('*') if p.is_file() and p.name not in ['ARTIFACT_MANIFEST.json', 'CURRENT_WORK_STATE.json'] and ('__pycache__' not in str(p)) and ('.mpl' not in str(p))}
    write('ARTIFACT_MANIFEST.json', manifest)
    state('R4_ROUND_COMPLETE_PENDING_INDEPENDENT_REVIEW' if passed else 'R4_GATE_FAILED', 'Source/range gates pass; strict full-chain gate fails as expected' if passed else 'Inspect VALIDATION/CLI_TESTS', 'Independent review, then matched region/state pilot R4C or source acquisition for67 scalar range debts')
    print(json.dumps({'demo_pass': passed, 'summary': result['summary'], 'wall_seconds': time.perf_counter() - start, 'strict_expected_exit': p.returncode}))
    return 0 if passed else 1
if __name__ == '__main__':
    raise SystemExit(main())
