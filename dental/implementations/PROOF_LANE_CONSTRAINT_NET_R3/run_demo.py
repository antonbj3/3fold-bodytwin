import os, subprocess, time, resource
from datetime import datetime, timezone
from common import *

def run():
    start = time.perf_counter()
    steps = []
    commands = [['python3', 'prepare_inputs.py'], ['python3', 'build_net.py'], ['python3', 'validate.py', '--integrity-only', '--mutations'], ['python3', 'lab_plan.py'], ['python3', 'stress_map.py'], ['python3', 'check_controls.py'], ['python3', 'make_figure.py']]
    for cmd in commands:
        t = time.perf_counter()
        p = subprocess.run(cmd, cwd=HERE, check=False, capture_output=True, text=True)
        steps.append({'argv': cmd, 'exit_code': p.returncode, 'wall_seconds': time.perf_counter() - t, 'stdout': p.stdout, 'stderr': p.stderr})
        if p.returncode:
            dump('RUN_FAILURE.json', {'utc': datetime.now(timezone.utc).isoformat(), 'steps': steps})
            raise SystemExit(p.returncode)
    v = read(HERE / 'VALIDATION.json')
    controls = read(HERE / 'CONTROLS.json')
    assert v['integrity_pass'] and controls['all_rejected']
    p = subprocess.run(['python3', 'validate.py'], cwd=HERE, capture_output=True, text=True)
    dump('STRICT_CLI_RECEIPT.json', {'argv': ['python3', 'validate.py'], 'exit_code': p.returncode, 'stdout': p.stdout, 'stderr': p.stderr, 'expected': 'Reject every remaining OPEN/UNKNOWN without numeric gap or physical obligation; five known nonphysical exceptions.'})
    dump('VALIDATION.json', v)
    if p.returncode != 1:
        raise RuntimeError('STRICT_GATE_DID_NOT_REJECT')
    cost = {'measured_at_utc': datetime.now(timezone.utc).isoformat(), 'software_wall_seconds_before_reporting': time.perf_counter() - start, 'peak_child_rss_kib': resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss, 'steps': steps, 'threads': int(os.environ['OMP_NUM_THREADS']), 'fit_seconds': 0, 'physical_acquisition': 'NOT_RUN', 'discovery_manual_labour': 'UNKNOWN; PREREG precedes exact verification/counting', 'lab_time': 'X71 intervals are PHENOMENOLOGICAL planning assumptions; all new extension time/cost UNKNOWN', 'validation': 'source identity, net mutations, proposal mutations, join mutations and direct count control included', 'queries': 'finite set groups included in lab_plan step', 'fallback': 'none; five nonphysical exceptions retained', 'resolution_level': 'PHENOMENOLOGICAL; timings are this software run only'}
    dump('COST_RUN.json', cost)
    from report import main
    main()
    dump('CURRENT_WORK_STATE.json', {'milestone': 'R3B_COMPLETE', 'current': 'Delivery and independent-review handoff', 'latest_gate': 'integrity PASS; 75/75 fault controls reject; strict availability FAIL for five nonphysical exceptions', 'next_operation': 'Independent review of measurement contracts; obtain matched pilot data and code/reference/rights evidence for five exceptions. No new physical observations available in this run.', 'review_state': 'PENDING_INDEPENDENT_REVIEW'})
    print(json.dumps({'integrity': v['integrity_pass'], 'strict_availability': v['availability_contract_complete'], 'numeric': v['numeric_open_unknown_edges'], 'physical_blocked': v['physically_blocked_edges'], 'nonphysical': v['unhandled_open_unknown_edges'], 'controls': controls['count']}))
    return 1
if __name__ == '__main__':
    raise SystemExit(run())
