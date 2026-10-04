import copy
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
ROOT = Path(__file__).resolve().parent

def compare(got, expected, tol=1e-06):
    if isinstance(got, bool) or isinstance(expected, bool):
        return got is expected
    if isinstance(got, (int, float)) and isinstance(expected, (int, float)):
        return abs(got - expected) <= tol * max(1, abs(expected))
    if isinstance(got, dict) and isinstance(expected, dict):
        return got.keys() == expected.keys() and all((compare(got[k], expected[k], tol) for k in got))
    if isinstance(got, list) and isinstance(expected, list):
        return len(got) == len(expected) and all((compare(g, e, tol) for (g, e) in zip(got, expected)))
    return got == expected

def main():
    rows = []
    for x in json.loads((ROOT / 'LEGACY_SOURCE_MANIFEST.json').read_text()):
        assert hashlib.sha256((ROOT / x['local']).read_bytes()).hexdigest() == x['sha256']
    for name in ['adjres_end_to_end_demo', 'decidability_abstention_atlas', 'decidability_budget_aware_trigger', 'margin_tolerance_decidability', 'z24_sigma_governor_flagship', 'atlas_v2_ieff_reachability']:
        start = time.perf_counter()
        p = subprocess.run([os.environ.get('PYTHON', 'python3'), 'scripts/' + name + '.py'], cwd=ROOT / 'legacy', capture_output=True, text=True)
        (ROOT / 'raw' / ('replay_' + name + '.log')).write_text(p.stdout + '\nSTDERR\n' + p.stderr)
        report = ROOT / 'legacy/reports' / (name + '.json')
        expected = ROOT / 'legacy/expected' / (name + '.json')
        gates = json.loads(report.read_text())['gates'] if report.exists() else None
        exp = json.loads(expected.read_text())['gates']
        rows.append({'name': name, 'exit_code': p.returncode, 'seconds': time.perf_counter() - start, 'gates': gates, 'exact_gate_match': gates == exp, 'within_frozen_tolerance': compare(gates, exp)})
        print('legacy', name, 'exit', p.returncode, 'match', compare(gates, exp), flush=True)
    (ROOT / 'LEGACY_REPRODUCTION.json').write_text(json.dumps(rows, indent=2) + '\n')
    assert all((r['exit_code'] == 0 and r['within_frozen_tolerance'] for r in rows))
    atlas = json.loads((ROOT / 'legacy/reports/decidability_abstention_atlas.json').read_text())
    assert atlas['gates']['G3_naive_misfire_count'] == 3
    assert not compare(dict(atlas['gates'], G3_naive_misfire_count=2), atlas['gates'])
    fault = ROOT / 'raw/legacy_source_fault'
    (fault / 'reports').mkdir(parents=True, exist_ok=True)
    for name in ['s2_modal_resolution_frontier', 'z24_damage_decidability_governor', 'fatigue_resolution_frontier_ct', 'lpbf_hagb_veto_decidability_governor', 'stent_slice1_resolution_budget']:
        shutil.copyfile(ROOT / 'legacy/reports' / (name + '.json'), fault / 'reports' / (name + '.json'))
    fn = fault / 'reports/s2_modal_resolution_frontier.json'
    d = json.loads(fn.read_text())
    d['gates']['verdict'] = 'FAIL'
    fn.write_text(json.dumps(d, indent=2) + '\n')
    p = subprocess.run([os.environ.get('PYTHON', 'python3'), str(ROOT / 'legacy/scripts/decidability_abstention_atlas.py')], cwd=fault, capture_output=True, text=True)
    (fault / 'run.log').write_text(p.stdout + '\nSTDERR\n' + p.stderr)
    outcome = json.loads((fault / 'reports/decidability_abstention_atlas.json').read_text())
    receipts = {'changed_count_rejected': not compare(dict(atlas['gates'], G3_naive_misfire_count=2), atlas['gates']), 'failed_source_rejected': p.returncode == 1 and outcome['gates']['verdict'] == 'FAIL', 'injected_source_result': outcome['gates']}
    (ROOT / 'LEGACY_MUTATION_CHECK.json').write_text(json.dumps(receipts, indent=2) + '\n')
    assert receipts['changed_count_rejected'] and receipts['failed_source_rejected']
if __name__ == '__main__':
    main()
