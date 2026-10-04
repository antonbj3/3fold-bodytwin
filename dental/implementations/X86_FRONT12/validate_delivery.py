"""Check review bindings, complete execution identities and falsifiable reader gates."""
import copy, json
from pathlib import Path
from prepare_replay import HERE, PACKAGE, sha, write, execution_manifest
from receipt_contract import validator_pass

def check_review(r):
    if not r.get('path') or not r.get('sha256') or (not Path(r['path']).is_file()):
        return False
    if sha(r['path']) != r['sha256']:
        return False
    j = json.loads(Path(r['path']).read_text())
    return j.get('result_sha256') == r.get('result_sha256') and r.get('result_sha256') is not None and (j.get('graph_decision') in ['ACCEPT', 'ACCEPT_WITH_CORRECTION']) and (j.get('demo_decision') in ['DEMO_READY', 'DEMO_READY_WITH_CORRECTION'])

def complete_coverage(records, expected):
    ids = [r['demo_id'] for r in records]
    return len(ids) == len(set(ids)) and set(ids) == set(expected) and all(('exit_code' in r and 'wall_seconds' in r and ('log' in r) for r in records))

def passed_run(r):
    receipt = r.get('validator_receipt')
    if r.get('exit_code') != 0 or receipt is None:
        return False
    return validator_pass(receipt, r['demo_id']) and (not r.get('mismatches'))

def validate(require_complete=True):
    reviews = json.loads((HERE / 'raw/REVIEW_DELTA.json').read_text())['records']
    eligible = [r for r in reviews if r['eligible']]
    assert len(reviews) == 8 and len(eligible) == 5 and all((check_review(r) for r in eligible))
    probes = []
    for r in eligible:
        for (field, value) in [('path', None), ('sha256', '0' * 64), ('result_sha256', '1' * 64)]:
            bad = copy.deepcopy(r)
            bad[field] = value
            probes.append({'control': 'review_' + field, 'injected_value': value, 'rejected': not check_review(bad)})
    expected = [r['id'] for r in json.loads(execution_manifest().read_text())['executions']]
    inventory = json.loads((HERE / 'raw/EXECUTABLE_INVENTORY_AUDIT_R3.json').read_text())
    assert inventory['pass'] and (not inventory['omitted']) and (not inventory['unexpected'])
    ledger = json.loads((HERE / 'raw/ALL_REPLAY.json').read_text())
    if require_complete:
        assert complete_coverage(ledger['records'], expected)
    for rec in ledger['records']:
        if rec['status'] == 'PASS':
            assert passed_run(rec)
            bad = copy.deepcopy(rec)
            bad['exit_code'] = 1
            probes.append({'control': 'exit_code_' + rec['demo_id'], 'injected_value': 1, 'rejected': not passed_run(bad)})
            bad = copy.deepcopy(rec)
            bad['mismatches'] = [{'expected': 8, 'actual': 999}]
            probes.append({'control': 'mismatch_' + rec['demo_id'], 'injected_value': 999, 'rejected': not passed_run(bad)})
    if ledger['records']:
        fake = copy.deepcopy(ledger['records'])
        fake[0]['demo_id'] = 'ABSENT_DEMO'
        probes.append({'control': 'executable_identity', 'injected_value': 'ABSENT_DEMO', 'rejected': not complete_coverage(fake, expected)})
    suff = json.loads((HERE / 'raw/SUFFICIENCY.json').read_text())
    assert suff['identity_error'] == suff['float64_identity_error'] == 0 and suff['downstream_difference'] == 1
    assert all((p['rejected'] for p in probes))
    report = {'review_bindings_pass': True, 'eligible_review_records': len(eligible), 'rejected_review_records': 3, 'full_identity_coverage': complete_coverage(ledger['records'], expected), 'expected_executables': len(expected), 'scientific_demos': 96, 'historical_administrative_controls': 1, 'completed_demos': len(ledger['records']), 'fault_injections': probes, 'all_injections_rejected': True, 'sufficiency': suff, 'scientific_status_promotions': 0, 'physical_measurements': 0, 'linear_sensitivity_claims': 0, 'rigorous_physical_or_rounding_enclosure': 'MISSING; original ideal-model enclosures keep their stated assumptions'}
    write(HERE / 'VALIDATION.json', report)
    return report
if __name__ == '__main__':
    validate('--partial' not in __import__('sys').argv)
