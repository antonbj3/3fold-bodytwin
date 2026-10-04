"""Read the existing heterogeneous receipt schemas without changing their gates."""

def validator_pass(receipt, ident):
    if isinstance(receipt, list):
        matched = [r for r in receipt if isinstance(r, dict) and r.get('id') == ident]
        return len(matched) == 1 and matched[0].get('returncode') == 0 and (matched[0].get('source_trees_hidden') is True) and (matched[0].get('network_disabled') is True)
    if not isinstance(receipt, dict):
        return False
    if 'status' in receipt and receipt['status'] != 'PASS':
        return False
    if receipt.get('pass') is False:
        return False
    if 'demo' in receipt:
        import re
        refs = ['results.json', 'REVIEWED_RESULTS.json', 'INDEPENDENT_VERIFICATION.json', 'REPLAY_LOCK.json']
        return ident == 'PATIENT360' and receipt.get('demo') == ident and (receipt.get('exit_code') == 0) and (receipt.get('decision') == 'DEMO_READY_WITH_CORRECTION') and (receipt.get('templates_unchanged') is True) and (receipt.get('physical_release') == 'ABSTAIN') and all((isinstance(receipt.get(n), dict) and receipt[n].get('path') and re.fullmatch('[0-9a-f]{64}', receipt[n].get('sha256', '')) for n in refs))
    matched = []

    def visit(x):
        if isinstance(x, dict):
            if x.get('demo_id', x.get('id')) == ident:
                matched.append(x)
            for value in x.values():
                visit(value)
        elif isinstance(x, list):
            for value in x:
                visit(value)
    visit(receipt)
    if matched:
        return all(((x.get('status') == 'PASS' or x.get('pass') is True) and ('status' not in x or x['status'] == 'PASS') and (x.get('pass') is not False) for x in matched))
    if ident in receipt.get('demos', []) or receipt.get('demo_id') == ident:
        return receipt.get('status') == 'PASS'
    return False
