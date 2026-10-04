"""Executed case/group/hash leakage checks, conditional on upstream identity."""
import copy
from collections import Counter

def violations(rows, donors):
    errors = []
    seen = {}
    keys = set()
    for r in rows:
        if r['case_key'] in keys:
            errors.append('duplicate case key')
        keys.add(r['case_key'])
        if r['dataset'] == 'Bits2Bites':
            if r['split'] != 'auxiliary':
                errors.append('auxiliary in primary split')
            continue
        if r['legacy_v2'] and r['split'] != 'dev':
            errors.append('legacy case outside development')
        for token in ['group:' + r['group']] + ['mesh:' + h for h in r['geometry_hashes']]:
            if token in seen and seen[token] != r['split']:
                errors.append('primary overlap ' + token)
            seen[token] = r['split']
    expected = {r['group'] for r in rows if r['dataset'] == 'Bite2Text' and r['split'] == 'train'}
    if set(donors) != expected:
        errors.append('training donors differ from primary train')
    return errors

def run(bundle):
    rows = bundle.json('payload/private/COHORT.json')['payload']['cases']
    donors = bundle.json('payload/private/TEMPLATE_FIT.json')['groups']
    positive = not violations(rows, donors)
    train = next((r for r in rows if r['split'] == 'train'))
    test = next((i for (i, r) in enumerate(rows) if r['split'] == 'test'))
    checks = []
    for (name, change) in [('shared patient group', lambda r: r.update(group=train['group'])), ('shared arch geometry', lambda r: r['geometry_hashes'].__setitem__(0, train['geometry_hashes'][0])), ('legacy test case', lambda r: r.update(legacy_v2=True))]:
        bad = copy.deepcopy(rows)
        change(bad[test])
        checks.append(dict(name=name, rejected=bool(violations(bad, donors)), positive=positive))
    checks.append(dict(name='test donor leakage', rejected=bool(violations(rows, donors + [rows[test]['group']])), positive=positive))
    expected = {r['case_key'] for r in rows if r['split'] != 'train'}
    actual = {p.rsplit('/', 1)[-1][:-5] for (p, _) in bundle.tasks()}
    coverage = expected == actual
    if not positive or not coverage or (not all((c['rejected'] for c in checks))):
        raise AssertionError('cohort leakage/coverage gate')
    return dict(all_pass=True, checks=checks, exact_case_coverage=coverage, primary_split_counts=dict(Counter((r['split'] for r in rows if r['dataset'] == 'Bite2Text'))), case_resolution='PER_ARCH', uncertainty_resolution='POPULATION', patient_identity='Conditional on Bite2Text one-case-per-patient documentation; cross-dataset and pretraining identities UNKNOWN')
