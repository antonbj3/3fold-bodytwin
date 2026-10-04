"""Default CLI fails every OPEN/UNKNOWN edge without a finite sourced gap.

Integrity and completeness are deliberately separate: --integrity-only does
not certify the user-requested numeric coverage. All reported failures persist.
"""
import argparse
import copy
import math
from fractions import Fraction
from common import *
from build_net import CORRECTED

def number(x):
    return isinstance(x, (int, float)) and (not isinstance(x, bool)) and math.isfinite(x)

def gap_errors(g, label):
    errs = []
    if not number(g.get('value')) or g.get('value', -1) < 0:
        return [label + ': MISSING_FINITE_NUMERIC_GAP']
    if not isinstance(g.get('unit'), str) or not g['unit'].strip():
        errs.append(label + ': MISSING_UNIT')
    ops = g.get('operands', [])
    if not ops:
        return errs + [label + ': UNSOURCED_GAP']
    vals = []
    for r in ops:
        try:
            got = ev_value(r)
            if got != r['value']:
                errs.append(label + ': SOURCE_VALUE_MISMATCH')
            if r['unit'] != g['unit']:
                errs.append(label + ': OPERAND_UNIT_MISMATCH')
            vals.append(Fraction(str(got)))
        except (ValueError, KeyError, FileNotFoundError, TypeError) as exc:
            errs.append(label + ': SOURCE_ERROR ' + str(exc))
    if len(vals) != len(ops):
        return errs
    mode = g.get('operation')
    if mode == 'difference' and len(vals) == 2:
        r = vals[0] - vals[1]
    elif mode == 'reported_residual' and len(vals) == 1:
        r = vals[0]
    elif mode == 'maximum':
        r = max((abs(v) for v in vals))
    elif mode == 'complement100':
        r = Fraction(100) - sum(vals)
    else:
        return errs + [label + ': INVALID_RESIDUAL_OPERATION']
    if float(abs(r)) != g['value']:
        errs.append(label + ': WRONG_NUMERIC_GAP')
    if float(r) != g.get('signed_value'):
        errs.append(label + ': WRONG_SIGNED_RESIDUAL')
    try:
        if Fraction(g['exact_decimal']) != abs(r):
            errs.append(label + ': WRONG_EXACT_DECIMAL')
    except (KeyError, ValueError):
        errs.append(label + ': MISSING_EXACT_DECIMAL')
    scale = g.get('normalization_scale')
    if scale:
        try:
            if ev_value(scale) != scale['value'] or not number(scale['value']) or scale['value'] <= 0 or (scale['unit'] != g['unit']):
                errs.append(label + ': INVALID_NORMALIZATION_SCALE')
        except Exception:
            errs.append(label + ': INVALID_NORMALIZATION_SOURCE')
    return errs

def validate(net):
    integrity = []
    missing = []
    seen = set()
    numeric = 0
    inherited = {e['id']: e for e in source('results/PROOF_LANE_CONSTRAINT_NET_DENTAL/CONSTRAINT_NET_DENTAL.json')['edges']}
    lock = read(HERE / 'INPUT_LOCK.json')['files']
    for e in net['edges']:
        label = e['id']
        if label in seen:
            integrity.append(label + ': DUPLICATE_ID')
        seen.add(label)
        if label in inherited and e['status'] != inherited[label]['status']:
            integrity.append(label + ': INHERITED_STATUS_PROMOTION')
        if any((v not in net['variables'] for v in e['between'])):
            integrity.append(label + ': UNKNOWN_ENDPOINT')
        if e['timescale'] not in ['SIMULTANEOUS', 'HANDOVER']:
            integrity.append(label + ': TIMESCALE')
        if e['decision_count'] != len(set(e['decision_ids'])):
            integrity.append(label + ': DECISION_COUNT')
        if label in CORRECTED:
            if e['resolution_level'] != 'POPULATION' or e['evidence_resolution_level'] != 'POPULATION' or any((r['resolution_level'] != 'POPULATION' for r in e['evidence'])):
                integrity.append(label + ': AGGREGATION_PROMOTION')
            if e.get('anatomical_support_level') not in ['PER_TOOTH', 'PER_SURFACE_REGION']:
                integrity.append(label + ': ANATOMICAL_SUPPORT_MISSING')
        for r in e.get('evidence', []):
            try:
                if ev_value(r) != r['value']:
                    integrity.append(label + ': EVIDENCE_VALUE_MISMATCH')
            except Exception as exc:
                integrity.append(label + ': EVIDENCE_ERROR ' + str(exc))
        if e['status'] in ['OPEN', 'UNKNOWN']:
            g = e.get('gap', {})
            errors = gap_errors(g, label)
            if any(('MISSING_FINITE' in x for x in errors)):
                missing.extend(errors)
            else:
                numeric += 1
                integrity.extend(errors)
        for (i, g) in enumerate(e.get('residual_field', []) + e.get('gap_components', [])):
            integrity.extend(gap_errors(g, label + ':component:' + str(i)))
        if e.get('r2_origin') not in ['INHERITED_R1', None]:
            if e.get('source_review_state') != 'REVIEWED' or not any((r['hash_match'] and r['decision'] in ['ACCEPT', 'ACCEPT_WITH_CORRECTION'] for r in e['review'])):
                integrity.append(label + ': UNREVIEWED_MAIN_EDGE')
            for review in e['review']:
                try:
                    actual = source(review['record'])
                    if lock[review['record']]['sha256'] != review['record_sha256']:
                        integrity.append(label + ': REVIEW_HASH')
                    if actual['result_sha256'] != lock[actual['result_file']]['sha256']:
                        integrity.append(label + ': REVIEW_RESULT_HASH')
                    if actual.get('graph_decision', actual.get('decision')) != review['decision']:
                        integrity.append(label + ': REVIEW_DECISION')
                except Exception as exc:
                    integrity.append(label + ': REVIEW_ERROR ' + str(exc))
    if len([e for e in net['edges'] if e['id'] in CORRECTED]) != 10:
        integrity.append('MISSING_CORRECTED_EDGE')
    for (path, h) in read(HERE / 'R1_READ_ONLY_MANIFEST.json').items():
        if sha(path) != h:
            integrity.append('R1_CHANGED: ' + path)
    return {'integrity_pass': not integrity, 'numeric_gap_complete': not missing, 'pass': not integrity and (not missing), 'numeric_open_unknown_edges': numeric, 'missing_numeric_open_unknown_edges': len(missing), 'integrity_errors': integrity, 'completeness_errors': missing, 'clinical_or_physical_admission': False}

def mutations(net):
    reports = []

    def check(name, change, code):
        n = copy.deepcopy(net)
        change(n)
        out = validate(n)
        es = out['integrity_errors'] + out['completeness_errors']
        reports.append({'fault': name, 'rejected': any((code in s for s in es)), 'expected_error': code, 'observed_errors': [s for s in es if code in s]})

    def target(n):
        return next((e for e in n['edges'] if e['id'] == 'D-E-TANDLAST-PARTIAL-REFERENCE'))
    check('actual gap +1 percentage point', lambda n: target(n)['gap'].__setitem__('value', 8.0), 'WRONG_NUMERIC_GAP')
    check('remove numeric gap', lambda n: target(n).pop('gap'), 'MISSING_FINITE_NUMERIC_GAP')
    check('null residual', lambda n: target(n)['gap'].__setitem__('value', None), 'MISSING_FINITE_NUMERIC_GAP')
    check('boolean masquerading as gap', lambda n: target(n)['gap'].__setitem__('value', True), 'MISSING_FINITE_NUMERIC_GAP')
    check('infinite residual', lambda n: target(n)['gap'].__setitem__('value', float('inf')), 'MISSING_FINITE_NUMERIC_GAP')
    check('wrong unit', lambda n: target(n)['gap'].__setitem__('unit', 'N'), 'OPERAND_UNIT_MISMATCH')
    check('actual source operand +1', lambda n: target(n)['gap']['operands'][0].__setitem__('value', 1 + target(n)['gap']['operands'][0]['value']), 'SOURCE_VALUE_MISMATCH')
    check('corrupt source SHA', lambda n: target(n)['gap']['operands'][0].__setitem__('source_sha256', '0' * 64), 'SOURCE_HASH')
    check('reverse residual sign', lambda n: target(n)['gap'].__setitem__('signed_value', -7.0), 'WRONG_SIGNED_RESIDUAL')
    check('inflate decision multiplicity', lambda n: target(n).__setitem__('decision_count', 999), 'DECISION_COUNT')
    check('false individual group observation', lambda n: next((e for e in n['edges'] if e['id'] == 'D-E-CONTRAST-03')).__setitem__('resolution_level', 'PER_POINT'), 'AGGREGATION_PROMOTION')
    check('admit unreviewed candidate', lambda n: n['edges'].append(copy.deepcopy(n['pending_unreviewed_edges'][0])), 'UNREVIEWED_MAIN_EDGE')
    return {'controls': reports, 'all_rejected': all((r['rejected'] for r in reports)), 'count': len(reports)}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--net', default='CONSTRAINT_NET_DENTAL_R2.json')
    ap.add_argument('--integrity-only', action='store_true')
    ap.add_argument('--mutations', action='store_true')
    a = ap.parse_args()
    n = read(HERE / a.net)
    result = validate(n)
    if a.mutations:
        result['fault_controls'] = mutations(n)
    dump('VALIDATION.json', result)
    print(json.dumps({k: v for (k, v) in result.items() if not isinstance(v, (list, dict))}))
    ok = result['integrity_pass'] if a.integrity_only else result['pass']
    if a.mutations:
        ok = ok and result['fault_controls']['all_rejected']
    raise SystemExit(0 if ok else 1)
if __name__ == '__main__':
    main()
