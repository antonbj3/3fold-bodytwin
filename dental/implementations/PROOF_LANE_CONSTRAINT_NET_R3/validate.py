"""R3 validator: require sourced numeric gap or explicit physical-observation contract.
Nonphysical exceptions remain strict failures; normalization is source-bound (batch20 correction)."""
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

def validate_base(net):
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
        if e['id'].startswith('D-E-X70-FORCE-') or any((k in e for k in ['comparison_error', 'comparison_tolerance'])):
            for required in ['comparison_error', 'comparison_tolerance']:
                if required not in e:
                    integrity.append(label + ': MISSING_NORMALIZATION_FIELD:' + required)
        for field in ['comparison_error', 'comparison_tolerance']:
            if field not in e:
                continue
            row = e[field]
            try:
                actual = ev_value(row)
                if not number(row.get('value')) or actual != row['value']:
                    integrity.append(label + ': NORMALIZATION_SOURCE_VALUE_MISMATCH:' + field)
                if row.get('unit') != '1':
                    integrity.append(label + ': NORMALIZATION_UNIT:' + field)
                if field == 'comparison_tolerance' and (not number(actual) or actual <= 0):
                    integrity.append(label + ': INVALID_COMPARISON_TOLERANCE')
                if field == 'comparison_error' and (not number(actual) or actual < 0):
                    integrity.append(label + ': INVALID_COMPARISON_ERROR')
            except Exception as exc:
                integrity.append(label + ': NORMALIZATION_SOURCE_ERROR:' + field + ':' + str(exc))
        if e['status'] in ['OPEN', 'UNKNOWN']:
            g = e.get('gap', {})
            errors = [] if g.get('kind') == 'BLOCKED_ON_PHYSICAL_MEASUREMENT' else gap_errors(g, label)
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
    return {'integrity_pass': not integrity, 'numeric_gap_complete': not missing, 'pass': not integrity and (not missing), 'numeric_open_unknown_edges': numeric, 'missing_numeric_open_unknown_edges': len(missing), 'integrity_errors': integrity, 'completeness_errors': missing, 'clinical_or_physical_admission': False}

def validate(net):
    result = validate_base(net)
    errors = result['integrity_errors']
    contracts = read(HERE / 'MEASUREMENT_CONTRACTS.json')
    original = source('results/PROOF_LANE_CONSTRAINT_NET_R2/CONSTRAINT_NET_DENTAL_R2.json')
    orig = {e['id']: e for e in original['edges']}
    got = {e['id']: e for e in net['edges']}
    if set(got) != set(orig):
        errors.append('EDGE_CENSUS')
    if net.get('measurement_groups') != contracts['groups']:
        errors.append('MEASUREMENT_GROUP_CONTRACT')
    for (id, e) in got.items():
        if id not in orig:
            continue
        o = orig[id]
        if e['status'] != o['status']:
            errors.append(id + ': INHERITED_STATUS_PROMOTION')
        for f in ['between', 'inputs', 'outputs', 'timescale', 'decision_ids', 'decision_count']:
            if e.get(f) != o.get(f):
                errors.append(id + ': INHERITED_CONSUMER_OR_TIME:' + f)
        if o.get('gap', {}).get('value') is not None and e['gap'] != o['gap']:
            errors.append(id + ': FROZEN_NUMERIC_GAP_CONTRACT')
        if o.get('gap_components') and e.get('gap_components') != o.get('gap_components'):
            errors.append(id + ': FROZEN_COMPONENT_CONTRACT')
        if id in contracts['edges']:
            c = contracts['edges'][id]
            if e.get('measurement_obligation') != c:
                errors.append(id + ': MEASUREMENT_OBLIGATION_CONTRACT')
            if e.get('gap', {}).get('kind') != c['status'] or e.get('blocker_status') != c['status']:
                errors.append(id + ': BLOCKER_STATUS')
            if c['status'] == 'BLOCKED_ON_PHYSICAL_MEASUREMENT':
                for k in ['required_measurement', 'specimen', 'would_decide', 'measurement_ids']:
                    if not e['gap'].get(k):
                        errors.append(id + ': MISSING_PHYSICAL_' + k.upper())
                    ck = 'measurement' if k == 'required_measurement' else k
                    if e['gap'].get(k) != c.get(ck):
                        errors.append(id + ': PHYSICAL_GAP_METADATA_MISMATCH:' + k)
                if e['gap'].get('value') is not None:
                    errors.append(id + ': FABRICATED_PHYSICAL_GAP')
                for g in c['measurement_ids']:
                    if g not in contracts['groups']:
                        errors.append(id + ': UNKNOWN_MEASUREMENT_GROUP')
                for r in c['source_of_missingness']:
                    try:
                        if ev_value(r) != r['value']:
                            errors.append(id + ': MISSINGNESS_SOURCE_VALUE')
                    except Exception:
                        errors.append(id + ': MISSINGNESS_SOURCE_ERROR')
        for field in ['comparison_error', 'comparison_tolerance']:
            if field in o and e.get(field) != o[field]:
                errors.append(id + ': NORMALIZATION_FROZEN_CONTRACT:' + field)
        for context in e.get('r3_verified_context', []):
            try:
                row = next((r for r in read(HERE / 'PROPOSAL_AUDIT.json')['rows'] if r['id'] == context['proposal_id']))
                p = row['proposal']
                s = row['semantic_decision']
                ev = context['evidence']
                if ev_value(ev) != ev['value']:
                    errors.append(id + ': CONTEXT_SOURCE_VALUE')
                if ev['key'] != p['key_path'] or ev['unit'] != p['unit'] or ev['resolution_level'] != s['resolution_level'] or (context['supports_full_edge_physical_closure'] is not False):
                    errors.append(id + ': CONTEXT_ESTIMAND')
            except Exception:
                errors.append(id + ': CONTEXT_SOURCE_ERROR')
    if 'D-E-CANAL-RELEASE' in got:
        e = got['D-E-CANAL-RELEASE']
        s = source('results/LANE_X58_CANAL_WALL_SPREAD/raw/PER_SITE_LOCATED.json')
        f = e.get('located_release_field', [])
        if len(f) != len(s):
            errors.append('X58_LOCATED_FIELD_CENSUS')
        for (i, r) in enumerate(f):
            try:
                v = ev_value(r['evidence'])
                if not number(r['evidence']['value']) or v != r['evidence']['value'] or v != s[i]['old_minus_new_mid_mm']:
                    errors.append('X58_LOCATED_VALUE')
                if r['instance'] != {'case': s[i]['case'], 'fdi': s[i]['fdi'], 'pose': s[i]['pose']} or r['spatial_witness'] != {'old': s[i]['old'].get('witness'), 'new': s[i]['new'].get('witness')}:
                    errors.append('X58_LOCATED_JOIN')
                if r['evidence']['unit'] != 'mm' or r['observable_resolution'] != 'PER_TOOTH' or r['binding_resolution'] != 'PER_POINT':
                    errors.append('X58_LOCATED_UNIT_OR_RESOLUTION')
            except Exception:
                errors.append('X58_LOCATED_SOURCE')
        for g in [e['gap']] + e['gap_components']:
            if g['resolution_level'] != 'POPULATION' or g['unit'] != 'mm' or g.get('measured_physical_gap') is not False:
                errors.append('X58_SUMMARY_SCOPE')
        if e['gap']['operands'][0]['key'] != '/R3/max_located_old_minus_new_mm':
            errors.append('X58_ESTIMAND')
        if not e.get('physical_truth_obligation', {}).get('specimen'):
            errors.append('X58_TRUTH_DEBT')
    physical = [e['id'] for e in net['edges'] if e.get('gap', {}).get('kind') == 'BLOCKED_ON_PHYSICAL_MEASUREMENT']
    numeric = [e['id'] for e in net['edges'] if e['status'] in ['OPEN', 'UNKNOWN'] and number(e.get('gap', {}).get('value'))]
    unhandled = [e['id'] for e in net['edges'] if e['status'] in ['OPEN', 'UNKNOWN'] and e['id'] not in physical + numeric]
    result.update({'integrity_pass': not errors, 'pass': not errors and (not unhandled), 'numeric_gap_complete': not physical and (not unhandled), 'availability_contract_complete': not unhandled, 'numeric_open_unknown_edges': len(numeric), 'physically_blocked_edges': len(physical), 'missing_numeric_open_unknown_edges': len(physical) + len(unhandled), 'unhandled_open_unknown_edges': unhandled, 'physical_edges': physical, 'numeric_edges': numeric, 'clinical_or_physical_admission': False})
    return result

def mutations(net):
    reports = []

    def e(n, id):
        return next((x for x in n['edges'] if x['id'] == id))

    def check(name, fn, code):
        n = copy.deepcopy(net)
        fn(n)
        r = validate(n)
        err = r['integrity_errors'] + r['completeness_errors']
        reports.append({'fault': name, 'rejected': any((code in x for x in err)), 'expected': code, 'observed': [x for x in err if code in x]})
    t = 'D-E-TANDLAST-PARTIAL-REFERENCE'
    p = 'D-E-K04'
    norm = 'D-E-X70-FORCE-1'
    check('gap+1', lambda n: e(n, t)['gap'].__setitem__('value', 8), 'WRONG_NUMERIC_GAP')
    check('gap=null', lambda n: e(n, t)['gap'].__setitem__('value', None), 'MISSING_FINITE_NUMERIC_GAP')
    check('gap=true', lambda n: e(n, t)['gap'].__setitem__('value', True), 'MISSING_FINITE_NUMERIC_GAP')
    check('gap=inf', lambda n: e(n, t)['gap'].__setitem__('value', float('inf')), 'MISSING_FINITE_NUMERIC_GAP')
    check('gap source wrong', lambda n: e(n, t)['gap']['operands'][0].__setitem__('value', 99), 'SOURCE_VALUE_MISMATCH')
    check('gap hash wrong', lambda n: e(n, t)['gap']['operands'][0].__setitem__('source_sha256', '0' * 64), 'SOURCE_HASH')
    check('gap unit wrong', lambda n: e(n, t)['gap'].__setitem__('unit', 'N'), 'OPERAND_UNIT_MISMATCH')
    check('gap sign reversed', lambda n: e(n, t)['gap'].__setitem__('signed_value', -7), 'WRONG_SIGNED_RESIDUAL')
    check('decision count inflate', lambda n: e(n, t).__setitem__('decision_count', 99), 'DECISION_COUNT')
    check('population→point', lambda n: e(n, 'D-E-CONTRAST-03').__setitem__('resolution_level', 'PER_POINT'), 'AGGREGATION_PROMOTION')
    for (field, part, value) in [('comparison_error', 'value', 0), ('comparison_tolerance', 'value', 20), ('comparison_error', 'unit', 'N'), ('comparison_tolerance', 'source_sha256', '0' * 64)]:
        check('normalization ' + field + '/' + part, lambda n, f=field, k=part, v=value: e(n, norm)[f].__setitem__(k, v), 'NORMALIZATION_')
    for field in ['comparison_error', 'comparison_tolerance']:
        check('remove ' + field, lambda n, f=field: e(n, norm).pop(f), 'MISSING_NORMALIZATION_FIELD')
    for field in ['measurement', 'specimen', 'would_decide', 'measurement_ids']:
        check('remove physical ' + field, lambda n, f=field: e(n, p)['measurement_obligation'].pop(f), 'MEASUREMENT_OBLIGATION_CONTRACT')
    check('wrong nonempty physical specimen', lambda n: e(n, p)['gap'].__setitem__('specimen', 'different specimen'), 'PHYSICAL_GAP_METADATA_MISMATCH')
    check('X58 boolean masquerading as zero delta', lambda n: e(n, 'D-E-CANAL-RELEASE')['located_release_field'][0]['evidence'].__setitem__('value', False), 'X58_LOCATED_VALUE')
    check('empty physical specimen', lambda n: e(n, p)['gap'].__setitem__('specimen', ''), 'MISSING_PHYSICAL_SPECIMEN')
    check('nonphysical→physical', lambda n: e(n, 'D-E-K52')['gap'].__setitem__('kind', 'BLOCKED_ON_PHYSICAL_MEASUREMENT'), 'BLOCKER_STATUS')
    check('invent physical zero', lambda n: e(n, p)['gap'].__setitem__('value', 0), 'FABRICATED_PHYSICAL_GAP')
    check('remove unsatisfied edge', lambda n: n['edges'].remove(e(n, 'D-E-K52')), 'EDGE_CENSUS')
    check('remove X58 located field', lambda n: e(n, 'D-E-CANAL-RELEASE')['located_release_field'].pop(), 'X58_LOCATED_FIELD_CENSUS')
    check('X58 wrong pose', lambda n: e(n, 'D-E-CANAL-RELEASE')['located_release_field'][0]['instance'].__setitem__('fdi', 99), 'X58_LOCATED_JOIN')
    check('X58 aggregate promoted to point', lambda n: e(n, 'D-E-CANAL-RELEASE')['gap'].__setitem__('resolution_level', 'PER_POINT'), 'X58_SUMMARY_SCOPE')
    check('wrong physical protocol', lambda n: n['measurement_groups']['M01'].__setitem__('measurement', 'total mean only'), 'MEASUREMENT_GROUP_CONTRACT')
    return {'controls': reports, 'count': len(reports), 'all_rejected': all((x['rejected'] for x in reports))}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--net', default='CONSTRAINT_NET_DENTAL_R3.json')
    ap.add_argument('--integrity-only', action='store_true')
    ap.add_argument('--mutations', action='store_true')
    a = ap.parse_args()
    verify_frozen()
    n = read(HERE / a.net)
    r = validate(n)
    if a.mutations:
        r['fault_controls'] = mutations(n)
    dump('VALIDATION.json', r)
    print(json.dumps({k: v for (k, v) in r.items() if not isinstance(v, (list, dict))}))
    ok = r['integrity_pass'] if a.integrity_only else r['pass']
    if a.mutations:
        ok = ok and r['fault_controls']['all_rejected']
    raise SystemExit(0 if ok else 1)
if __name__ == '__main__':
    main()
