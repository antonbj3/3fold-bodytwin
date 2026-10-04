"""Freeze a fixed plan, then evaluate all planned IDs with conservative censoring."""
import argparse, datetime, json, math
from pathlib import Path
from common import sha, clean
from acceptance import plan, upper

def dump(p, o):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(clean(o), indent=2, allow_nan=False) + '\n')

def freeze(context, output):
    required = ['lot_id', 'configuration', 'protocol_locator', 'protocol_sha256', 'quantity', 'unit', 'proof_threshold', 'planned_ids', 'p_bad', 'p_good', 'power', 'alpha']
    if any((k not in context for k in required)):
        raise ValueError('Missing fixed-plan context fields')
    if context['quantity'] not in ['cycles', 'strength'] or context['unit'] != {'cycles': 'cycles', 'strength': 'MPa'}[context['quantity']]:
        raise ValueError('Quantity/unit mismatch')
    if not math.isfinite(context['proof_threshold']) or context['proof_threshold'] <= 0:
        raise ValueError('Invalid proof threshold')
    for k in ['protocol_locator', 'protocol_sha256', 'lot_id', 'configuration']:
        if not context[k]:
            raise ValueError('Missing protocol/lot identity')
    if len(context['protocol_sha256']) != 64 or any((c not in '0123456789abcdef' for c in context['protocol_sha256'])):
        raise ValueError('Protocol SHA256 required')
    p = plan(context['p_bad'], context['p_good'], context['power'], context['alpha'])
    ids = context['planned_ids']
    if len(ids) != p['n'] or len(set(ids)) != len(ids) or any((not isinstance(i, str) or not i for i in ids)):
        raise ValueError('Supply exactly n unique string IDs before testing')
    if Path(output).exists():
        raise ValueError('Do not overwrite a frozen plan')
    o = {'frozen_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'claim_type': 'capability', 'plan': p, 'context': context, 'physical_eligibility': 'LAB_MUST_VERIFY_REGIME', 'optional_acceptance': False}
    dump(output, o)
    Path(str(output) + '.sha256').write_text(sha(output) + '\n')
    return o

def evaluate(planfile, records):
    if sha(planfile) != Path(str(planfile) + '.sha256').read_text().strip():
        raise ValueError('Frozen plan hash changed')
    p = json.loads(Path(planfile).read_text())
    ctx = p['context']
    pred = p['plan']
    seen = {}
    bad = good = unknown = 0
    for row in records:
        if not row.get('locator'):
            raise ValueError('Each record needs source locator')
        ident = row['specimen_id']
        if ident not in ctx['planned_ids'] or ident in seen:
            raise ValueError('Unplanned or duplicate specimen')
        for k in ['lot_id', 'configuration', 'protocol_sha256', 'unit']:
            if row.get(k) != ctx[k]:
                raise ValueError('Lot/configuration/protocol/unit mismatch')
        v = row['value']
        e = row['failure']
        if e not in (0, 1) or not isinstance(v, (int, float)) or (not math.isfinite(v)) or (v <= 0):
            raise ValueError('Invalid event/value')
        seen[ident] = row
        if e and v <= ctx['proof_threshold']:
            bad += 1
        elif v >= ctx['proof_threshold']:
            good += 1
        else:
            unknown += 1
    missing = pred['n'] - len(seen)
    worst = bad + unknown + missing
    if missing:
        verdict = 'INCOMPLETE_FIXED_PLAN'
    elif worst <= pred['failures_allowed']:
        verdict = 'ACCEPT_RISK_CRITERION_ONLY'
    else:
        verdict = 'NOT_CERTIFIED'
    return {'plan_sha256': sha(planfile), 'verdict': verdict, 'bad': bad, 'known_survivors': good, 'unknown_endpoints': unknown, 'missing_ids': missing, 'conservative_bad_count': worst, 'upper_at_fixed_n': upper(pred['n'], worst, ctx['alpha']), 'standard_conformity': 'UNKNOWN; this is a probability criterion for the stated lab population', 'resolution': 'POPULATION', 'optional_stopping': 'No acceptance before all planned IDs have records; unresolved endpoints treated as bad'}

def main():
    a = argparse.ArgumentParser(description=__doc__)
    s = a.add_subparsers(dest='cmd', required=True)
    p = s.add_parser('freeze')
    p.add_argument('--context', required=True)
    p.add_argument('--output', required=True)
    p = s.add_parser('evaluate')
    p.add_argument('--plan', required=True)
    p.add_argument('--measurements', required=True)
    args = a.parse_args()
    if args.cmd == 'freeze':
        o = freeze(json.loads(Path(args.context).read_text()), args.output)
    else:
        o = evaluate(args.plan, json.loads(Path(args.measurements).read_text()))
    print(json.dumps(o, indent=2, allow_nan=False))
if __name__ == '__main__':
    main()
