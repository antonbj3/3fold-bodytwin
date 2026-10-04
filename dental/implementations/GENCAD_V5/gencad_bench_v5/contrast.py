from common import *
import collections, itertools

def is_reversal(a, b):
    if a['key'] != b['key'] or a.get('status') != 'SCORED' or b.get('status') != 'SCORED':
        return False
    qa = a['function']['0']
    qb = b['function']['0']
    if qa['status'] != 'GEOMETRIC_ONLY' or qb['status'] != 'GEOMETRIC_ONLY':
        return False
    return a['reconstruction_p95_mm'] + 0.01 < b['reconstruction_p95_mm'] and qa['penetration_area_mm2'] > qb['penetration_area_mm2'] + 0.05

def run():
    rows = read(ROOT / 'raw/FUNCTIONAL_ROWS.json')
    groups = collections.defaultdict(list)
    for r in rows:
        if r['status'] == 'SCORED':
            groups[r['key']].append(r)
    pairs = []
    n = 0
    for (key, rs) in groups.items():
        for (a, b) in itertools.permutations(rs, 2):
            n += 1
            if is_reversal(a, b):
                pairs.append(dict(key=key, a=a['participant'], b=b['participant'], track_a=a['track'], track_b=b['track'], a_p95_mm=a['reconstruction_p95_mm'], b_p95_mm=b['reconstruction_p95_mm'], a_penetration_mm2=a['function']['0']['penetration_area_mm2'], b_penetration_mm2=b['function']['0']['penetration_area_mm2'], resolution='PER_TOOTH'))
    result = dict(claim_type='capability', candidate_ordered_pairs=n, reversals=len(pairs), cases_with_reversal=len({r['key'] for r in pairs}), rows=pairs, claim='Shape resemblance is insufficient to order geometric interference' if pairs else 'No preregistered ordering reversal found', no_best_crown_claim=True)
    dump(ROOT / 'raw/R5_CONTRAST.json', result)
    return result
if __name__ == '__main__':
    run()
