#!/usr/bin/env python3
"""At what scale does each quantity in the net live, and which edges cross scales without saying so?

The measured state, 2026-10-05. The net has 166 edges and 148 variables. Only 36 variables declare a
resolution level at all; 112 say UNKNOWN, and four of the declared ones use lowercase spellings of a
level the rest write in capitals. The consequence is structural, not cosmetic: only 21 edges have a
known level at BOTH ends, so 145 cannot be checked for a scale jump even in principle.

Of the 21 that can be checked, 8 connect quantities at different levels. Three of those are declared
TIGHT:

    T-E1   charge inventory (MOLECULE) -- ionic swelling pressure (TISSUE)
    T-E14  C3b recycling fraction (MOLECULE) -- self surface growth rate (CELL)
    T-E15  C3b recycling fraction (MOLECULE) -- activator/host discrimination (CELL)

A TIGHT edge across a scale jump asserts that a molecular quantity determines a tissue or cell
quantity. That assertion needs a stated averaging step -- how many molecules, over what volume, with
what assumption about their arrangement -- and none of the three states one. The edge may well hold;
what it cannot do is hold TIGHTLY without the step that carries it being written down.

So this audit does three things and refuses a fourth. It normalises the spellings. It reports level
coverage, because a twin that cannot say what scale three quarters of its own quantities live at
cannot tell a tissue average from a molecular count. And it flags every edge that crosses a level
without declaring the homogenisation that bridges it. It does NOT guess a level from a unit: MPa is a
tissue stress in one edge and a molecular interaction energy density in another, and a guessed level
would be worse than a missing one, because the cross-scale check would then run on invented data.
"""
from __future__ import annotations

import collections
import json
from pathlib import Path

W = Path(__file__).resolve().parents[2]
OUT = W / 'results/ASSEMBLY_RESOLUTION'
# Coarsest last. The order is what makes "crosses a scale" and "how far" answerable.
ORDER = ['MOLECULE', 'FIBRIL', 'CELL', 'TISSUE', 'ORGAN', 'PHENOMENOLOGICAL']
# A field that would state the bridge. Any one of them is enough to count as declared.
BRIDGE_KEYS = ('homogenisation', 'homogenization', 'averaging', 'scale_bridge',
               'resolution_note', 'coarse_graining')


def walk(obj, key):
    if isinstance(obj, dict):
        if isinstance(obj.get(key), list):
            for e in obj[key]:
                if isinstance(e, dict):
                    yield e
        for v in obj.values():
            yield from walk(v, key)
    elif isinstance(obj, list):
        for v in obj:
            yield from walk(v, key)


def main() -> int:
    nets = ['CONSTRAINT_NETS.json', 'data/CONSTRAINT_NET_TISSUE.json']
    normalised = 0
    for f in nets:
        p = W / f
        d = json.loads(p.read_text())
        def fix(o):
            nonlocal normalised
            if isinstance(o, dict):
                for k in ('level', 'resolution_level'):
                    v = o.get(k)
                    if isinstance(v, str) and v and v != v.upper() and v.upper() in ORDER:
                        o[k] = v.upper()
                        normalised += 1
                for v in o.values():
                    fix(v)
            elif isinstance(o, list):
                for v in o:
                    fix(v)
        fix(d)
        p.write_text(json.dumps(d, indent=1, ensure_ascii=False) + '\n')

    d = json.loads((W / 'CONSTRAINT_NETS.json').read_text())
    edges = list(walk(d, 'edges'))
    variables = list(walk(d, 'variables'))
    lvl = {v['id']: str(v.get('level') or v.get('resolution_level') or 'UNKNOWN').upper()
           for v in variables if v.get('id')}

    both_known, one_missing, crossing = 0, 0, []
    for e in edges:
        ends = e.get('between') or []
        levels = [lvl.get(x, 'MISSING') for x in ends]
        if ends and all(l in ORDER for l in levels):
            both_known += 1
            if len(set(levels)) > 1:
                span = max(ORDER.index(l) for l in levels) - min(ORDER.index(l) for l in levels)
                crossing.append({
                    'edge': e.get('id'), 'status': e.get('status'),
                    'between': ends, 'levels': levels, 'levels_apart': span,
                    'bridge_declared': any(e.get(k) for k in BRIDGE_KEYS),
                    'tight_without_a_bridge': (e.get('status') == 'TIGHT'
                                               and not any(e.get(k) for k in BRIDGE_KEYS)),
                })
        else:
            one_missing += 1

    report = {
        'variables': len(lvl),
        'variables_without_a_level': sum(1 for v in lvl.values() if v == 'UNKNOWN'),
        'level_counts': dict(collections.Counter(lvl.values())),
        'lowercase_spellings_normalised': normalised,
        'edges': len(edges),
        'edges_with_a_known_level_at_both_ends': both_known,
        'edges_unknown_at_one_end_or_more': one_missing,
        'edges_crossing_a_level': crossing,
        'tight_edges_crossing_a_level_without_a_declared_bridge':
            [c['edge'] for c in crossing if c['tight_without_a_bridge']],
        'why_no_level_is_inferred_from_units': ('MPa is a tissue stress in one edge and a molecular '
                                               'interaction energy density in another; a guessed '
                                               'level would make the cross-scale check run on '
                                               'invented data'),
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'AUDIT_V1.json').write_text(json.dumps(report, indent=1, ensure_ascii=False))

    print(f"  variabler {report['variables']}, without level {report['variables_without_a_level']} ({100 * report['variables_without_a_level'] / report['variables']:.0f} %), normaliserade stavningar {normalised}")
    print(f'  kanter {len(edges)}: level known at both ends {both_known}, unknown in at least one {one_missing}')
    print(f"  crosses a level: {len(crossing)}, varav TIGHT without the bridge declared: {len(report['tight_edges_crossing_a_level_without_a_declared_bridge'])}")
    for c in crossing:
        mark = 'TIGHT WITHOUT BRYGGA' if c['tight_without_a_bridge'] else ''
        print(f"   {c['levels_apart']} steg  {str(c['levels']):32} {c['status']:8} "
              f"{c['edge'][:44]:44} {mark}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
