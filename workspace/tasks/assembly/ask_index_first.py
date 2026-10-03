"""Before sending a job out to search, ask our own corpus. The graph asks; the index answers.

Why not RAG, measured rather than argued. A word lookup for the quantity an undetermined edge needs
returned 26261 hits on "angle" and 1406 on "mac", and the plausible ones were wrong in a specific way:
"contact" matched a quantity in seconds, "peak density" matched a volume ratio, "graph distance" matched
a timescale. Semantic retrieval over prose makes exactly that failure more likely, because it is built
to return text that reads as relevant. Tonight's repeated lesson is that a plausible match is worse than
no match, since it gets consumed.

What worked instead was dimensional: require the UNIT to agree, then the name. That turned seven word
matches, mostly name collisions, into six candidates that are at least the right kind of quantity. A
unit is a hard constraint that prose similarity is not.

So the two structures have different jobs and neither replaces the other. The net says which quantity a
question needs and carries the gap as a number. The index says whether that quantity, in that unit,
already exists somewhere in 13124 reports. This asks the second before paying for the first: of fifteen
jobs that went out and came back "not published", some of the quantities were already in our own corpus
under a different key.

MEASURED OUTCOME OF THIS VERSION, AND IT IS MOSTLY NOISE. 19 of 58 undetermined variables got a
dimensional hit, and reading them shows why that number means little: `pupil_diameter` in micrometres
matched a fibre-load report about a snow bunting, `refractive_power` in dioptres matched a transport
parameter, and `scatter_tail_fraction` as dimensionless matched 82775 rows. Unit agreement is
necessary and for the common units -- dimensionless, um, D -- it discriminates nothing at all, because
most numbers in a large corpus share a handful of dimensions.

So the conclusion is not that retrieval needs to be smarter. It is that an index reverse-engineered
from arbitrary key names cannot carry meaning the writer never recorded: 37 of 58 variables have no
dimension in their name at all. The fix belongs at WRITE time -- a report declaring its quantities in
a stated schema with unit, tissue and method -- not at read time. Retrofitting gives dimensionally
right and semantically wrong, which is the most expensive kind of match because it gets consumed.

The output is deliberately a RANKED SUGGESTION and not an answer. A unit match is necessary and nowhere
near sufficient -- `contact_area_mm2_median` is the right dimension and may still be the wrong tissue,
and that judgement is not automatable. The job is told what we already hold so it can start from it
rather than from nothing.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

W = Path('.')
NET = W / 'CONSTRAINT_NETS.json'
INDEX = W / 'results/ASSEMBLY_QUANTITY_INDEX/QUANTITY_INDEX_V1.json'
OUT = W / 'results/ASSEMBLY_ASK_INDEX_FIRST'

# Unit written into a variable name, so an edge can be matched dimensionally at all.
UNIT_HINT = [
    ('thickness', 'um'), ('diameter', 'um'), ('temperature', 'K'), ('energy', 'J'),
    ('mass', 'kg'), ('velocity', 'um/h'), ('friction', 'N'), ('concentration', 'nM'),
    ('pressure', 'Pa'), ('power', 'D'), ('cylinder', 'D'), ('astigmatism', 'D'),
    ('length', 'mm'), ('position', 'mm'), ('area', 'mm^2'), ('angle', 'deg'),
    ('dose', 'h'), ('conductance', '1'), ('fraction', '1'), ('ratio', '1'),
]


def unit_for(name: str) -> str | None:
    n = name.lower()
    for word, unit in UNIT_HINT:
        if word in n:
            return unit
    return None


def main() -> int:
    net = json.loads(NET.read_text())['bodytwin']['tissue_constraint_net']
    rows = json.loads(INDEX.read_text())['rows']
    by_unit: dict[str, list[dict]] = {}
    for r in rows:
        by_unit.setdefault(r['unit'], []).append(r)

    report = []
    for e in net['edges']:
        if e['status'] not in ('OPEN', 'UNKNOWN'):
            continue
        for var in e['between']:
            unit = unit_for(var)
            if not unit:
                report.append({'edge': e['id'], 'variable': var, 'unit': None,
                               'verdict': 'NO_UNIT_IN_NAME',
                               'note': 'the variable name carries no dimension, so it cannot be '
                                       'matched dimensionally; that is a defect in the name'})
                continue
            words = [w for w in re.split(r'[_\s]+', var.lower()) if len(w) > 3]
            cands = [r for r in by_unit.get(unit, [])
                     if any(w in r['quantity'].lower() for w in words)]
            cands.sort(key=lambda r: (0 if r['locator'] else 1, -len(r['quantity'])))
            report.append({
                'edge': e['id'], 'variable': var, 'unit': unit,
                'candidates_in_our_corpus': len(cands),
                'with_external_locator': sum(1 for c in cands if c['locator']),
                'verdict': 'ALREADY_HELD_CHECK_BEFORE_SEARCHING' if cands else 'NOT_IN_CORPUS',
                'top': [{'quantity': c['quantity'][:90], 'value': c['value'], 'unit': c['unit'],
                         'job': c['job'], 'locator': c['locator']} for c in cands[:3]],
            })

    held = [r for r in report if r['verdict'] == 'ALREADY_HELD_CHECK_BEFORE_SEARCHING']
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'ASK_FIRST_V1.json').write_text(json.dumps(
        {'why_not_rag': ('a word lookup returned 26261 hits on "angle" and the plausible ones were '
                         'dimensionally wrong; unit agreement is a hard constraint that prose '
                         'similarity is not'),
         'the_two_structures': ('the net says which quantity a question needs and carries its gap as a '
                                'number; the index says whether that quantity in that unit already '
                                'exists in 13124 reports'),
         'output_is_a_suggestion': ('a unit match is necessary and far from sufficient: the right '
                                    'dimension can still be the wrong tissue, which is not automatable'),
         'undetermined_variables_checked': len(report),
         'already_held': len(held),
         'review_state': 'PENDING_INDEPENDENT_REVIEW',
         'rows': report}, indent=1, ensure_ascii=False))

    print(f'  indefinite variables checked: {len(report)}')
    print(f'  already in our own corpus (the right dimension): {len(held)}')
    print(f"  without dimension in the name: {sum((1 for r in report if r['verdict'] == 'NO_UNIT_IN_NAME'))}")
    for r in held[:8]:
        t = r['top'][0]
        print(f"    {r['variable'][:30]:30s} [{r['unit']:>6s}] {r['candidates_in_our_corpus']:>5d} st  "
              f"t.ex. {t['quantity'][:48]} = {t['value']}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
