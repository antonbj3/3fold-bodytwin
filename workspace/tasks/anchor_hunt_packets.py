#!/usr/bin/env python3
"""Turn every orderable acquisition item into a swarm search job for a public dataset.

Why this exists. Five independent lanes converged on the same binding constraint during the night of
2–3 October: what stops us is measurements we do not have, not resolution and not method. The
resolution atlas found 3 of 108 quantities needing finer resolution while 51 were already adequate, and
96.08 % of the undecidable ones had no second quantity to blame. At the same time the acquisition demand
sat scattered as prose in twelve lanes. The harvester collected the well-formed items; this script sends
each one out as a search job with web access, so the swarm looks for a public dataset or a published
measurement in the stated unit.

The operator's rule for this round, in his words: raise resolution, stop at the stress point, bring data
anchors in bulk, take datasets where they exist. A stress point is a quantity whose value decides a
verdict — the eye chain's surface allocation spans 0.87 D between its extremes, the complement cell's
recycling fraction swings discrimination by a factor 98. Those get searched first.

Safety. These briefs leave the machine, so a job is skipped if its text mentions a private target
tissue or a device. Nothing else about the project is included: each brief carries one quantity, one
unit and what the measurement would decide.
"""
from __future__ import annotations

import json
import os
import re
import sys

W = '.'
HARVEST = os.path.join(W, 'notes/ACQUISITION_HARVEST.json')
TARGETS = os.path.join(W, 'notes/ACQUISITION_TARGETS.json')
QUEUE = os.path.join(W, 'tasks/lanes/bt_queue.txt')
RESULTS = os.path.join(W, 'results')

# Never leaves the machine. See tasks/build_night/LIT_REFS_FILTER.md.
BLOCK = re.compile(r'excluded_category|excluded_category|excluded_category|sinusoid|excluded_category|excluded_category|excluded_category|device dosimetry',
                   re.I)

# Quantities already identified as stress points: their value decides a verdict, with the span measured.
STRESS = {
    'surface': 'the eye chain spans 0.87 D between the extremes of this quantity, which is 3.5x the '
               '0.25 D clinical autorefractor threshold, and the sign of the refractive change flips',
    'recycl': 'complement discrimination between host and activator surface swings from 486x to 4.97x, '
              'a factor 98, depending on this quantity',
    'fluence': 'our chain sits 26.5 % above the published 1.15 J/cm^2 ablation threshold and the '
               'mismatch is in the source convention for this quantity',
    'allocation': 'two states with identical mean thickness, index field and hydration differ by '
                  '0.87 D downstream, and this quantity is the candidate missing variable',
}


def slug(s: str) -> str:
    s = re.sub(r'[^A-Za-z0-9]+', '-', s).strip('-').upper()
    return s[:44] or 'ITEM'


def brief(item: dict, stress_note: str | None) -> str:
    q = item.get('quantity', '')
    unit = item.get('unit', '')
    decides = item.get('decides') or item.get('what_it_decides') or ''
    terms = item.get('search_terms') or []
    cls = item.get('quantity_class', '')
    level = item.get('resolution_level', '')
    head = '# Find a public measurement of one quantity\n\n'
    body = [
        f'**Quantity.** {q}',
        f'**Unit required.** `{unit}` — a value in any other unit must be converted and the '
        f'conversion shown, or reported as not comparable.',
    ]
    if cls:
        body.append(f'**Quantity class.** {cls}')
    if level:
        body.append(f'**Resolution level.** {level}')
    if decides:
        body.append(f'**What the measurement decides.** {decides}')
    if stress_note:
        body.append(f'**This is a stress point.** {stress_note} Spend the budget here rather than on '
                    f'breadth.')
    body.append(
        '\n## Deliver a machine-readable table FIRST, then the prose\n'
        'Measured on the first 44 of these searches: the reports carry hundreds of values but only 55 '
        'were extractable, because the citation sits in one cell and the value in another, or in prose. '
        'So begin `RESULTS.md` with one table in exactly these columns, one row per value, and put the '
        'discussion after it:\n\n'
        '| quantity | value | unit | locator | validity | n |\n'
        '|---|---|---|---|---|---|\n'
        '| what was measured | the number alone | the unit alone | DOI or PMID | species/tissue/temp/method | sample size |\n\n'
        'One value per row, no ranges in the value cell (give low and high as two rows), and the locator '
        'in the same row as its value. A row whose unit differs from the one requested above must show '
        'the conversion in the quantity cell or be marked `NOT_COMPARABLE`.\n'
        '\n## What to deliver\n'
        '1. **A value with its unit**, or `VERIFIED_ABSENCE_OF_EVIDENCE` with the searches you ran. '
        'Absence, stated with what you searched, is a full result and is preferred over a guess.\n'
        '2. **Provenance.** DOI, PMID or a stable dataset URL, plus where in the source the number sits '
        '(table, figure, section). A number without a locator does not count.\n'
        '3. **Validity range.** Species, tissue, temperature, method, and `n`. A value measured outside '
        'the stated resolution level must say so — this is the single most common reason a record is '
        'unusable for us.\n'
        '4. **Prefer a primary table over a summary.** We measured that a printed summary table is not '
        'a sufficient statistic for the measurements behind it: 24 pairs with identical printed summary '
        'differed by up to 0.5 N. If only a summary is available, say so explicitly.\n'
        '5. **If a whole dataset exists**, report its access route, licence and size, not just one '
        'number. A dataset is worth more than a value.\n'
        '\n## Rules\n'
        '- Do not fit, tune or infer the value from a model. This job looks for a measurement.\n'
        '- Do not report a range as a value; report the range with its definition.\n'
        '- Status `PENDING_INDEPENDENT_REVIEW`. No claim of biological validation.\n'
    )
    if terms:
        body.append('\n## Starting search terms\n' + '\n'.join(f'- {t}' for t in terms))
    return head + '\n'.join(body) + '\n'


def main() -> int:
    items = []
    if os.path.exists(TARGETS):
        d = json.load(open(TARGETS))
        items += [dict(x, _src='ACQUISITION_TARGETS') for x in d.get('targets', [])]
    if os.path.exists(HARVEST):
        d = json.load(open(HARVEST))
        items += [dict(x, _src='harvest') for x in d.get('wellformed_orderable', [])]

    made, skipped_private, skipped_thin, q = [], 0, 0, []
    seen = set()
    for it in items:
        blob = json.dumps(it, ensure_ascii=False)
        if BLOCK.search(blob):
            skipped_private += 1
            continue
        qty, unit = it.get('quantity'), it.get('unit')
        if not qty or not unit or len(str(qty)) < 12:
            skipped_thin += 1
            continue
        key = re.sub(r'[^a-z0-9]+', '', str(qty).lower())[:80]
        if key in seen:
            continue
        seen.add(key)
        note = next((v for k, v in STRESS.items() if k in str(qty).lower()), None)
        nid = it.get('node_id') or slug(str(qty))
        jid = 'BT-ANCHOR-' + slug(nid.replace('bodytwin:ACQ-', ''))
        d = os.path.join(RESULTS, jid)
        if os.path.exists(os.path.join(d, 'RESULTS.md')):
            continue
        os.makedirs(d, exist_ok=True)
        open(os.path.join(d, 'BRIEF.md'), 'w').write(brief(it, note))
        open(os.path.join(d, 'ALLOW_WEB'), 'w').write('1\n')
        json.dump({'id': jid, 'kind': 'acquisition_search', 'category': 'acquisition',
                   'claim_type': 'information_link', 'quantity': qty, 'unit': unit,
                   'stress_point': bool(note), 'source_list': it.get('_src'),
                   'review_state': 'PENDING_INDEPENDENT_REVIEW'},
                  open(os.path.join(d, 'JOB.json'), 'w'), ensure_ascii=False, indent=1)
        made.append((jid, bool(note)))

    # Stress points first: they decide a verdict we already hold.
    made.sort(key=lambda x: (not x[1], x[0]))
    slots = ('A', 'B', 'C', 'D')
    for i, (jid, _s) in enumerate(made):
        q.append(f'{slots[i % len(slots)]} swarm {jid}')
    if q:
        with open(QUEUE, 'a') as f:
            f.write('\n'.join(q) + '\n')
    print(f'anchor jobs created: {len(made)}  '
          f'(stress points: {sum(1 for _j, s in made if s)})  '
          f'skipped private: {skipped_private}  skipped without quantity+unit: {skipped_thin}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
