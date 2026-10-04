#!/usr/bin/env python3
"""Fetch the references our own lanes already identified but never consumed.

The gap this closes. A lane cited 135 distinct DOIs in its own files while acquiring zero raw data
files: 31 dataset metadata descriptions and nothing behind them. Checked against everything we can
actually consume — the harvested anchor records and the literature collection — **131 of the 135 are
absent**. So our lanes identify external sources at a good rate and nothing harvests them. At least one
of the 131 is a data archive rather than a paper, which is the more valuable kind.

Why a fetch and not a search. The anchor searches so far were open-ended: name a quantity, find a
source. These jobs are the opposite and much cheaper — the source is already given, so the job is to
extract the values with their units and validity range from one named reference. An open search returned
33 partial results, 5 full hits and 5 verified absences out of 42; a fetch against a known DOI should do
better, and if it does not, that is worth knowing about the sources our lanes cite.

Ranking is by how many times our own files mention the DOI, on the assumption that a reference cited
eighteen times carries more of our reasoning than one cited once. That is an assumption, not a
measurement, and it is stated as such.
"""
from __future__ import annotations

import json
from pathlib import Path

W = Path('')
SRC = W / 'notes/UNCONSUMED_REFERENCES.json'
QUEUE = W / 'tasks/lanes/bt_queue.txt'
RESULTS = W / 'results'
TOP = 15


def brief(doi: str, mentions: int) -> str:
    return f"""# Extract the values from one named reference

**Reference.** `{doi}`

This reference is cited **{mentions} times** in our own lane files, and it is one of 131 that our lanes
identified while nothing harvested them. The source is given, so this is not a search. Open it and take
out what it actually measures.

## Deliver a machine-readable table first, then the prose

One row per value, exactly these columns, and the discussion after the table:

| quantity | value | unit | where in the source | validity | n |
|---|---|---|---|---|---|
| what was measured | the number alone | the unit alone | table or figure number | species/tissue/temp/method | sample size |

One value per row. Give a range as two rows, low and high, never as one cell. If the source reports a
quantity we would want but does not give a number for it, write the quantity with `NO_NUMBER` rather
than omitting it — an absence that is stated is usable and an absence that is silent is not.

## What matters beyond the numbers

1. **Is this a paper or a data archive?** If raw data is downloadable, report the access route, the
   licence and the size. A dataset is worth more than a value, because the next question can be asked
   of it too.
2. **State the validity range for every row.** Species, tissue, temperature, method, sample size. A
   value used outside its range is the single most common reason a record turns out unusable for us.
3. **Say whether the number is measured or derived in the source.** A value computed by the authors
   from a model is not a measurement, and it must not arrive here labelled as one.
4. **If the reference turns out not to support what our files cite it for**, say so plainly with the
   quote. That is a more valuable result than a table.

Status `PENDING_INDEPENDENT_REVIEW`. No claim of biological validation. Write `RESULTS.md` starting
with the job id, plus `results.json`.
"""


def main() -> int:
    data = json.loads(SRC.read_text())
    items = data['not_consumed'][:TOP]
    made = []
    for it in items:
        doi = it['doi']
        slug = doi.replace('/', '-').replace('.', '-').replace('(', '').replace(')', '')[:46].upper()
        out = RESULTS / f'BT-FETCH-{slug}'
        if (out / 'RESULTS.md').exists():
            continue
        out.mkdir(parents=True, exist_ok=True)
        (out / 'BRIEF.md').write_text(brief(doi, it['mentions']))
        (out / 'ALLOW_WEB').write_text('1\n')
        json.dump({'id': out.name, 'kind': 'reference_fetch', 'category': 'calibration',
                   'doi': doi, 'mentions_in_our_files': it['mentions'],
                   'claim_type': 'information_link',
                   'review_state': 'PENDING_INDEPENDENT_REVIEW'},
                  (out / 'JOB.json').open('w'), ensure_ascii=False, indent=1)
        made.append(out.name)

    slots = ('A', 'B', 'C', 'D')
    if made:
        with QUEUE.open('a') as f:
            f.write('\n'.join(f'{slots[i % 4]} swarm {j}' for i, j in enumerate(made)) + '\n')
    print(f'fetch jobs queued: {len(made)} of {len(items)} considered '
          f'({data["_meta"]["not_consumed"]} unconsumed in total)')
    for m in made[:8]:
        print('  ' + m)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
