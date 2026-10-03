"""Corrected brief: find a SECOND independent source for a quantity we already have.

Why the previous batch produced nothing, measured rather than guessed. 69 bounded jobs returned 0
values that any of our cells consume, and 15 of 15 jobs built from our own named gaps returned
VERIFIED_ABSENCE_OF_EVIDENCE. Fifteen identical answers is a result about the instruction, not about the
literature. Two defects in my brief caused it:

1. I made absence a PASSING answer -- "a stated absence is a result" -- which for a weaker model is
   always the cheapest exit. Search a little, declare nothing exists, be marked correct.
2. I ranked which references to fetch by how often they appear in our own files. The mentions came from
   a reference inventory I had built, so the ranking measured what was listed rather than what any cell
   needs. One extraction returned 898 well-formed records about body fat from circumference in soldiers.

And the pattern that DID produce value tonight was absent from all 69 briefs. Every real gain came from
finding a second independent source for a quantity we ALREADY hold: two instruments measured the same
eyes, and the toric result went from one device to three independent sources landing at 0.28, 0.29 and
0.31 D against practice's 0.63. The corneal model's own error was located the same way, 0.29 D against a
direct measurement of the quantity it assembles.

So the brief is inverted here. The job no longer hunts a quantity we lack. It takes a quantity we
compute and must find an independent way to obtain the SAME number, then report the disagreement.
Absence is no longer a passing answer on its own: a job that finds no second source must name the
nearest measurable substitute and where it is obtained, because that is what makes the finding
actionable instead of a dead stop.
"""
from __future__ import annotations

import json
from pathlib import Path

W = Path('.')
QUEUE = W / 'tasks/lanes/bt_queue.txt'
RESULTS = W / 'results'

# Quantities we already compute, each with the cell that computes it. The job must find a DIFFERENT
# route to the same number -- another instrument, another modality, another published cohort.
TARGETS = [
    ('CORNEAL-THICKNESS', 'central corneal thickness', 'um',
     'two devices in our own workbook disagree by 9.90 um with a spread of 6.26 on 89 eyes, i.e. a '
     'systematic bias; is there a third modality that measures the same thickness, and what does it '
     'disagree by'),
    ('LENS-POSITION', 'postoperative aqueous depth, i.e. where the implanted lens sits', 'mm',
     'we predict it to 0.114 mm from preoperative biometry and the chain converts 1 mm into 1.349 D of '
     'refraction; is this quantity measured directly by any other route, and to what repeatability'),
    ('POSTERIOR-CORNEAL-CYLINDER', 'posterior corneal cylinder', 'D',
     'two devices disagree by 0.384 D in magnitude with a systematic -0.343 D in the with/against-rule '
     'component; find a third independent source and its spread'),
    ('STROMAL-INDEX', 'refractive index of corneal stroma', 'dimensionless',
     'our chain uses one value from our own hydration model; find an independently measured value with '
     'its hydration state and temperature, because every corneal power in the chain scales with it'),
    ('IOL-POWER-LABEL', 'the true power of a labelled intraocular lens', 'D',
     'we treat the label as exact; manufacturing tolerance is published by standard and may exceed the '
     '0.5 D grid step, which would place a floor under every power decision we make'),
    ('AXIAL-LENGTH', 'axial length of the eye', 'mm',
     'one instrument only in our data; find an independent modality and its disagreement, because this '
     'is the largest single term in the refraction calculation'),
]

BRIEF = """# Find a SECOND independent source for a number we already have

**Quantity.** {quantity} [{unit}]

**Why this one.** {why}

## What counts as an answer, and absence alone does not

This job is NOT a search for a quantity we lack. We already compute this number. Your task is to find a
**different route to the same number** -- another instrument, another imaging modality, another
published cohort measuring the same thing -- and report what it disagrees by.

A result is one of these two, and nothing else:

1. **A second source, with its disagreement.** Give the value, the unit, the locator as a DOI or PMID,
   the population and method, and the spread. If it disagrees with ours, that disagreement IS the
   result and is more valuable than agreement.
2. **No second source exists, AND the nearest measurable substitute named.** If nothing independent
   measures this quantity, you must name what could be measured instead and where it is obtained.
   A bare "no evidence found" is NOT a result here and will be rejected. The previous batch returned
   fifteen such answers out of fifteen jobs, which told us about the brief rather than about the world.

## Required form

| quantity | value | unit | locator (DOI/PMID) | population and method | spread | agrees with ours? |
|---|---|---|---|---|---|---|

A bare URL is not a locator. Give a DOI or a PMID.

State whether each number is MEASURED or DERIVED in the source. A value the authors computed from a
model is not a measurement and must not arrive labelled as one.

Status `PENDING_INDEPENDENT_REVIEW`. No claim of biological validation. Write `RESULTS.md` starting
with the job id, plus `results.json`.
"""


def main() -> int:
    made = []
    for tag, quantity, unit, why in TARGETS:
        jid = f'BT-2ND-{tag}'
        d = RESULTS / jid
        if (d / 'RESULTS.md').exists():
            continue
        d.mkdir(parents=True, exist_ok=True)
        (d / 'BRIEF.md').write_text(BRIEF.format(quantity=quantity, unit=unit, why=why))
        (d / 'ALLOW_WEB').write_text('1\n')
        json.dump({'id': jid, 'kind': 'second_independent_source', 'quantity': quantity, 'unit': unit,
                   'category': 'calibration', 'claim_type': 'information_link',
                   'absence_alone_is_not_a_result': True,
                   'review_state': 'PENDING_INDEPENDENT_REVIEW'},
                  (d / 'JOB.json').open('w'), ensure_ascii=False, indent=1)
        made.append(jid)

    slots = ('A', 'B', 'C', 'D')
    if made:
        with QUEUE.open('a') as f:
            f.write('\n'.join(f'{slots[i % 4]} swarm {j}' for i, j in enumerate(made)) + '\n')
    print(f'queued {len(made)}: {", ".join(m.replace("BT-2ND-", "") for m in made)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
