# LANE_DETAIL_LAYER_EDGES — decisions outside the eye, from the detail layer

Result directory `results/LANE_DETAIL_LAYER_EDGES/`.

## Measured state, the night of 3–4 October
The readable repo `source_repository/` carries 313 MECHANISM documents with their own published
reference answers and their own falsifiers, and 9 688 scripts with numbers and units. Until tonight **one** of
the documents had a node in the constraint network; the harvest added 99 edges, and the network now has 155 edges
where 149 have evidence that `tasks/assembly/net_staleness.py` can read (6 are declared unresolved).

That entire layer is still **text**. `results/LANE_READ_CREATIVE/CONSUMABLE.json` lists 162
consumable observables from 55 papers, 118 of them outside the eye, and exactly one is consumed
(the colon's gas regime, `tasks/assembly/gas_regime_decision.py`).

## The obstacle
A document carrying a reference answer is not a decision. The eye chain took an entire night to get there and is
therefore the only part of the twin that chooses anything. What is missing is the same step outside the eye:
preoperative or pre-intervention quantities in, **the choice** out, measured against the source's own reference answer.

## The operation
Choose the observable whose source already carries BOTH a published number and a falsifier, and write an
executable file in `tasks/assembly/` that returns the choice, not the description of the choice. Two inputs are
the strongest, both verified by me tonight against the source file:

| input | source's number | state |
|---|---|---|
| meniscus, contact area and peak pressure at 1000 N | intact 1150 mm²/3 MPa, total meniscectomy 520 mm²/6 MPa, Baratz partial −10 %/+65 %, Rivarola FE 110±8 mm²/1,2±0,2 MPa | `MECHANISM_MENISCUS_LOAD_DISTRIBUTION.md` §2, ten network edges point to the table rows |
| disc, nucleus pressure from axial load | `P = k·F/A`, `k = 1,3–1,5`, `A = 1800 mm²`; in-vivo reference answer 0,53–0,65 MPa (Wilke 1999) | `MECHANISM_INTERVERTEBRAL_DISC.md`, already carries an executed FAIL: the disc model alone does not reproduce measured facet load |

## The strongest control
Equally informed: the same inputs, but the population mean instead of the relation. The decision
counts only if it beats the control on the source's own reference answer. `claim_type` determines the comparison — for
`information_link` the control is practice WITHOUT the measured number, not another method.

## The falsifier
Print before the run which outcome defeats the decision. Two traps from tonight apply here:
a composed standard deviation is not a floor (50 of 89 eyes lay below their own), and a number on
the lens plane is not the same number on the spectacle plane (factor 0,6864 in the eye chain). The source's flag is
not a gate: in Kim 2008 table 2, 26 of 28 deviations above 2 % are unflagged, largest unflagged 20 %.
Evidence must point to a file and line or JSON pointer and pass `net_staleness.py`.
