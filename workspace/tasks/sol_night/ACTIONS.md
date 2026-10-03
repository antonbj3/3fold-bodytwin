# Action queue

The night log is a journal: a row there proves that something was found, not that something was done. That led
to findings being reported to the operator instead of worked on, which he pointed out on 2026-10-03.
This file is the queue. A row moves to DONE with the number that shows it is done, never on a
claim that it is done.

## OPEN

| # | what | evidence needed to close | owner |
|---|---|---|---|
| A3 | 162 consumable observables (1 consumed: the colon's gas regime, see D9) from 55 papers, 118 outside the eye, in `LANE_READ_CREATIVE/CONSUMABLE.json` | through the unit and range gate, then admitted edges with numbers in the text | me |
| A4 | 13 orderable measurement targets and 7 dataset routes from the same round | those that can be retrieved without ordering anything: retrieved and consumed | the swarm |
| A6 | the detail layer: harvest DONE — 312 inventoried, 99 edge proposals, 3 082 external references, 244 excluded records removed | edge proposals through the unit and range gate, then admitted; reference list crossed against our quantities | me |
| A7 | `native_glucose.py` using `np.trapz` which disappeared in numby 2.4.1; `myofascial_transmission.py` requires opensim | executable import, or a line stating that the dependency does not exist | jag |

## KLART

| # | what | the number that closed it |
|---|---|---|
| D1 | `from_vector` rotated every axis 90° | round trip exact across 36 axes; toric 0,2817 → 0,2808 |
| D2 | the net's tracked copy 34 edges behind | 56 of 56 in both copies, sync in the hourly cycle |
| D3 | the queue 92 % completed work | 13 919 → 1 047 rows, pruning in the cycle |
| D4 | Sol lanes without web search | `tools.web_search=true` in all three drivers |
| D9 | first harvested observable consumed into a decision | gas-regime decision outside the eye: threshold 76 mL/6h separates two regimes, single-coefficient control 77,8 % high in every row |
| D8 | three edges reported as UNCHECKED | all three already carried `evidence_unresolved` with a note; the error was in my check that did not read the field. Now **UNRESOLVED_DECLARED** (searched for, does not exist) is distinguished from **UNCHECKED** (unchecked), and UNCHECKED is zero |
| D7 | the laser refutation was mine, and it was wrong | quantities, denominators and observation operators differ; and I wrote the unit **µm⁻¹** where the source table says **µM**, thus read a concentration as an inverse length. The spread 2033× is real and does not affect the witness. Edge back to UNKNOWN, original annotation preserved verbatim |
| D6 | Q005 stored result 1000× wrong in both concentrations | the code's own run gives Met 0,20734 and Cr 0,31685 against stored 207,342 and 316,849; CL_Met identical 29,5852, thus only the concentration conversion. Old preserved as `.stale_20261003_factor1000`. Assay constant untouched. **The cell's own `unit_check` still passes** — the gate does not check what was wrong |
| D5 | memory reservation double the consumption | 1500 → 1100 for our rows, 1500 for dental's after their p95 1319 |
