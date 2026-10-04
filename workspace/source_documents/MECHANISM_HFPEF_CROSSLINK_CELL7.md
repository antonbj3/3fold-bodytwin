# HFpEF cell 7 — stiffness per unit collagen: the last named route

Cell: `data/hfpef_diastolic_v1/`. Pre-registration: `PREREG_CELL7.md` (+2 addenda), committed before
any cell-7 anchor was in hand. Code: `scripts/tissuetwin/hfpef_cell7_crosslink.py`.
Evidence: `cell7_evidence.json`. Anchors: `anchors/cell7_*.json` (3 files, 19 records).

## What was left to explain
Cell 1 left a passive **material** stiffening of **3.058x** in preserved-EF after removing both
geometric levers. Cell 2's G5 bounded the collagen **content** route at **2.300x** — the largest
reading any preserved-EF record states, table-grade, and generous by construction since it is the
ceiling at `w = 1`. So

    X = 3.058 / 2.300 = 1.330x

must come from **stiffness per unit collagen**. Three of the four named routes were already measured
shut: titin isoform (cell 2 G4), titin phosphorylation (cell 5), collagen content (cell 2 G5). This
is the fourth.

## The verdict: `ANIMAL_ONLY_DOES_NOT_CLOSE_THE_HUMAN_BUDGET`

| | count |
|---|---|
| records harvested across 3 families | 19 |
| clear G1 (a stated number, not narrative) | 2 |
| clear G1 **and** G3 (normalised to collagen) | 1 |
| of those, in the **stiffness** family (the only family that can carry M) | 1 |
| of those, **human** | **0** |

The one surviving record is **Norton 1996** — streptozotocin-diabetic **rat** LV, crosslink
fluorescence per µg hydroxyproline **1.667x**, stiffness slope **1.741x**, both collapsing toward
control under aminoguanidine, with collagen content verbatim confirmed unchanged. It is the only
record in the searched literature that varies crosslinking and measures the mechanical response at
constant collagen content. `1.741x > X = 1.313x`, so the mechanism is **sufficient in a rodent** —
and that is exactly as far as it goes. `n_eff = 1`, one lab, one species, not HFpEF.

No human record pairs a crosslink number with a mechanical number in the same hearts. The comparison
the cell exists to make **cannot be made in humans**, and that is the result.

## Two defects in my own driver, both caught after it first ran green

The first run returned **`CROSSLINKING_CAN_CLOSE_THE_BUDGET`** on a single admitted record. Both
halves of that were wrong:

1. **A dose was pooled with a response.** The admitted number was Beaumont's CCL ratio (3.0x) — a
   crosslink **content** ratio, insoluble:soluble collagen. That is the *dose*; M is the stiffness
   *response* per unit collagen. Two different quantities sharing a units-free "x" were summarised
   together — the **fourth** time this cell family has been caught doing exactly that. Families are
   now partitioned by source file and only the stiffness family can decide M.
2. **`+/-` is a plus-minus.** The reader accepted `+` and `±` only, so `47+/-4 vs control 27+/-3`
   matched nothing and Norton — the one genuinely usable record — was excluded as
   `NARRATIVE_NO_NUMBER`. The fail-**closed** twin of the fail-open defect cell 2 was caught in.

Fixing both flipped the verdict from a category error to an honest negative.

Separately, the driver's pool fell back to animal records when no human record survived
(`pool = admitted_h or rows`), which would have let a rodent multiplier decide a human budget with
no verdict string saying so. Removed; recorded in PREREG addendum 2, which is labelled **post-harvest**.

## P1 — what the excluded records measured (POST-HOC, not a test)

None of these clears G1+G3, so none enters M. They are reported because a measured direction that
contradicts the route under test should not be quietly dropped:

- **LeWinter** — the most direct human LV test found: CML **45.6 / 45.8 / 49.3 per µm²** across
  control / HTN / HTN+DM, **ANOVA P = 0.85**, and reported as *intracellular* with minimal
  association with ECM collagen;
- **van Heerebeek** — AGE deposition higher in **reduced-EF** diabetics (24.1 vs 8.8, P = 0.005)
  than in normal-EF diabetics (15.7 vs 8.2, P = NS) — the opposite of a preserved-EF-specific AGE
  route;
- **Gunja-Smith** — the only per-collagen pyridinoline number found **falls** with disease,
  2.07 → 1.0 mol/mol collagen, opposite in sign to the increase-in-crosslinking records, none of
  which carries a magnitude.

Every genuinely HFpEF-labelled crosslink record in the harvest reports the quantity as narrative
("enhanced", "increased P<0.001") or as a correlation coefficient. Not one states a magnitude with a
denominator.

## Where this leaves the residual
Four named routes, four negatives — three measured shut, one unmeasured in humans. The 3.058x is
real, geometric explanations are exhausted, and no molecular route currently carries a number that
closes it. **The single measurement that would decide it:** a human LV myocardial specimen with a
collagen-normalised crosslink assay *and* a passive stiffness measurement at a stated sarcomere
length, at constant collagen content — Norton's design, run on human tissue.
