# Steering r19 — coverage has not moved in two rounds; raise it or close

## The obstacle, in your own numbers
R17 and r18 report **identical** counts: `silent_total` 181, `rows_with_computed_dimension` 55,
`CONDITIONAL_PHYSICAL_LIFT` 10, `not_reconstructed` 92, `pending_job_dispositions` 110. Two rounds without
a single new row is a state and not a coincidence. R18 instead worked on the precision of
the differential law through O(M) cohort weights, which is an internal improvement of something that already worked.

## The operation this round — one of two, choose explicitly
**A. Raise coverage.** 92 rows are not reconstructed and 110 dispositions are pending. Take thirty of the 92
and reconstruct them. Report new `rows_with_computed_dimension` and new lifts. Coverage is about
one third today; the target is for the number to move.

**B. Close with what you have.** Say what the 55 computed dimensions and the 10 conditional physical lifts
are useful for, what is required for the 92, and leave it as a named acquisition item. That is a clean outcome.

What is NOT an option is a third round with the same five numbers.

## Control and falsifier
- **Control:** a single reporter on the same held-out target without extra references — your own control row.
  Keep the fact that conventional equally informed moment and bin controls match; no algorithm requirement.
- **Falsifier:** if the thirty new rows give zero new lifts, the silent set is silent for a reason that
  is not rank deficiency, and that reason must then be named — it is a sharper finding than more lifts.
- **Forbidden:** another precision improvement of the differential law; reporting unchanged
  coverage numbers as progress.

## ADDENDUM r20 — your falsifier did not fire, and one measurement closes the nullspace
r19 chose A and moved the numbers: `rows_with_computed_dimension` **55 → 83**, `not_reconstructed`
**92 → 62**, `CONDITIONAL_PHYSICAL_LIFT` **10 → 18**, with 180 source forward checks and **0 validation
failures**. My falsifier was that thirty new rows giving zero new lifts would mean the silent set is
silent for some other reason than rank deficiency. You got **8 new lifts out of 30**, so rank deficiency
is a real cause in about a quarter of cases and the line is worth continuing.

Two things to carry into r20 beyond the next thirty:

1. **You found an exact sufficiency failure and should name it as one.** 100 grid pairs with
   `grid_summary_error_exact_binary64 = 0` and a transit gap of **exactly 1/300 h in every pair**. That is
   a summary that is identical in binary64 and still leaves a fixed downstream difference — the cleanest
   instance of the pattern we have, because the gap is the same rational number in all 100 cases. Say what
   the summary is and what the 1/300 h represents physically.
2. **The nullspace is one-dimensional and one measurement closes it.** `grid_regular_cell_rank = 3` with
   nullity 1, and `grid_rank_with_independent_transit = 4`. So an independent transit measurement removes
   the deficiency entirely. **Write that as an orderable acquisition item** — quantity, unit, what it
   decides, and the rank it restores — in `results/LANE_NULLSPACE_NAMES_IT/ACQUISITION_TARGETS_V1.json`.
   A one-dimensional nullspace with a named closing measurement is worth more than thirty more
   classifications, because it is actionable.

Keep the nominal-boundary rank as UNDEFINED with the step-artifact note. Do not let a finite-difference
artifact become a reported rank.
