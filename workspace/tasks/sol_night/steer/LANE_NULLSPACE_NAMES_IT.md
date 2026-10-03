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
