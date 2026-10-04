# Steering LANE_DISTRIBUTION_GATE_43 — new lane 2026-10-04

## The state, from our own RESULTS_2026-10-03.md §5
Of **43 cells**, **21 pass a scalar point-or-band test**, and only **7 carry any
distribution statistics at all**. Eight had no distribution check. When one was built and run
on the eight, **3 of 8 verdicts reversed** — they pass their current gate and fail against the distribution.

The case that started it: a hemostasis cell passes its point-in-band test while **7,305 %** of its
predictive mass never closes, stable at 0,069 across time steps 0,020, 0,010 and 0,004. Stable
across three steps means structural, not numerical.

Thus: eight are tested, thirty-five are not. If the 3-of-8 frequency holds, about five
reversed verdicts await among the remaining scalar-gated ones. That is an expectation, not a finding.

## The task
Extend the distribution test to all 43. One row per cell: cell name, gate type (scalar point, band,
distribution), scalar verdict, distribution verdict, reversed yes/no, and the fraction of predictive mass that
does not close if it can be calculated.

A reversed cell is booked as **FAIL**, not as nuance. That is the whole point of the test.

## The tool exists
`tasks/assembly/distribution_gate.py` is generic: `resample_verdict(point_entries, dispersioner,
verdiktfunktion)` resamples and returns agreement fraction against the point verdict, requirement 95 %. You only need
to write the verdict function per cell. It is ten lines and has already run on seven decisions.

## The trap that must be avoided, and it was committed in the gate itself
The first version of my gate set **per-unit spread as uncertainty on a GROUP MEAN**
and wrongly failed a decision in both readings. The right rule:
- the source's ± is an SD → uncertainty on the mean is SD/√n
- the source's ± is already an SEM → uncertainty is the stated value itself

For the mPTP case that decided everything: `verified_dispersion = SEM`, n = 10, true SD = 4 × √10 = 12,649.
Read every cell's dispersion measure before setting any spread. If it is not stated in the source, the cell is
UNDECIDED, not passed.

If the 3 of 8 that already reversed have the same mix, some of them may **not** be reversed. Retest them
with the correct scaling before building on their outcomes.

## The strongest control
The cell's current scalar gate. It is what the chain uses today, and every reversal is counted against it.

## Falsifier
Write out before the run how many reversals you expect out of 35. If you find zero the test is too
weak — the known 3-of-8 frequency says zero is wrong. If you find 35 the threshold is too strict.

## Rules
The dispersion measure in every field name: `_sd`, `_sem` or `_unknown`, never just a number. DOI or
PMID where the source is published. Everything PENDING_INDEPENDENT_REVIEW.
