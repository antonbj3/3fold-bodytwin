# Styrning LANE_SPINE_LIGAMENT_STIFFNESS — ny lane 2026-10-04

## Five edges in the network stand UNKNOWN for the same reason
HARVEST-E0057 to E0061. The conditions say it plainly: `linear_stiffness` for ALL 350, PLL 250,
LF 200, SSL 100, ISL 60 is "a representative CHOSEN force-per-unit-strain coefficient; not N/mm and
not an extracted published measurement"Source: `source_documents/MECHANISM_SPINE_LIGAMENTS.md`.

The document itself records two things you must read before calculating:
- line 240: the numbers are NOT from Pintar et al. 1992's own table.
- line 265: behavior is controlled by geometry, by the placement of attachment points, and NOT by
  the `linear_stiffness` magnitude. "A stiffer or softer spring at the same geometric placement..."

## What I measured and what you should build on
Converted to N/mm with reference lengths I ASSUMED and that you must fetch from the document instead:
ALL 14,0 against published 33,0 ± 15,7 N/mm (ratio 0,42), PLL 10,0 against 20,4 ± 11,9 (ratio 0,49).
An earlier review found ISL 4,00 N/mm against measured 42 and 51 N/mm, PMID 26726784, thus 10–13×.
The direction is consistent: the chosen numbers are too soft. MY reference lengths are guessed, however —
25/25/15/15/20 mm — and the whole conversion depends on them. Fetch the real ones from `setSlackLengthFromReferenceStrain`
and the slack-length section (lines 132–137) before trusting the ratios above.

## The task is not to replace numbers. It is to test the document's own claim.
Question: is the outcome sensitive to the `linear_stiffness` magnitude or not?

Run the model with the chosen numbers and with the published ones, same geometry, same pose, and measure what
changes in the outcome the chain actually consumes. Two possible results, both valuable:
- The outcome is INSENSITIVE. Then the document's claim is confirmed, the five edges can be closed as
  "magnitude not load-bearing" instead of UNKNOWN, and no measurement is needed. That is a result.
- The outcome is SENSITIVE. Then the chosen numbers are a defect, the direction is known (too soft), and the five
  edges should carry published values with a locator.

## Strongest control
Same geometry with the chosen numbers. That is what the chain uses today, and the gain is counted against it.

## Falsifier
Write before running which outcome change counts as sensitivity. A number, not a judgment.
If the outcome changes less than its own numerical resolution, it is insensitive regardless of how large
the stiffness differences are.

## Rules
Unit in every field name: N per unit strain and N/mm are different quantities and have already been confused here.
The reference length used for conversion must be written out. DOI or PMID for every published
number. Everything PENDING_INDEPENDENT_REVIEW.
