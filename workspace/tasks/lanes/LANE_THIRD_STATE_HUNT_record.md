# LANE_THIRD_STATE_HUNT — find a THIRD state in the same summary fiber

Results directory `results/LANE_THIRD_STATE_HUNT/`. Decisions outside the eye.

## Why this lane exists
Six chains in the network carry a summary with identity error **exactly zero** next to a decision-relevant
gap. `tasks/assembly/sufficiency_ledger.py` + edge `T-E33` collects them, and proof_lane's bound says that the
minimum number of extra fields is `⌈log_q M⌉` where `M` is the largest number of distinct decisions in **one** fiber.

The bound is sharp. Our data is not: **every lane has exhibited exactly one witness pair**, so `M ≥ 2`
everywhere and the bound gives the same answer — one field — for all six. It cannot distinguish a chain that
hides two decisions from one that hides twenty.

## The goal, in order of priority
**The adhesion chain first.** `LANE_TISSUE_ADHESION_CONTACT r4` has
`sufficiency_summary_identity_error = 0,0` **and** `sufficiency_normal_field_identity_error = 0,0`
while `sufficiency_peak_force_ratio = 10,0`. A tenfold difference in peak force behind two
identical fields makes it the one of the six most likely to hide more than two decisions. Its
contact radius is 244,95 µm, tangent stiffness 14 691,34 N/m, and 16 of 48 initial
normal configurations are unstable.

Then the complement window, where the gap is **32,36×** the entire decision interval
(17,60797843701942 nM versus [0; 0,5441770553588867]).

## The operation
Construct a **third** state in the same fiber — thus identical summary and
identical normal field — that gives a **third** decision, distinct from both known ones. Return
the state, not the argument that it should exist.

If you succeed, `M ≥ 3` and the bound jumps to two fields for that chain. Then continue: every
additional state you exhibit raises `M` and makes the ledger sharper.

## Strongest control
The two known states. A third that is merely an interpolation between them does not count — it must
be a state that gives a decision **neither** of the two gives.

## Falsifier
If you can **prove** that the fiber contains exactly two decisions, then `M = 2` is established as more than a lower
bound, and one field is then provably enough for that chain. It is as valuable a result as a
third state — state which of the two you deliver.

And the negative case, which proof_lane named: if the fiber contains an **infinite** family of decisions, there is
no finite repair. Say so with the certificate, not with a guess.

## Regler
`PENDING_INDEPENDENT_REVIEW`, ingen utsaga om biologisk validering, inget excluded_category material, inget som
names the collaboration. Each published source with PMID eller DOI.
