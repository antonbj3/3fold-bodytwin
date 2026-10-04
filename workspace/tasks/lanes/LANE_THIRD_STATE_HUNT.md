# LANE_THIRD_STATE_HUNT — find a TREDJE state in the same summary fiber

Results catalog `results/LANE_THIRD_STATE_HUNT/`. Decisions outside the eye.

## Why this lane exists
Six chains in the network carry a summary with identity error **exactly zero** next to a decision-relevant one
gap. `tasks/assembly/sufficiency_ledger.py` + border `T-E33` collects them, and the proof lane's border says that it
minimum number of extra fields is `⌈log_q M⌉` where `M` is the largest number of distinct decisions in **one** fiber.

The border is sharp. Our data is not: **each lane has exhibited exactly one witness pair**, so `M ≥ 2`
everywhere and the limit gives the same answer — a field — for all six. It cannot distinguish a chain that
hides two decisions from one that hides twenty.

## The goal, in order of priority
**The adhesion chain first.** `LANE_TISSUE_ADHESION_CONTACT r4` has
`sufficiency_summary_identity_error = 0,0` **and** `sufficiency_normal_field_identity_error = 0,0`
at the same time as `sufficiency_peak_force_ratio = 10,0`. A tenfold difference in peak power behind two
identical fields is the one of the six most likely to hide more than two decisions. Its
contact radius is 244,95 µm, tangent stiffness 14 691,34 N/m, and 16 of 48 initial
normal configurations are unstable.

Then the complement window, where the gap is **32,36×** the entire interval of the decision
(17,60797843701942 nM vs [0; 0,5441770553588867]).

## The operation
Construct a **third** state that lies in the same fiber — thus identical summary and
identical normal field — and which gives a **third** decision, separate from both known ones. Return
the condition, not the argument that it should exist.

If you succeed, `M ≥ 3` and the limit jumps to two fields for that chain. Then continue: each
additional states you exhibit raise `M` and make the book sharper.

## Strongest control
The two known conditions. A third that is only an interpolation between them does not count — it should
be a condition that yields a decision **neither** of the two yields.

## Forgers
If you can **prove** that the fiber holds exactly two decisions, then `M = 2` is fixed and not a lower
limit, and then one field is provably enough for that chain. It is as valuable a result as one
third state — specify which of the two you supply.

And the negative case, which the proof lane named: if the fiber accommodates an **infinite** family of decisions, exists
no final repair. Say it with the certificate, not with a guess.

## Rules
`PENDING_INDEPENDENT_REVIEW`, no statement of biological validation, nothing like
names an individual. Any published source with PMID or DOI.
