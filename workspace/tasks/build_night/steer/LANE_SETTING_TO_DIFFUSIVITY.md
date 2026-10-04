# Steering LANE_SETTING_TO_DIFFUSIVITY — after r7

## You established that the coupling cannot be approximated away
`forcing_only_bound_violations = **178**` against `spatial_max_numerical_bound_violation_K = **0,0**` and
`modal_max_numerical_bound_violation_K = 2,38e−14`. Thus: the spatial, coupled treatment violates
the bound **never**, while treating the forcing in isolation violates it **178 times**. That is a clean
result and it closes a shortcut.

And isolation is more costly than geometry: `isolated_summary_max_gap_K = **33,98863138521645**` against
geometry's 27,614871866185098 K in r6 — ratio **1,2308**. Isolation is thus worse than not
measuring the geometry at all.

## The obstacle
Mode acceptance is 35 of 64 in thirteen of sixteen records, with two outliers at 36 and 37 — thus
**54,7 %** of the modes are needed, and `minimum_basis_proved = False`: the smallest passing basis is 8 but
not proved minimal. The localised basis differs from the axial one by **5,893052747634478 K**,
so basis choice carries almost twice as much as the drainage question did.

`held_max_center_residual_K = 34,954280749785894` across 24 held-out records with 8 overlaps. That is the same
order of magnitude as the isolation gap, so held-out data does not yet distinguish the treatments.

## Changed operation
Prove the minimality of the basis, or state why it cannot be proved. One number: the smallest basis that
meets the bound with proof. That determines whether 54,7 % of the modes is a property of the problem or of
the construction.

## The strongest control
Forcing alone, which gives 178 violations. Every basis or mode saving must be measured against how many violations it
allows, not against how few modes it uses.

## Falsifier
If the held-out centre residual of 34,95 K cannot be pushed below the isolation gap 33,99 K, then
held-out data does not distinguish the two treatments, and the 178 result is then purely numerical and not empirical.
Say which.
