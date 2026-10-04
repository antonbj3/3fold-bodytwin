# Styrning LANE_EYE_OPTICAL_TWIN — efter r58

## You answered the decomposition, and one of the answers is a warning
The axial-length conversion carries part of the remainder: `original89_absolute_AL_slope_reduction_fraction =
0,45051175527931675` and MAE falls 0,36347869813018774 → 0,32601026339332045 D. I computed the gain
myself: **0,037468434736867284 D**, that is, 10,31 % relatively, and the count within 0,25 D goes 35 → 38 of 89.

**But the short eyes get worse.** `short20_paired_MAE_before_mm = 0,029999999999999895` →
`after = 0,05796520500000035`, a deterioration by a factor of **1,9321735**. A mean gain that almost
doubles the error in a subgroup is not a gain until the subgroup is reported.

## The obstacle, exactly
`physical_cause = NOT_IDENTIFIED` and `fixed_shift_span_D = 0,6553160292970219`. So you have a
correction that works on average, not a mechanism. And `population_significance = False` despite
empirical p = 0,0001 — which is the right abstention, since the assignment null is 0,0023351694417790836 D
and it tests the assignment, not the population.

## Changed operation
Split the 89 eyes by axial length into the groups your own correction treats differently, and report MAE
per group **before and after**, not just overall. If the correction harms a group, it must be conditioned on
the group or withdrawn. State the boundary in mm where it stops helping.

## Strongest control
The fixed shift alone, that is, `fixed_shift_span_D = 0,6553` applied without
axial-length dependence. If the groupwise conditional correction does not beat the fixed shift in every
group, the axial-length dependence is not load-bearing.

## Falsifier
If the deterioration in the short eyes remains under every conditioning boundary you test, then
axial length does not carry the remainder but another quantity that covaries with it — and then `physical_cause`
must remain NOT_IDENTIFIED and the correction must not be used.
