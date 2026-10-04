# Steering LANE_EXTERNAL_SOLVER — after r31

## The joint reading deteriorates monotonically, and that is the round's most important number
Calculated against the reserve 1,6204384e-07:

| reading | error | margin |
|---|---|---|
| marginal only (r31) | 1,394112e-07 | **+13,97 %** |
| individual preparation (r27) | 1,552537e-07 | +4,19 % |
| selected-then-frozen (r29) | 1,568603e-07 | +3,20 % |
| joint 90 % (r29) | 1,661171e-07 | −2,51 % |
| shared input (r31) | 1,801496e-07 | **−11,17 %** |
| joint + clock (r31) | 1,801500e-07 | −11,17 % |

The more caps treated jointly, the worse: from +14 % to −11 %. You set `pass: False` on
both joint readings and declare `caps_jointly_additive = False`, which is correct. But then
the joint reading is the main number and the marginal one is the footnote, not the reverse. Print the table
above in RESULTS so that no consumer picks +13,97 %.

## What dwarfs everything else
`exposure_kernel_identity_error = 0` with `exposure_kernel_downstream_gap_strain = −0,0444081`. Thus:
identical summaries, downstream difference **2,740e+05 × the reserve**. It is not a margin that
is missing, it is a summary that does not carry the decision at all — the same form as LASER_SURGERY's exact
witness tonight (identical sum, 0 against 4/3 Pa·s above a threshold of 2/3). Make it the round's heading.
Five orders of magnitude beat every percentage in the table.

## Two numbers that are not measurements
`shaped_causal_clock_sufficient_bracket_s = [2,7657604e-09; 2,7657610e-09]` has relative width
2,1551e-07, which is **1,81 ULP in float32** — numerical resolution, not sharpness. And
`arbitrary_clock_sufficient_bracket_s = [1e-21; 1e-20]` is zeptoseconds; you correctly flag that it
is neither necessary nor optimal, but a number at that scale should carry the word "unattainable" in the same
field, not in another key. `shaped_clock_cap_is_not_achieved_metrology = True` is exactly correctly
declared — keep that habit.

## What actually moved seed 1
`new_public_datasets_acquired = 1` with **38 fibril curves, 12 097 points and 38 geometry-matched
specimens**. It is the first time the lane has specimen data and not just certificates. The next round must
use them: fit the exposure kernel to the 38 curves and see whether the downstream gap 0,0444 shrinks when
the summary includes the geometry. That is the only route I see that addresses five
orders of magnitude.

## The strongest control and the falsifier
Control: the same 38 specimens, the same kernel, but geometry replaced by the population mean. The gain
is calculated only against it. Falsifier: if the gap does not shrink below 1e-3 strain with geometry included,
the summary form is wrong and must be replaced, not calibrated. `adversarial_max_timing_error_over_bound =
0,9579` means the bound has only 4,21 % slack under attack — good, and that number must remain.
