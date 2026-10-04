# Steering LANE_DETAIL_LAYER_EDGES — after r4

## You answered everything, and two of the answers overturn your own result
1. **The boundary.** `deterministic_affine_error_crossings_cycle = [432,5110612; 365,0592860]` — the same two
   numbers I calculated by hand from your four points. And you are right to refuse to put intervals on them:
   `regression_residual_df_per_group = 0`, `NOT_ESTIMABLE_FROM_TWO_POINTS`. Two points per group give
   zero residual degrees of freedom, and that is the honest answer to my request, not a shortfall.
2. **The falsifier TRIGGERED.** I wrote that if the boundary falls below 250 cycles, the decision is not
   valid even on the available data. Your conditional joint-95 route gives
   `conditional_joint95_scaffold_first_crossing_cycle = 112,2594`. Thus: under that treatment,
   the decision does not support a decision at 250 cycles, which is the measurement point itself.
3. **The third arm changes the size of the gain.** `three_arm_mae_mm` = frozen_power 0,1980, pooled
   repair mean **0,255**, exchangeability 0,51. I recalculated the pooled one: the means 3,365
   and 4,585 give MAE exactly 0,255. The model's advantage is thus **0,057 mm against the right control**,
   not 0,312 mm against the exchangeability control — the latter is a straw man, because SWAPPING the groups
   is strictly worse than pooling them.

## The obstacle, now without circumlocution
The decision rests on published group means, and three independent gates point to the same wall:
the rounding window ±0,05 (r3), the non-estimable boundary (r4) and the joint-95 crossing at 112
cycles. No model improvement moves any of them. `physical_choice = ACQUIRE_LOCAL_HISTORY` is
the right conclusion and must now be pursued as the only track.

## Changed operation
Look for per-specimen data for the same preparation, not more summaries. You noted yourself that
CrackMNIST (DOI 10.5281/zenodo.18454958) is aluminium and only metadata — keep that discipline, it
is right. If no per-specimen dataset is openly available, that is the RESULT: write that the chain is limited
by publication format and not by the model, with the three numbers above as evidence.

## The strongest control
The pooled repair mean, 0,255 mm. The exchangeability control may remain but never as
the main comparison. Add a fourth arm if you want an upper bound on what group structure can give:
the one that sees the group label but not the cycle count.

## The falsifier
If per-specimen data arrives and the model's advantage against the pooled mean is below
0,057 mm on it, the group structure does not support the decision and that must be written plainly.
`deliberate_contract_faults = 19` is the right direction — keep them and report that they fail.
