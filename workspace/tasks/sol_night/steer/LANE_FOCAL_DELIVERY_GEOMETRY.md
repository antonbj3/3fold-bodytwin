# Styrning LANE_FOCAL_DELIVERY_GEOMETRY — efter r4

## The round overturns my own decision at a scale I did not consider
My `focal_delivery_decision.py` sets the minimum allowed point spacing to **3,130 µm**, calculated from two
independent measurement pairs that give the same diffusivity within 1,40×. Your round gives
`fast_branch_sign_distance_um = 0,004122347`, thus **4,1 nm**, and `gain_sign_flips_bool = True`.

The ratio is **759×**. The decision operates in micrometres while the SIGN of the gain reverses within four
nanometres of registration error. A decision whose direction is determined 759 times more finely than its own
quantity scale is not robust at the scale it claims to apply to, and that is the round's result.

You also honestly report that `practical_gain_gate = FAIL_FIXED_FAMILY_ALL_PHASES`: the targeted
design gives **0,43 to 0,63 %** against the target 10 %. And `nominal_registration_summary_gate = FAIL` with
`summary_identity_error_rad = 0,0` — identical summaries, different downstream outcomes, the seventh
time in the project.

## Changed operation
1. **State the registration requirement as a number in the decision's own unit.** If the sign reverses at 4,1 nm,
   the registration requirement is 4,1 nm, and the question is whether any measurement method reaches it. Say which and with what
   uncertainty, or say that none does — then the decision cannot be made and that is the result.
2. **Separate the two scales in the report.** The diffusion length 3,13 µm is how far the carrier travels.
   The sign boundary 4,1 nm is how close two paths are. They are different quantities and must never appear in
   the same table without saying so.
3. **The gain 0,43–0,63 % against 10 %** — say whether 10 % was an arbitrary target or derived. If it is
   arbitrary, replace it with what would actually change a choice.

## The strongest control
The uniform bath at the same total dose, which my decision already has, and which is IDENTICAL to point delivery
below the diffusion length. Every gain is calculated against it.

## Falsifier
If no registration method reaches 4,1 nm, the decision's direction is indeterminate, and the decision should then stand as
a necessary condition without direction rather than as a choice. Write that before the next round.
