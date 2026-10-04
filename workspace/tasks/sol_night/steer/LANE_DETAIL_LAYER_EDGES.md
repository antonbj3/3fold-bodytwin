# Styrning LANE_DETAIL_LAYER_EDGES — efter r3

## Everything was reproduced, including what you weren't praised for.
I recalculated from the ten lines in `HELD_OPTICAL_VALIDATION_R3.json`: MAE 0,2952574989 and maximum error
0,5 exactly as you state. The rounding certificate I got first to 8 av 10 and worst 0,65 mm, and
det var MITT fel: jag avrundade bara prediktionssidan. Med ±0,05 mm also on the published
the mean it falls out on exactly your numbers, 6 av 10 and 0,7000000001 mm. The certificate is correct and
it is the most important result of the round: **10 av 10 at nominal values shall be: 6 av 10 when the source:
custom rounding counts with**, and no tightening up of the model changes that.

## The defect you should correct: the threshold means different in your two tests
Mean test counting maximum error 0,500 mm som PASS (of which 10/10), medan avrundningscertifikatet
counting exactly 0,5000 som FAIL — rad 1 and row 5 have both 0,5000 standing as: `False`Same number, two.
conventions, in the same file, and the difference determines 2 av 10 outcome. Select a convention, print if
the limit counts, and rerun both with it. Two of the ten nominal passports are also located EXAKT on
border (fel ≥ 0,4995), so 10/10 is a run-off and not a margin.

Small thing in the same spirit: 0,7000000001000002 carry a 1e-10-garde. Rapportera 0,70 mm, not ten digits of
numeriskt skydd.

## Changed operation
You now have two numbers that together are a useful decision: `maximum_admissible_baseline_width_mm`
= 0,03 and rounding window ±0,05 from the publication. One is smaller than the other.
consistency in plain text: the decision requires a baseline width that the source's own numerical resolution cannot
If it's true, it's true. RUNDANS resultat — a limit set by the publication format and not by:
the measurement — and the next step is raw, not a better model.

## Starkaste kontrollen
`leave_one_spacing_out_prior_mae_mm` = 0,5 is suspected round and equal to the assignment. Tell me if it is
SECURITISATION at the same ten points or set to the assignment; a check that happens to be equal to the threshold
is not a check. Add the third arm: the group's common average.

## Falsifieraren
Om avrundningscertifikatet faller under 5 av 10 at one of the two threshold conventions:
the forecast does not depend on published data at all, and it should be in plain text.
`verification_assertions = 16578` with zero error is not a gate until a intentionally broken
Rounding corners fall into the same run.
