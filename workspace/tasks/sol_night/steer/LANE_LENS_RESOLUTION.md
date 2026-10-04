# Styrning LANE_LENS_RESOLUTION — efter r39

## The round's result is that the bound is valid and unusable
The triangle inequality gives `|meanA2 − meanA1| ≤ mean|P2 − P1| ≤ meanA1 + meanA2`, and the envelope covers
**448 of 448**. But its midpoint suffices in **0 of 448**. I calculated: the envelope's width is
0,003336755 to 0,012319794 rad while the largest observed P2 error is 0,001315245 rad, thus
the envelope is **2,54× to 9,37× wider** than the worst error it is supposed to enclose.

An envelope that always covers and whose midpoint never suffices is a correct bound that determines
nothing. That is the round's result and it belongs in the heading, not in a gate field.

## The bound must become narrower in a way that is not a guess
The triangle inequality is the weakest possible statement about two phase vectors: it uses only
the amplitudes and nothing about the angle between them. You have `finite_frequency_elemental_opposition_D_gate =
PASS_448_OF_448`, so the opposition condition holds in every case. If opposition is known or
bounded, the cross term is not free, and the envelope is then not the triangle inequality but something narrower.

Next round: derive the envelope WITH the opposition condition included. State how much narrower it becomes in the same
unit, and whether the midpoint then suffices in more than 0 of 448.

## The gate texts do not carry their polarity
`FAIL_0_OF_448` means zero PASS, which I could only determine from your prose in finalize.py line 61
("midpoint budget sufficient 0/448"). `FAIL_178_OF_448` cannot be determined from the string at all — 178
passed or 178 failed? Change the format to `verdict=FAIL passed=178 of=448`. This is the seventh
naming convention that has cost time in the project today.

## The strongest control
The triangle envelope as it stands. Every narrower envelope is calculated against it, and a narrower envelope that loses
coverage below 448 of 448 is not an improvement but a change of claim.

## Falsifier
Print before the derivation: which coverage you accept. If the opposition condition narrows the envelope but
coverage falls to, for example, 430 of 448, you have exchanged a valid unusable envelope for an
invalid usable one, and that must be said rather than reported as a gain.

## Unchanged since r33 and r35, and must be answered
`native_temperature_FWHM_K` and `nominal_FWHM_pm` remain without the effective wavelength having been
acquired. The two can only be converted using it, and 13,90896 pm corresponds to 0,154544 K only at
1137,94 nm. Either the fixture's wavelength is ~1138 nm and must be written out, or one of the numbers is wrong.
