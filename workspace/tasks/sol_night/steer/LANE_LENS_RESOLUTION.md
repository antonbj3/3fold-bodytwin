# Steer, round 9: build the total error budget for the decision chain

## Why you are being turned
Round 7 reported its own result honestly: 72 asymmetric diagnostic pairs flip the sign of the effect
despite identical phase-gradient means, and "no clinical improvement shown: still 7 of 20 within
0.25 D". The obstacle is unacquired registered kernel and irradiance data that is not on this machine.
The sufficiency finding stands and stays; the line cannot go further tonight.

## What is needed instead, and why it is the highest-value thing left
Tonight produced three results on the surgical decision chain and they are currently separate facts:

- The power decision runs prospectively: mean miss 0.535 D against the implanted lens's 1.130 D on 89
  eyes, closer in 36 against 22. `results/ASSEMBLY_PROSPECTIVE_POWER/PROSPECTIVE_V1.json`
- The toric decision holds across three independent sources for its key input: 0.2817, 0.2934 and
  0.3112 D against practice's 0.6264. `results/ASSEMBLY_TORIC_CROSS_DEVICE/CROSS_DEVICE_V1.json`
- The lens label is permitted to deviate by more than the grid step we choose on, in 89 of 89 eyes.
  Section 13 of `RESULTS_2026-10-03.md`

Separately they are three findings. Added up they are one checkable claim: how much of the remaining
error is irreducible given the inputs, and how much is still the model's. Nobody has computed that.

## The operation
One error budget for the prospective power decision, in dioptres at the spectacle plane, one row per
input, each row carrying where its uncertainty number came from:

  lens label tolerance        ISO 11979-2 scheme, and the measured deviation (0.18 to 0.24 D)
  lens position prediction    0.114 mm MAE x 1.349 D/mm, both measured here
  corneal power               between-device disagreement, measured on the same 89 eyes
  axial length                NO SECOND SOURCE EXISTS; carry it as unmeasured and say so
  corneal model assembly      0.292 D overstatement against a direct measurement of the same quantity

Compose them twice and report both: in quadrature if independent, and summed if they are not. The
difference between those two numbers is the value of the independence assumption, and it is not ours to
assume — `graph_engine.stage_decorrelation_verifier` is installed and `covariance_aware_margin` returns
exactly that bracket. Use it rather than writing your own.

Then the number that matters: our measured mean miss is 0.5345 D. How much of it does the budget
account for? If the budget exceeds the miss, something is double-counted and you must find it. If it
falls well short, the gap is the model's own error and that is where remaining work belongs.

## Strongest control
The measured miss itself, 0.5345 D on 89 eyes, is the facit. A budget that cannot be compared to it is
not a budget. The secondary control is the implanted lens at 1.130 D: any claim that our error is
irreducible must also explain how the surgeon's was twice as large.

## Falsifier
If the budget's independent composition exceeds the measured miss, the budget is wrong, not the
measurement, and you say so. Also: axial length has no second source, so any budget that produces a
confident total without it is overclaiming. State the unmeasured term explicitly and give the total as
a bound, not a value.

## Rules
Measurement and model level only, no clinical framing. Status PENDING_INDEPENDENT_REVIEW. No
breakthrough without an equally informed control. Internal data stays on the machine.
