# Steer, round 25: leave the spectral line and build the one connection that has data and no owner

## Why you are being turned
Round 24 gated FAIL with the obstacle stated as unacquired registered spectral stimulus and
same-specimen component response. That is an honest stop: the line cannot proceed without data that is
not on this machine. The CIE tables you acquired are real external facit with DOI metadata and they
stay; the spectral question moves to the swarm, which is the half of the division of labour that
fetches external references. Do not spend round 25 on it.

## The obstacle that is actually binding
The decision chain this lane feeds now makes two decisions scored against real outcomes, and both were
strengthened today: the implant power runs prospectively at 0.535 D mean miss against the implanted
power's 1.130 D on 89 eyes, and the toric cylinder holds against a second instrument at 0.2934 D
against 0.2817 D. What neither can do is say WHICH eyes it should not be trusted on. A forecast-
verification attempt failed at exactly that: Brier skill +0.038 and -0.085 against the cohort rate,
resolution 0.0034, i.e. the per-patient confidence carried almost no information.

The reason is identifiability, and it is already diagnosed: an additive stage decomposition of a total
observed only once is not identifiable, so a per-eye uncertainty cannot be fitted from outcomes alone.
It needs a per-eye measurement of the instrument's own spread.

## That measurement exists and nothing consumes it
`results/LANE_CORNEA_SHAPE/r3/REPEATABILITY_V1.json` holds 900 exams over 300 eyes, three per eye, with
the per-eye repeat range and sample spread for KrSEQ, Kf3rSEQ, Kf6rSEQ and the zone difference. Measured
from it: the repeat range exceeds the 0.25 D clinical threshold in 43 of 300 eyes and 0.50 D in 6.

## The operation for round 25
Propagate a per-eye instrument spread through the chain to a per-eye spread on the DECISION, not on the
prediction. Concretely: for an eye with a measured repeat spread in corneal power, what is the induced
spread in the recommended implant power on the 0.5 D manufacturing grid, and in the recommended toric
cylinder? The chain's sensitivities are already measured — 1.349 D/mm to lens position, and `dR/dP` per
eye is in `results/ASSEMBLY_IOL_DECISION/DECISION_V1.json`.

Deliver the fraction of eyes where the induced spread crosses a grid step, because that is the decision
changing, not the number moving.

## Strongest control
The cohort base rate. Always answering "decidable" is right in 257 of 300 eyes, so accuracy is
worthless here; the question is whether the eyes you flag carry a materially higher rate of crossing a
grid step than the cohort does. State the lift over the base rate, not the accuracy.

## Falsifier
If the induced spread crosses a grid step in a fraction of eyes indistinguishable from the cohort rate,
or if the flagged and unflagged groups have overlapping confidence intervals on that fraction, say so
and stop. A second falsifier is available and should be run: the 300-eye cohort is a DIFFERENT
instrument and cohort from the 89 surgical eyes, so if your propagation only works when the two are
treated as one population, that is a failure and must be reported as one rather than assumed away.

## Rules
Status PENDING_INDEPENDENT_REVIEW. No clinical framing, measurement and model level only. No
breakthrough without an equally informed control. Internal data stays on the machine.
