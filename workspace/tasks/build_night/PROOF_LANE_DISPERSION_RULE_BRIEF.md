# PROOF_LANE_DISPERSION_RULE — what do you do with a published ± that doesn't say what it is?

## The trap, and what it has cost us in terms of speech
A published `±` with no stated convention has ruled here three times, and I've counted each
time itself:

- **mPTP ceiling.** With `±` read as SD, the separation **2,622950819672131** is pooled SD and that is enough
  with **3 animals per group** at 80 % strength. Read as SEM at n = 7 is true SD
  **10,583005244258363 pp**, the separation falls to **0,9913822242864975** SD and the group size
  rises to **16** — factor **5,33**. The decision itself carries the flag
  `dispersion_convention_unresolved = True`.
- **My own distribution gate.** I used per-animal SD as the uncertainty in a **group mean**
  and thus reversed the mPTP arms in both reads. The same trap, inside the tool that would catch it.
- **The coefficient of friction.** `0,295 ± 0,056` static and `0,255 ± 0,086` dynamic — you pointed out that
  the dispersion **not** is a rigorous uncertainty limit, so σ distances calculated from it are descriptive.
  I had counted 1,69–1,88 σ out of it and booked it as if it bound something.

## The question
Given a published central value and a `±` **without** stated convention, what is
decision theory correct to do? I want a rule I can code, not a recommendation.
Particularly:

- is there a **reading that is safe under both conventions** — that is, a decision that gives the same
  outcome if `±` is SD as if it is SEM — and what condition determines when it exists?
- when it does not exist: which of the two readings is the **conservative**, and is "conservative"
  even well-defined when the decision can fall to both hold-outs (which mPTP's does)?
- how many samples `n` are needed to **determine the convention from the data** instead of from the text, if
  the published speeches contain enough information;
- and a criterion for when a `±` should be treated as **no information at all** about dispersion, which
  is the position of the coefficient of friction.

## Strongest control
To always read `±` as SEM, which always gives the larger true the dispersion and thus the
most cautious conclusion. It is always sure and discards every decision that is true but close.
Indicate what it costs, in rejected decisions, on the three cases above.

## Falsifier
Construct a case where your rule produces a different outcome than both conventions do individually.
If you succeed, the rule is not a reading of the data but a third convention, and then it should be rejected.

## Rules
`PENDING_INDEPENDENT_REVIEW`, no biological validation statement, no excluded_category material, nothing like
names an individual. The endpoint rule applies to any precision.
