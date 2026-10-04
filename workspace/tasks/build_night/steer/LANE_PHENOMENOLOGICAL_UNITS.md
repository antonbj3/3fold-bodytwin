# Steering r19 — the census is closed; now the sweep that r18 did not get to

## r18 closed the census correctly
Final counts, 1 491 families: **FYSISK 1 270, FENOMENOLOGISK 153, NORMALISERAD 18, UNKNOWN 50**. Of the
62 open families, 4 resolved physical, 8 phenomenological and **50 are `UNCLASSIFIABLE_WITHOUT` one
source-native reference definition**, so `census_disposition = CLOSED_WITH_ONE_ACQUISITION_ITEM` is the
right verdict and the acquisition item is the single thing standing between us and an exact fraction.

So the headline moved from a guess to a bound: **10.26 % measured, at most 13.62 % if all 50 fall
phenomenological.** Report it that way and stop widening it.

## This round is the sweep, and nothing else
My addendum reached you after r18 had started, so it was never run. It is now the whole task.

A swarm job compared our epithelial renewal front velocity to a held-out external measurement: ours is
**30.0 µm/h against a measured 9.00 ± 0.465 µm/h**, ratio 3.3333, `refutes_us: true`. The damning part is
not the factor but the provenance — the source states the constant was *"chosen to give 0.03 mm/h"*. A
constant back-solved from a desired output is the purest phenomenological debt there is, and **the census
did not flag it**.

1. **Find out why it was invisible.** Does the classifier look only at declared units and provenance
   fields, and never at whether a value was chosen to produce a target? Say which, in one sentence.
2. **Add the class `TUNED_TO_TARGET` and sweep all 45 cells for the signature.** Phrases next to an
   assignment of the form "chosen to give", "set so that", "calibrated to match", "to reproduce",
   "tuned", "adjusted until". Report how many constants carry it and in how many cells.
3. **The number I want is this one:** how many of our REPORTED OUTPUTS depend on a constant that was
   chosen to produce a desired value. Use affine sensitivity on each such constant and count the outputs
   whose value moves. That bounds how much of the model is self-referential, which nothing else we have
   measured does.
4. **And check the three surgical cells and the eye chain first**, because their numbers are the ones
   now being compared to clinical thresholds. A tuned constant inside a chain that reports 0.426 D would
   make that dioptre figure circular.

## Control and falsifier
- **Control:** the closed census, which reports 10.26 % and contains none of these cases.
- **Falsifier:** if no cell outside the front-velocity case carries the signature, it is one defect and
  not a class. Say so plainly — the sweep is then cheap insurance and the census stands as written.
- **Forbidden:** rerunning the census; searching for the airway reference definition (that is the
  acquisition item, not this round's work); reporting 10.26 % as exact while 50 families are open.
