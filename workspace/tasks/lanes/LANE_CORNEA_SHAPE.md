# LANE_CORNEA_SHAPE — the one quantity that decides the direction of a surgical outcome

Result directory `results/LANE_CORNEA_SHAPE/`.

## Why this is the highest-value construction left
The assembled twin now makes the surgical decision. On 20 patients it chooses the implant power and
lands closer to the hindsight-correct power than the implanted lens in 11 cases against 2, mean miss
0.688 D against 1.070 D. On 69 toric patients, using the **measured** posterior cornea instead of the
population estimate that practice uses cuts predicted residual astigmatism from **0.626 D to 0.282 D**
against measured refraction, better in 50 of 69.

What the decision cannot yet do is predict the **direction** of a corneal procedure's effect, and the
reason is measured precisely:

| measurement | result |
|---|---|
| two corneas with identical mean thickness, index field, hydration AND surface allocation | differ by **0.3056 D**, 4 of 4 pairs, all identity errors exactly **0.0** |
| span of the refractive change over the surface-allocation parameter alone | **0.8724 D** |
| the same span measured as a sufficiency gap | 0.87239 D — the two coincide to four decimals |
| rotationally symmetric control with the same mean thickness | **−1.36e-08 D**, so the entire effect is non-symmetric structure |

So thickness, index and hydration are jointly insufficient however finely they are resolved, and the
missing quantity is the **surface shape field with a common coordinate registration** — both corneal
surfaces as a height map in one frame, not a thickness summary.

**And unlike every other missing measurement found tonight, this one is routinely acquired**: a corneal
topographer measures exactly it. The clinical dataset already in use carries anterior and posterior
radii with their axes per patient, which is a two-parameter projection of that field.

## Do this
1. **Represent both surfaces as height fields in one frame**, with the registration explicit. Start from
   what the dataset gives — anterior radii and axis, posterior radii and axis, central thickness — and
   state exactly which degrees of freedom of the true field those four numbers fix and which they leave
   free. That statement is the deliverable even if nothing else lands.
2. **Then measure how much of the 0.306 D gap the dataset's projection closes.** Two corneas that agree
   on all four dataset numbers: how far apart can their refraction still be? That is the residual
   ambiguity of current clinical input, and it bounds what any calculator built on radii can achieve.
3. **Carry it into the decision.** Re-run the power and toric decisions with the shape field in place of
   the mean-curvature scalar, on the same 89 patients, and report the same statistics: power miss
   against hindsight, residual cylinder against measured. If the numbers do not improve, the shape field
   is not the binding quantity for these particular decisions and that is worth knowing precisely.
4. **And state the direction result.** With the shape field, is the sign of a simulated procedure's
   refractive change determined? The twelve-condition sweep flipped sign between +0.4658 D and
   −0.8618 D; say whether the flip survives, and if it does, name what else is missing.

## Control and falsifier
- **Control:** the mean-curvature scalar cornea as the chain uses today, on the same patients. Equally
  informed by construction — same data, same chain, only the corneal representation differs.
- **Falsifier:** if the dataset's four numbers already fix the refraction to inside 0.25 D, then the
  shape field adds nothing for these decisions and the 0.306 D gap lives in degrees of freedom that
  clinical input does not vary. Report that plainly; it would redirect the whole eye track.
- **Forbidden:** fitting any shape parameter to the measured refractions — they are the held-out
  outcome; simulating a procedure on a rotationally symmetric cornea; reporting an improvement without
  the same-patient control.

## Delivery
`PORT.json` with the registration statement, the degrees of freedom fixed and left free by the four
dataset numbers, the residual refractive ambiguity they permit, the re-run decision statistics on the
same 89 patients against the scalar control, and the verdict on sign determinacy.

No internal data. Everything PENDING_INDEPENDENT_REVIEW.
