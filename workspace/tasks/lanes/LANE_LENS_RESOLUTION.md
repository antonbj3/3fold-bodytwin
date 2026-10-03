# LANE_LENS_RESOLUTION — raise the resolution where the eye chain's error actually sits

Result directory `results/LANE_LENS_RESOLUTION/`.

## Why this lane exists, measured
The eye chain was compared to real post-operative refraction outcomes from a published dataset
(DOI 10.6084/m9.figshare.22736687.v1) in `LANE_EYE_OPTICAL_TWIN` r15. It lands inside the clinical
0.25 D tolerance in **7 of 20 cases** and inside 0.50 D in 12 of 20. That is not yet usable — but the
error is structured, and the structure points at one component:

| treatment | bias | MAE | RMSE | max error | within 0.25 D |
|---|---|---|---|---|---|
| thin lens | **−0.3865 D** | 0.5301 | 0.6591 | 1.2980 | 7/20 |
| thick lens | **+0.1455 D** | 0.4697 | 0.5902 | 1.4350 | 7/20 |
| index-convention control | +0.1643 D | 0.4754 | 0.5951 | 7/20 |

**The systematic bias is almost entirely a lens-model artefact: it is −0.39 D with a thin lens and
+0.15 D with a thick one.** The index-convention control reproduces the thick-lens figure, so the
difference is not a unit or convention artefact. And the cornea side was separately shown not to be the
missing piece: thickness, index, hydration and surface allocation together still leave 0.306 D
(EYE r14, all identity errors exactly 0.0).

So this is the one place tonight where raising resolution is measured to pay, which is the operator's
standing instruction applied where it bites rather than everywhere.

## Do this
1. **Replace the lens with a resolved one** — finite thickness, both surfaces, and a graded index rather
   than a single refractive index. Report the dioptre prediction per eye against the same 20 held-out
   rows, with bias, MAE, RMSE, max error and the count within 0.25 D and 0.50 D. The comparison must be
   the same rows, or it is not a comparison.
2. **Report the convergence order in the lens's own discretisation.** Number of surface samples and
   index layers against the change in predicted power. The coarsest sufficient setting is the deliverable
   — not the finest one you can run. If the prediction is invariant beyond a handful of layers, say so;
   that is a cheap and valuable result.
3. **Separate the two error sources explicitly.** How much of the remaining MAE is the lens and how much
   is everything else? Hold the lens fixed at its best setting and vary the rest, then the reverse. Two
   numbers.
4. **Then state what is still missing and in what unit.** If the residual needs a measurement we do not
   have, write it to `results/LANE_LENS_RESOLUTION/ACQUISITION_TARGETS_V1.json` with quantity, unit and
   what it decides — a clinical lens-geometry measurement is routinely made, so say exactly which one.

## Control and falsifier
- **Control:** the thin-lens treatment at its current resolution, on the same 20 rows. That is an equally
  informed control by construction — same data, same chain, only the lens representation differs — which
  is why a win here can be claimed at all.
- **Falsifier:** if the resolved lens does not improve the count within 0.25 D beyond 7 of 20, the lens is
  not the binding component and the error lives in the geometry we do not measure. That would be the more
  important finding and must be reported as the headline, not buried.
- **Forbidden:** fitting any lens parameter to the 20 rows — they are held out, and a fitted prediction is
  not a prediction; reporting a mean error without the within-tolerance count; changing the cornea model
  in the same round, because then the two effects cannot be separated.

## Delivery
`PORT.json` with the five error statistics against the same 20 rows for both lens treatments, the
convergence table with the coarsest sufficient setting, the two separated error contributions, and any
acquisition item in the standard schema.

No internal data. Everything PENDING_INDEPENDENT_REVIEW.
