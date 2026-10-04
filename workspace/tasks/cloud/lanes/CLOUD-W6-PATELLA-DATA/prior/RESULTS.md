# RESULTS — CLOUD-W4-PATELLA-RANGE


## Preregistered criterion (PREREG.md, sha256 9450d437…a6e, verified unchanged)

| Criterion | Result |
|---|---|
| Extrema and virtual-work derivatives agree within 1% on 0–120° | **PASS** (worst relative difference 1.5e-9) |
| Sensitivity separates endpoint and interior peaks | **PASS** |
| Empirical peak reproduction | **UNKNOWN**: there are no digitized, matched curves |

Interpretation fixed before viewing the sweep results: "peak angle within 1%" means within 1% of the 120° range, i.e. 1.2°.

## Sources and what they support

- Krevolin JL, Pandy MG, Pearce JC. *Moment arm of the patellar tendon in the human knee.* J Biomech 2004;37(5):785–788. DOI **10.1016/j.jbiomech.2003.09.010**.
  - **Direct access failed.** The egress proxy returned HTTP 403 for doi.org, api.crossref.org, pubmed.ncbi.nlm.nih.gov (PMID 15047009), Europe PMC, OpenAlex and findanexpert.unimelb.edu.au. I did not fetch the verbatim abstract.
  - What I saw came only from a web-search result summary: peak patellar-tendon moment arm of 4–6 cm across 6 cadaver knees, maximal near 45° flexion, measured about the finite screw axis of the tibia relative to the femur. Treat this as secondary and unverified text.
  - The abstract supplies no digitized curves, and I obtained none. **No empirical values are reproduced or compared.**

## Model (original, `patella_range.py`, Python 3 standard library only)

The model is sagittal-plane, in millimetres and degrees, in a femur-fixed frame.

- **Condyle and tibia.** The femoral condyle is a circle of radius R. The tibia flexes by θ. Seen from the tibia, the condyle centre rolls back by s·R·θ, where s is the rollback fraction: s = 0 is a hinge and s = 1 is pure rolling.
- **Tuberosity.** The tibial tuberosity T is fixed in the tibia frame.
- **Patella.** The patella apex P lies on a trochlear circle (radius Rp, centre (cx, cy)). Its position is set by an inextensible patellar tendon, |P−T| = L0.
- **Quadriceps.** The quadriceps runs as a straight line from a fixed point Q to P.

**Baseline parameters** (chosen a priori, not fitted): R 22, s 0.5, T (32, −40), Rp 30, trochlea centre (8, 12), γ0 35°, Q (10, 250).

**Moment-arm definitions and their paired checks:**

1. **Patellar-tendon arm, two ways.**
   - `r_geo` is the perpendicular distance from the numerically computed instant centre to the tendon line.
   - `r_vw` is the virtual-work value, +∂L_PT/∂θ with P held fixed.
2. **Effective quadriceps arm, two ways.**
   - `r_eff_vw` is +dL_Q/dθ over the whole chain.
   - `r_eff_eq` is `r_geo`·F_PT/F_Q, where the force ratio comes from tangential equilibrium of a frictionless point patella. This is the more independent of the two checks.

**Command** (deterministic; results.json sha256 ef8b7300…3f5 on two runs):

```
python3 patella_range.py          # writes results.json; Python 3.11, no dependencies
```

## Baseline results (synthetic)

| θ (°) | 0 | 15 | 30 | 45 | 60 | 75 | 90 | 105 | 120 |
|---|---|---|---|---|---|---|---|---|---|
| r_geo = r_vw (mm) | 32.32 | 32.30 | 31.69 | 30.62 | 29.23 | 27.65 | 26.02 | 24.48 | 23.16 |
| r_eff (mm) | 34.61 | 36.14 | 33.75 | 29.07 | 23.27 | 17.19 | 11.37 | 6.12 | 1.54 |

- **Tendon arm peak:** 32.40 mm at 6.93°, in both methods (difference 2e-11 relative). It is only 0.07 mm above the 0° value, so the interior peak is weak.
- **Effective arm peak:** 36.21 mm at 12.24°. The equilibrium and virtual-work methods differ by 1.5e-10 relative. The worst pointwise difference is 1e-8.
- **Instant-centre check:** the two-point instant-centre estimates agree to 1.2e-8 mm.

## Sensitivity: interior vs endpoint peaks (r_geo, 2° grid)

- **Rollback s controls the regime.**
  - s = 0 and s = 0.25: the peak is at the endpoint (0°) for every tuberosity offset Tx from 15 to 45 mm.
  - s = 0.5: the peaks are interior but marginal, at 2.7–7.0° and only 0.01–0.07 mm above r(0).
  - s = 0.75: robust interior peaks at 20–55° (margin 0.8–1.8 mm).
  - s = 1.0: robust interior peaks at 40–101° (margin up to 4.9 mm).
- **Boundary where dr/dθ at 0° is zero (bisection):** s* = 0.473, 0.423, 0.408 and 0.423 for Tx = 15, 25, 35 and 45 mm. Below s* the peak is at the 0° endpoint.
- **Tuberosity offset:** at fixed s ≥ 0.75, a larger Tx moves the peak toward extension. Interior peaks near 45° occur for (s 0.75, Tx 20 mm → 47.7°) and (s 1.0, Tx 35 mm → 47.3°). These synthetic parameter sets happen to land near the 45° figure from the search summary. **That is not evidence that the model matches cadaver data**, and the magnitudes (25–40 mm) sit below the 4–6 cm range from the same summary.
- **Patella start angle γ0:** 10° and 20° give an endpoint peak; 35°, 50° and 65° give interior peaks at 6.9°, 18.5° and 28.4°.
- **Across the sweep:** in all 34 feasible cases, the peak type and angle from r_geo and r_vw agree, with worst relative difference 1.5e-9.

## Uncertainty, failure cases and counterexamples

- **Weak independence.** The r_geo vs r_vw agreement at about 1e-10 is close to an identity for a planar rigid body: both come from the same kinematic map, and only the finite-difference and instant-centre numerics differ. The equilibrium vs virtual-work check on the effective arm is more independent, but it still shares the geometry. Passing either check shows the code is self-consistent. It says nothing about whether the model is anatomically correct.
- **Numerical precision.** The peak angle is quantized by a 0.5° grid (baseline) or 2° grid (sweep), with parabolic refinement. The finite-difference step is 1e-4°.
- **Counterexample: degenerate hinge.** With s = 0 and a trochlea concentric with the condyle, r is exactly constant, so "peak location" is undefined. An earlier exploratory version showed this, and it is the reason the trochlea centre is offset. This was an exploratory model correction made before the final run; the criterion was not changed.
- **Counterexample: fragile classification.** At s = 0.5 the interior peaks have margins under 0.1 mm. Any real measurement noise would make interior vs endpoint undecidable there.
- **Failure: infeasible geometry.** s = 1.0 with Tx = 45 mm has no solution: the patella constraint cannot be satisfied at θ = 118°.
- **Missing case:** no endpoint peak at 120° appeared inside the swept bounds (s = 1, Tx = 15 comes closest, with its peak at 101°).
- **Planar vs 3-D axis.** Krevolin et al. used a 3-D finite screw axis; this model uses a planar instant centre. Both the definitions and the reference frames differ.

## Limitations

- The model has a single planar degree of freedom: a circular condyle and trochlea, a point patella, no patellar tilt or flexion law, no soft tissue, and no 3-D effects.
- The parameters are illustrative and were not fitted to any subject.
- Primary sources were inaccessible, so no empirical claim is made. The empirical criterion stays **UNKNOWN**.
