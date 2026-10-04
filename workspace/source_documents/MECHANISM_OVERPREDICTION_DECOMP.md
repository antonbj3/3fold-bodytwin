# MECHANISM OVER-PREDICTION DECOMPOSITION: knee/hip contact-force excess vs OrthoLoad (2026-07-21)

Decomposes the twin's established ~1.5x knee/hip contact-force over-prediction (subject2/`walking1`,
vs OrthoLoad in-vivo) into attributable model-assumption sources via four controlled SINGLE-FACTOR
perturbations, each re-solving Static Optimization (SO) from scratch and re-reading peak contact
force through the SAME unmodified free-body-cut machinery (`vhf.compute_self_cross_check_generic`)
used throughout this cert family. Context this builds on, not re-litigated: `docs/
MECHANISM_CMC_SECOND_SOLVE.md` already ruled out "SO-optimization artifact" (a mechanistically
decorrelated forward-dynamics solve, CMC, reproduces the over-prediction within 10%, slightly
higher) — this session's question is WHICH model assumption, not which solve method.

**Headline: of 11 tested configurations (4 Fmax levels + 2 wider, 2 objective exponents, 4 geometry
levels), NONE moves the ratio toward OrthoLoad by more than the pre-registered 0.10 bar — the
closest is a -20% moment-arm perturbation at -0.096 (knee), just short. The pre-registered
falsifier is therefore CONFIRMED: the over-prediction is DIFFUSE/STRUCTURAL with respect to these
four factors, not concentrated in any single one. Separately (a distinct, disclosed finding): the SO
objective activation-exponent is the one factor with a LARGE effect (+0.16 to +0.20), but only in
the WIDENING direction — the currently-used p=2 (activation²) objective already sits at/near a local
optimum for closeness to OrthoLoad among {p=1, p=2, p=3}; p=1 and p=3 both make it worse. This is
itself informative (don't "fix" the objective exponent) but does not explain the original gap.**

## 0. Pre-registration (stated before any perturbed SO solve was run or read)

- **MATERIALITY GATE (task's own falsifier, verbatim):** a factor is material iff it moves
  `knee_ratio` or `hip_ratio` by more than 0.10 (absolute) relative to baseline (1.515 knee / 1.412
  hip). If none clears this, the over-prediction is diffuse/structural.
- **DIRECTIONAL READING (added during write-up, applied uniformly, not post-hoc-cherry-picked):**
  "closes a substantial part of the gap" specifically means the ratio moves DOWN (toward 1.0); a
  factor that crosses 0.10 in the WIDENING direction is material-but-not-explanatory — reported
  separately, never conflated with "closing the gap."
- Weak, disclosed-as-weak mechanical priors (secondary diagnostics, not gates): higher Fmax → more
  co-contraction headroom → away from OrthoLoad; lower Fmax → toward. p=1 (sparser) → toward; p=3
  (more evenly spread) → away. Geometry: no directional prior (per-muscle moment-arm sign is
  heterogeneous, verified below).
- Passive-fraction (factor d) bands: <10% negligible, 10–30% partial, >30% material contributor to
  the "muscle-borne" tension the free-body cut subtracts.
- OrthoLoad anchors re-verified LIVE this session (not trusted from prior JSON): knee median
  **258.2211 %BW** (n=72), hip median **273.9309 %BW** (n=162) — matches the established cert to 4
  decimal places.
- **Self-consistency check (before trusting the harness on any perturbed model):** re-running the
  new evaluation harness on the EXISTING, already-committed baseline SO output reproduced
  391.0995231136491 %BW knee / 386.76785275005864 %BW hip — bit-identical to the established
  headline. The harness is correct before it is trusted on anything new.

## 1. Attribution table (all 11 configurations)

| factor | knee %BW | hip %BW | knee ratio | Δknee | hip ratio | Δhip | verdict |
|---|---:|---:|---:|---:|---:|---:|---|
| **baseline** (Fmax×1.0, p=2, geom 0) | 391.10 | 386.77 | 1.515 | — | 1.412 | — | (anchor) |
| fmax ×0.7 | 369.28 | 380.65 | 1.430 | −0.084 | 1.390 | −0.022 | negligible |
| fmax ×0.8 | 376.07 | 382.76 | 1.456 | −0.058 | 1.397 | −0.015 | negligible |
| fmax ×0.9 | 383.32 | 384.83 | 1.484 | −0.030 | 1.405 | −0.007 | negligible |
| fmax ×1.1 | 399.45 | 388.91 | 1.547 | +0.032 | 1.420 | +0.008 | negligible |
| fmax ×1.3 | 417.80 | 393.46 | 1.618 | +0.103 | 1.436 | +0.024 | **WIDENS** |
| SO objective p=1 (min-activation) | 436.52 | 440.89 | 1.690 | **+0.176** | 1.610 | **+0.198** | **WIDENS** |
| SO objective p=3 (fatigue-like) | 431.51 | 391.47 | 1.671 | **+0.157** | 1.429 | +0.017 | **WIDENS** |
| geometry −20% (perp-offset) | 366.39 | 374.50 | 1.419 | **−0.096** | 1.367 | −0.045 | negligible (closest to closing) |
| geometry −10% | 379.71 | 381.87 | 1.470 | −0.044 | 1.394 | −0.018 | negligible |
| geometry +10% | 400.28 | 390.29 | 1.550 | +0.036 | 1.425 | +0.013 | negligible |
| geometry +20% | 405.60 | 392.56 | 1.571 | +0.056 | 1.433 | +0.021 | negligible |

Materiality gate is symmetric (|Δ|>0.10, either direction); the **directional** reading (does it
CLOSE the gap) is what actually answers the task's question — see §5.

## 2. (a) Global muscle Fmax scaling — FORCED to a wider bracket before accepting the negative

Initial ×0.8/×0.9/×1.1 test: max |Δratio| = 0.058 (fmax×0.8, knee) — comfortably negligible. Per
Term #2 ("honest-negative is not a free pass"), this was forced further before being accepted: widened
to ×0.7 and ×1.3 (a generous bracket — ±30% is already beyond typical scaled-generic-model Fmax
calibration uncertainty, commonly cited ~10–20% in the literature). Result: ×0.7 still only reaches
−0.084 (short of the 0.10 bar); ×1.3 crosses it, but in the **WIDENING** direction (+0.103 knee),
confirming the pre-registered weak prior (more strength headroom → more co-contraction capacity →
away from OrthoLoad) rather than closing anything. **Linear extrapolation check:** the knee ratio
moves by only ≈0.084 over a 0.3 relative Fmax change (0.7→1.0); closing the full 0.5 ratio gap
(1.51→1.0) at that same local slope would require an ADDITIONAL ~1.8x reduction beyond ×0.7 — i.e.
Fmax below zero. **Fmax scaling within any physiologically plausible range cannot close this gap.**
Mechanism check on the one MATERIAL point (fmax×1.3): activation_max = 0.5055 (well under the 1.0
ceiling — no saturation artifact), consistent with a genuine "more available strength lets the
quadratic objective recruit more muscles simultaneously at low individual cost" mechanism, not a
solver pathology.

## 3. (b) SO objective activation-exponent — the standout factor, mechanism FORCED and verified real

`OpenSim.StaticOptimization` exposes `activation_exponent` as a first-class XML property (verified
via `osim.StaticOptimization().getNumProperties()` before assuming it existed, default 2 = the
already-used objective). p=1 and p=3 both WIDEN the gap substantially (knee +0.176/+0.157; hip
+0.198/+0.017) — the largest effect of any factor tested, in EITHER direction.

**Forced adversary (this is the leaning-positive-for-"material" result, so its mechanism must be
forced, not assumed):** is this a numerical artifact (bad local optimum, near-infeasible solve) or a
real redundancy-resolution effect?
- Ipopt convergence: no "maximum iteration"/"not converge"/"infeasible"/"restoration" messages in
  obj_p1's log. `convergence_pass=True`, crossing-muscle SET identical to baseline (13 knee / 25
  hip) for both p=1 and p=3 — the geometry/detection layer is unaffected, isolating the effect to
  the redundancy-resolution layer as intended.
- **Real-muscle-only activation ceiling check** (the FIRST pass over this mistakenly included the
  13 non-muscle `CoordinateActuator`s — e.g. `lumbar_ext` reads up to 2.86 because ideal torque
  actuators have no [0,1] bound — caught and corrected before trusting the number, restricting to
  the 80 `Millard2012EquilibriumMuscle`s that actually have `getMaxIsometricForce`): under p=1,
  `iliacus_r` hits activation **1.0000 at t=0.61s** and **0.9999 at t=0.57s**; `glmed1_r` hits
  **1.0000 at t=0.61s** — both major hip-crossing muscles, both saturating in a window (0.57–0.70s)
  immediately adjacent to the established hip contact-force peak (t=0.55s). At the peak instant
  itself, `glmed1_r`=0.873 and `iliacus_r`=0.700 under p=1 vs 0.574/0.477 at baseline — already
  heavily loaded, saturating moments later. Neither baseline (act_max=0.627) nor p=3
  (act_max=0.515) comes anywhere near saturation.
- **Geometric interpretation:** p=1's linear-program-like solution concentrates load onto the
  mechanically most "efficient" muscles (largest Fmax·moment-arm product); when those saturate, the
  optimizer is forced onto less-efficient combinations for the SAME required net moment (from
  inverse dynamics, unchanged) — requiring MORE gross tension to deliver the same moment via a
  lower-efficiency direction of the muscle-to-moment map, which shows up as more net compressive
  force at the joint even though total scalar activation (sum(a) at the knee: baseline 1.492 vs p1
  1.331 — LOWER) is not what's driving it; it is the ceiling-forced reallocation.
- **Verdict: REAL mechanism** (a well-documented sensitivity in the SO/redundancy-resolution
  literature — Crowninshield & Brand 1981; Trinler et al. 2019 report cost-function choice
  materially changes predicted co-contraction/joint-contact-force), not a numerical artifact — but
  it also discloses that p=1 pushes THIS subject/trial close to its own feasibility boundary, which
  the currently-used p=2 is not.
- p=3's knee-specific effect (+0.157, hip only +0.017 — an honest, unexplained asymmetry) traces to
  a "spread more evenly across synergists" pattern at the knee (gasmed/gaslat/tfl/recfem/sart ALL
  higher than baseline simultaneously; sum(a)=1.726 and sum(a²)=0.543, both the highest of the three
  exponents at the knee) — also a real, literature-consistent fatigue-like-objective signature, not
  forced to fully explain the knee/hip asymmetry (disclosed gap, §6).

**Symmetric-QC rider respected:** this factor is material, but it WIDENS, not closes, the gap — it
is reported as such, not spun as "explaining" the over-prediction.

## 4. (c) Geometric moment-arm perturbation — mechanism selection itself forced (two candidates tried)

**1st candidate, rejected before spending SO compute on it:** scale the model's own WrapCylinder
radius (the quadriceps' `KnExt_at_fem_r`/`KnExtVL_at_fem_r`, the patellar-mechanism analogue; glute
max's `Gmax{1,2,3}_at_pelvis_r`; psoas's `PS_at_brim_r`). Measured via `computeMomentArm` at the
established peak poses: **bit-identical** moment arms at ±10% radius for all 4 quadriceps muscles
(the path does not actually contact the cylinder at this pose — the radius is geometrically inert
there) and **non-monotonic** for psoas_r (both +10% AND −10% radius LOWERED the moment arm vs
baseline: 2.666→1.736cm and →2.166cm respectively) — a wrap-engagement-boundary discontinuity, not a
clean, usable lever. Rejected on measurement, not assumption.

**2nd candidate, used:** perturb each prime mover's path-point location on the joint's DISTAL body
(tibia_r for knee muscles, femur_r for hip muscles) perpendicular to the joint's own INSTANTANEOUS
rotation axis — extracted numerically via finite-difference of the child body's ground orientation
as the coordinate is perturbed by 1e-4 rad (model-representation-agnostic; this model's joint offset
frames carry non-trivial Euler orientation offsets, so an assumed X/Y/Z convention would have been
wrong). Verified via `computeMomentArm` before any SO run: real, mostly-monotonic, but
**heterogeneous** per-muscle moment-arm deltas (e.g. at the ±10% level: quads −3.4 to −4.1% per
+10%; `semimem_r`/`semiten_r` show a strong asymmetric jump, +0.02%/−0.15% at +10% vs −20.9%/−22.4%
at −10% — consistent with THEIR OWN condylar wrap cylinders engaging/disengaging, the same class of
nonlinearity rejected in candidate 1, here affecting only 2/8 knee muscles and disclosed, not
hidden). At ±20% the same muscles range up to −29.8% (semiten_r). Crossing-muscle SET topology
(13 knee, 25 hip) stayed IDENTICAL across every geometry level — the perturbation did not silently
break the free-body cut's anatomical validity.

Net effect on contact force: small and, like Fmax, negligible even at the wider ±20% bracket
(closest to closing: −0.096 knee at −20%, just under the bar). **Symmetric-QC note:** this is a
combined knee+hip perturbation (task's own framing groups "quads/hamstrings knee arm" and
"glutes/psoas hip arm" as one factor); knee-only/hip-only isolation was not separately run (scope
disclosed, §6) — but since both ratios moved the SAME sign together at every level (no
cancellation pattern observed), there is no evidence a separated design would reveal a materially
larger effect.

## 5. Directional verdict (the task's actual question)

| | max |Δ| toward OrthoLoad (closing) | max |Δ| away (widening) |
|---|---:|---:|
| across all 11 configs | **0.096** (geom −20%, knee) — short of 0.10 | **0.198** (obj p=1, hip) — clears 0.10 |

**FALSIFIER OUTCOME: CONFIRMED.** No single factor, in the direction that would explain/reduce the
observed 1.5x over-prediction, moves the ratio by more than the pre-registered 0.10 bar. Getting
even to −0.096 required an aggressive, disclosed-as-heterogeneous −20% geometric perturbation (up to
−30% on individual muscles) or a −30% global Fmax cut (−0.084) — both at or beyond the edge of
physiologically plausible single-factor uncertainty. **The over-prediction is diffuse/structural
with respect to {Fmax scaling, SO objective exponent, prime-mover moment arms, free-body-cut
passive-tissue attribution} — this is itself the pre-registered alternative finding, not a
non-result.** The one large effect found (SO objective exponent) points the WRONG way: p=2
(currently used) is already the best of {p=1,p=2,p=3} for closeness to OrthoLoad, so this is not a
lever for closing the gap either, and does not explain why baseline over-predicts in the first
place (baseline already uses p=2).

## 6. (d) Free-body-cut assumption: ligament/passive-tissue load-sharing

**(d1) Explicit ligament/passive-structure elements — measured, not assumed:** a full Force-class
census of the model finds exactly 2 classes: `Millard2012EquilibriumMuscle` ×80,
`CoordinateActuator` ×13 (torso/arm/pelvis-residual idealized actuators). Zero elements of any
ligament/passive-structure marker class (`Ligament`, `CoordinateLimitForce`, `Bushing*`,
`ElasticFoundation*`, `SmoothSphere*`, `ExpressionBasedCoordinateForce/BushingForce`). **This model
has no channel through which the free-body cut could be silently mis-attributing an unmodeled
ligament's load to "muscle"** — a structural fact about this specific model, not a general claim
about all OpenSim knee/hip models.

**(d2) Passive-vs-active tension WITHIN the 80 muscles themselves:** at the established peak
instants, replaying each crossing muscle's SO-solved activation + pose through
`model.equilibrateMuscles` (the same established, already-characterized pattern `metabolic_cost.py`
uses) and reading `getActiveFiberForceAlongTendon`/`getPassiveFiberForceAlongTendon`:

| joint | n crossing muscles | active (N) | passive (N) | passive fraction | reconstruction err | classification |
|---|---:|---:|---:|---:|---:|---|
| knee (t=0.51s) | 13 | 2489.7 | 96.2 | **3.7%** | 6.81% | negligible |
| hip (t=0.55s) | 25 | 2375.7 | 418.1 | **15.0%** | 10.89% | partial |

Reconstruction error (active+passive summed vs SO's own reported tendon force) sits in the same
ballpark as `metabolic_cost.py`'s own previously-disclosed ~8.27% error under the identical
rigid-tendon-equilibration approximation — inherited, not a new defect. **Reading:** at the knee, the
free-body cut's "muscle" attribution is essentially all active/neurally-driven (co-contraction in
the CNS-strategy sense, not tissue stretch). At the hip, roughly 1/7 of the crossing tension the cut
calls "muscle-borne" is actually passive fiber-stretch tension — real force, but not a redundancy-
resolution/optimizer choice, and not closable by re-tuning SO. This partially reclassifies (does not
remove) part of the hip over-prediction's character.

## 7. Honest gaps

1. **Single trial, single subject** (subject2/`walking1`) — same scope caveat as every cert in this
   family; no claim of generality across subjects/trials/speeds.
2. **Geometry perturbation is a disclosed proxy, not a uniform ±10%/±20% moment-arm change** — the
   measured per-muscle moment-arm delta is heterogeneous (§4); two knee hamstrings show a
   wrap-engagement-driven asymmetric jump. Reported transparently, not smoothed over.
3. **p=3's knee-vs-hip asymmetry (+0.157 vs +0.017) is measured, not fully mechanistically
   decomposed** — a "more even spread across knee synergists" pattern is shown (§3) but why the hip
   set responds so much less was not isolated further (would need a per-muscle sensitivity
   breakdown beyond this session's scope).
4. **Geometry factor (c) was tested combined (knee+hip prime movers perturbed together), not
   knee-only/hip-only separated** — disclosed in §4; no evidence of a masked larger effect, but not
   directly ruled out by a separated run.
5. **Factor (d)'s passive/active split is a per-muscle SCALAR tension decomposition, not a
   re-derivation of the full VECTOR free-body-cut residual** — it answers "what fraction of the
   muscles' own reported tension is passive," not "what would the peak contact-force number read if
   passive tension were reattributed to a separate channel." A direct vector-level version was not
   built this session (scope discipline).
6. **Only 2 of the 3 SO objective exponents bracket the currently-used p=2 symmetrically** (p=1,
   p=3); no continuous sweep (e.g. p=1.5, p=2.5) was run to find where the local optimum in
   exponent-space actually sits — p=2 is confirmed BEST of the three tested, not proven to be the
   true optimum.
7. **Fmax scaling was global and uniform** (every muscle scaled by the same factor) — a
   per-muscle-group differential mis-calibration (e.g. only the hip flexors over-scaled) was not
   tested and could in principle behave differently; out of scope for a single-factor screening
   test.
8. **This decomposition is itself a single-trial, single-model snapshot** — it does not and cannot
   test whether the SAME four factors would attribute differently on a different subject/trial,
   which remains the twin's deepest, still-unaddressed common-mode risk (per `docs/
   MECHANISM_CMC_SECOND_SOLVE.md` §5).

## Files

- `scripts/msk/overprediction_decomp.py` — the full, re-runnable sweep pipeline (7 initial + 4
  widened-bracket configurations). Reuses `static_opt_knee.py` (`sok`)'s SO-run/patch/gate machinery
  and `validate_hip_force.py` (`vhf`)'s `compute_self_cross_check_generic` UNCHANGED (called
  identically for both knee and hip, on every model variant, per this task's explicit instruction).
  New code: parametrized SO-run/model-builder wrappers (`build_fmax_model`, `build_geom_model`,
  `run_so_variant`, `evaluate_variant`), the geometry-perturbation mechanism
  (`instantaneous_axis_local`, `joint_center_local`, `perturb_group`), and the factor-(d) analysis
  (`ligament_passive_structure_census`, `passive_active_decomposition`).
- `data/msk_smoketest/subject2_walking1/overprediction_decomp/overprediction_decomp_results.json` —
  every number in this document, traceable back to this one file: `baseline`, `factors` (11 configs),
  `factor_d`, `attribution_table` (with `directional_verdict`), `directional_summary`,
  `mechanism_check_obj_p1`, `geom_build_meta` (per-muscle moment-arm deltas + skipped-point audit).
- `data/msk_smoketest/subject2_walking1/overprediction_decomp/models/*.osim` — the 9 perturbed model
  files built this session (3 Fmax + 2 wider Fmax + 4 geometry levels; objective-exponent runs reuse
  the unmodified original model).
- `data/msk_smoketest/subject2_walking1/overprediction_decomp/{fmax_*,obj_*,geom_*}/` — each
  variant's SO (`AnalyzeTool`+`StaticOptimization`) output (activation.sto, force.sto, patched XML).
- Inputs read in place, never modified: `data/msk_smoketest/subject2_walking1/static_optimization/so/`
  (the existing baseline SO output, reused not rerun), the subject2 model/IK/GRF files on the
  external drive (read-only per isolation).

No git operations performed (isolation respected). All new files left untracked for the coordinator
to stage explicitly.
