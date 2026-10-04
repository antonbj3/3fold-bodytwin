# MECHANISM MSK — elastic-band posterior-chain loading demo (2026-07-21)

Executes `docs/MECHANISM_MSK_BUILD_PLAN.md` §5 step 2, elaborated per the operator's explicit ask ("resistance band
to load, posterior chain ground floor and the force continues through it" -- a resistance band anchored to
the floor, loading the posterior chain, with the force propagating up the kinetic chain). Every number below
is machine-measured this session; where a number is an assumption (not measured/sourced), it is flagged
explicitly. Script: `scripts/msk/attach_band.py` + `scripts/msk/band_config.json` (the hotswappable dial).
Run: `source_repository/.venv-msk/bin/python3 scripts/msk/attach_band.py`.

## Status: band attached YES. Numbers below are VERIFIED, not the tool's raw output --
## a real OpenSim InverseDynamicsTool bug was found and worked around (see §3).

---

## 1. Band attachment (yes, and how)

One `OpenSim::PathSpring` ("band"), added to the scaled `LaiArnoldModified2017_poly_withArms_weldHand`
model (`/media/anton/8838D60F38D5FBDE/mechanism_data/LabValidation_withVideos/subject2/OpenSimData/Mocap/Model/
LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim`, the same Step-1-verified model):

- **Ground end:** a path point on the model's `ground` body, positioned at the (X,Z) midpoint of the real
  `r_calc`/`L_calc` heel markers (i.e. genuinely "under the feet", not guessed) at floor height, plus a
  0.10 m anterior offset (config: `anchor.ground_anchor_anterior_offset_m`) so the straight-line Tier-1 path
  does not graze the shin.
- **Body end:** a path point on the `torso` body, at the local-frame midpoint of the model's own
  `R_Sternum`/`L_Sternum` markers -- a real anatomical landmark already in the model, not an invented offset.
  Represents a band held in both hands close to the chest during a band-resisted hip-hinge / RDL (the task's
  own "hands/torso" option).
- **Routing:** Tier-1 straight-line 2-point path, no `WrapCylinder`/`WrapSphere`. Hotswappable per the build
  charter: add a wrap object later for a curved/limb-following route; zero new dependency, same `PathSpring`.
- Verified live: `PathSpring` is a `Force` but **not** an `Actuator`/`ScalarActuator`
  (`issubclass(opensim.PathSpring, opensim.Actuator) == False`) -- confirms it behaves as a passive element
  automatically included in the system dynamics (see §2), not something the ID tool's actuator-exclusion step
  touches.

## 2. Real band data (verified LIVE, not invented)

**Product:** Thera-Band Blue ("Heavy" tier, one of 8 color-coded resistance levels).

**Primary source, fetched and read directly (not a search snippet):** Uchida MC, Nishida MM, Sampaio RAC,
Moritani T, Arai H. "Thera-band elastic band tension: reference values for physical activity." *J Phys Ther
Sci.* 2016;28(4):1266-1271. Table 1 gives measured tension (kgf, mean of 10 trials) at 10 elongation levels
(25%-250%) for a 0.30 m resting-length test sample. Transcribed directly from the PDF and independently
re-verified via a second read-through:

| Elongation | 25% | 50% | 75% | 100% | 125% | 150% | 175% | 200% | 225% | 250% |
|---|---|---|---|---|---|---|---|---|---|---|
| Force (N), converted via ×9.80665 | 11.28 | 16.67 | 20.99 | 24.62 | 27.95 | 30.30 | 32.56 | 34.52 | 37.76 | 39.03 |

**Machine cross-check of the transcription:** re-fit the full 10-point curve (`numpy.polyfit`, kgf vs
dimensionless elongation fraction Δ) gives slope=1.2143, intercept=1.1413 -- matching the paper's own
independently-published Table 6 regression (`y = 1.14 + 1.216·Δ`, r²=0.958) to <0.3%. This confirms the
transcription is correct, not by eyeballing but by reproducing the source's own published statistic from the
raw table.

**Second, independent source (cross-check, not blindly substituted):** Performance Health "Thera-Band Color
Progression" chart (doc P05358, © 2011): Blue = 5.8 lbf @100% (25.80 N), 8.6 lbf @200% (38.25 N) -- **+4.8%
and +10.8% above the measured values respectively**, same direction and order of magnitude as the primary
paper's own headline finding ("manufacturer reference values are overestimates"). A third source embedded in
the primary paper itself (Page 2003, their Table 2) gives Blue=3.22 kgf@100% (31.58 N, **+28%** vs measured)
-- noted in `band_config.json` as the "adversary" that was deliberately NOT used (it would silently inflate
the reported tension by more than a quarter versus the directly-measured curve).

**Calibration method (Tier-1 is a zero-intercept linear law; the real curve is concave with a non-zero
effective intercept, so a single choice must be made and disclosed):** stiffness is set via a **secant through
the origin, exact at 100% elongation** -- i.e. `stiffness = F(100%) / stretch_at_100%`, so the PathSpring's
predicted tension is **numerically exact** at the demo's representative pose, at the cost of a known,
quantified deviation elsewhere on the curve (reported below, not hidden).

| Elongation | 25% | 50% | 75% | **100%** | 125% | 150% | 175% | 200% | 225% | 250% |
|---|---|---|---|---|---|---|---|---|---|---|
| Measured (N) | 11.28 | 16.67 | 20.99 | **24.62** | 27.95 | 30.30 | 32.56 | 34.52 | 37.76 | 39.03 |
| Tier-1 linear (N) | 6.15 | 12.31 | 18.46 | **24.62** | 30.77 | 36.92 | 43.08 | 49.23 | 55.38 | 61.54 |
| Deviation | -45% | -26% | -12% | **0%** | +10% | +22% | +32% | +43% | +47% | +58% |

Honest read: Tier-1's zero-intercept line under-predicts at low stretch and over-predicts at high stretch --
exactly the expected signature of forcing a straight line through a concave curve. Tier 2
(`Blankevoort1991Ligament`, quadratic-toe-region + linear, already confirmed present in this venv) is the
documented upgrade path if the sweep beyond ~50-150% elongation matters for a future use case.

**Dissipation: 0.02 s/m, ASSUMED, NOT sourced** -- no published TheraBand damping/hysteresis coefficient was
found. Flagged explicitly rather than invented-and-hidden. Irrelevant to the results below regardless (the
demo pose is static, velocity=0, so the dissipation term is inactive).

**Installed geometry (computed from real model forward-kinematics, not assumed):** path length at the chosen
pose = 1.5848 m -> resting_length = 0.7924 m -> installed stiffness = **31.06 N/m** -> stretch = 0.7924 m
(100.0% elongation, exact by construction) -> **tension = 24.615 N** (matches the source table's 100% value
exactly, as designed).

**Zero-tension control** (re-run every script execution, not a one-off check): setting `resting_length` to 2×
the path length gives tension = `0.000000000` N -- confirms the `PathSpring` law (`tension=0` below rest
length, linear above) is wired correctly, straight from the class's own `getStretch()`/`getTension()`.

## 3. Mechanism verification -- why a ground-reaction-force pair is load-bearing, not decoration

Verified live (joint-tree dump): this model's `ground_pelvis` is a free 6-DOF root, and **`hip_r`, `hip_l`,
and `back` are all parented directly off `pelvis_offset`** -- i.e. the legs and the torso are **parallel
branches**, not serial. A path from a fixed ground point to the torso/sternum therefore has **zero direct
geometric coupling** to hip/knee/ankle:

```
computeMomentArm(hip_flexion_r) = 0.00000 m   (exact)
computeMomentArm(knee_angle_r)  = 0.00000 m   (exact)
computeMomentArm(ankle_angle_r) = 0.00000 m   (exact)
computeMomentArm(lumbar_extension) = -0.18603 m   (nonzero -- back IS on the ground->torso path)
```

This was cross-validated two ways: (a) directly via `PathSpring.computeMomentArm`, and (b) via a
with-band-vs-without-band `InverseDynamicsTool` run on a torso-only attachment with **no GRF at all** --
hip/knee/ankle deltas came back **exactly 0.0000**, lumbar's delta matched the moment-arm prediction
(`-T·momentArm`) to machine precision. **This means the ankle→knee→hip propagation this build demonstrates is
real but INDIRECT** -- it is the physically-correct whole-body-equilibrium effect (exactly how a standing
person's feet transmit a band's pull through the legs via the floor), closed here by a symmetric left/right
ground-reaction-force pair (`PrescribedForce`, force-only equilibrium: ΣF=0 over gravity + band + 2×GRF,
applied at the real `r_calc`/`L_calc` marker locations, split 50/50 under the sagittal-symmetry assumption).
Moment equilibrium about the pelvis is **not** independently enforced by construction; it is measured as a
residual (§5) and reported honestly, not hidden.

## 4. Bug found and fixed: OpenSim 4.6 `InverseDynamicsTool` corrupts `knee_angle_r/l` + `hip_flexion_r/l`

**This is the main technical finding of this build step.** Initial results showed `knee_angle_r` moments of
**-1065 to -1288 N·m** and `hip_flexion_r` moments of **+1371 to +5741 N·m** -- physically absurd (a real
max-effort knee-extension moment tops out around 300 N·m; hip similar). The coordinator flagged this as almost
certainly a bug, not a real result, and directed a forced re-diagnosis rather than reporting it. Diagnosis
(each step a separate, falsifiable, machine-checked test):

1. **Ruled out spline/differentiation artifact:** the corrupted values are **bit-for-bit identical** across
   coordinate-file `dt` in [0.001, 0.1] s and row count in [10, 50] -- a numerical differentiation artifact
   would change with sampling; this doesn't, so it is a deterministic (wrong) answer, not noise.
2. **Independent first-principles cross-check:** generalized force at static equilibrium equals `dPE/dq`
   (finite difference on total system potential energy, computed via pure forward-kinematics + mass
   properties -- zero ID-tool machinery). Result for a neutral standing pose, gravity only:

   | coordinate | verified dPE/dq | ID reported | agreement |
   |---|---|---|---|
   | ankle_angle_r | 1.22 N·m | 1.17 N·m | matches |
   | knee_angle_r | **-0.69 N·m** | **-1252.62 N·m** | **1815× too large** |
   | hip_flexion_r | **1.31 N·m** | **1337.12 N·m** | **1021× too large** |
   | lumbar_extension | -7.59 N·m | -7.61 N·m | matches |

   **⚠ PERSISTENCE NOTE (2026-07-21):** this specific sub-table had no persisted artifact anywhere on
   disk (only the separately-reported hip-hinge-pose propagation table, Sec.5 below, was
   reproducible from a raw `.sto` file) — flagged in `docs/MECHANISM_TRUST_LEDGER.md` §10 item 7. Now
   fixed at source: re-run and persisted at
   `data/msk_smoketest/elastic_band/neutral_pose_idbug_results.json`
   (`scripts/msk/elastic_band_neutral_pose_idbug_check.py`, reuses this script's own `MODEL_PATH`/
   `generalized_force`/`total_effective_pe`/`read_sto_row`/`find_col` unchanged, at the model's
   real file-default pose — exactly this section's own stated falsifier). Fresh reproduction matches
   the table above almost to the reported digit: knee_angle_r verified=-0.6858, ID=**-1252.616**
   (1826x); hip_flexion_r verified=1.3133, ID=**1337.124** (1018x); ankle_angle_r verified=1.2229,
   ID=1.17025 (matches); lumbar_extension verified=-7.5911, ID=-7.61169 (matches) — and, as a bonus
   over-determination the original table didn't report, the LEFT side reproduces the same corruption
   independently: knee_angle_l verified=-0.8024 vs ID=-1238.532 (1544x), hip_flexion_l verified=1.2973
   vs ID=1340.171 (1033x). Full disposition: `docs/MECHANISM_TRUST_LEDGER_REMEDIATION.md`.

   Ankle and lumbar match the independent calculation to <5%; knee and hip are wrong by three orders of
   magnitude. All 6 `ground_pelvis` coordinates and `hip_adduction_r`/`hip_rotation_r` were separately checked
   and also match (all near-zero and correct) -- **the corruption is narrowly confined to exactly
   `knee_angle_r/l` and `hip_flexion_r/l`.**
3. **Ruled out the patellofemoral constraint:** disabling both `CoordinateCouplerConstraint`s
   (`patellofemoral_knee_angle_{r,l}_con`) changes nothing -- identical corrupted output.
4. **Structural correlate (most likely root cause, not fully proven at the C++ level):** `walker_knee_r/l` is
   the only joint in this model whose `CustomJoint` `SpatialTransform` couples **all six axes** (3 rotation +
   3 translation) to a **single** coordinate (`knee_angle_r`) via nonlinear `PolynomialFunction`/
   `MultiplierFunction`s -- the standard "rolling knee" model (verified by reading the raw XML). `hip_r`, by
   contrast, has three independent coordinates each driving exactly one axis (like the unaffected `back`
   joint) -- yet `hip_flexion_r` is ALSO corrupted, while `hip_adduction_r`/`hip_rotation_r` are not. Both
   corrupted coordinates are the "rotation1"/X-axis slot in their respective joint -- consistent with the
   knee's corrupted reaction propagating proximally along that specific axis through the standard recursive
   inverse-dynamics computation, while the orthogonal hip axes and all distal/other-branch coordinates
   (ankle, subtalar, mtp, lumbar, the contralateral analysis is symmetric) are untouched.
5. **Not a modeling error on our part:** this exact DOF (`knee_angle_r`) was independently IK-validated to
   0.049-0.271° median RMSE in the Step-1 smoke test -- the joint's forward kinematics are correct; only the
   ID Tool's generalized-force attribution for it is broken.

**Workaround (used for every number reported here):** every joint moment in this document is computed via a
from-scratch virtual-work method (`generalized_force()` in `attach_band.py`) -- finite difference on an
"effective potential" = gravity PE + the band's elastic PE (`0.5·k·stretch²`) + the GRF's virtual-work
potential (`-F·x`, exactly analogous to gravity's `-F·x = mgy`), evaluated directly via forward kinematics.
This method was cross-validated against the (uncorrupted) ID columns for ankle, lumbar, and all 6 pelvis
coordinates (agreement <0.5-5%) before being trusted for knee/hip, where no ID cross-check is available. The
`InverseDynamicsTool` is still run and its `.osim`/`.sto` artifacts are kept on disk (useful for the pelvis
residual QC, §5, and as a smoke-test of model validity), but its `knee_angle_r/l` and `hip_flexion_r/l`
columns are never read or reported. **A concrete, self-contained falsifier for anyone who wants to check
this**: run ID on this exact model at its own file-default pose (all coordinates 0) with gravity only, no
added forces at all -- `knee_angle_r_moment` and `hip_flexion_r_moment` come back at -1252.6 and +1337.1 N·m,
which cannot be correct for a person's leg hanging straight down under only its own ~13 kg weight.

## 5. Force propagation table (VERIFIED numbers -- ankle → knee → hip → lumbar)

Representative pose: hip-hinge/RDL (pelvis_tilt 25°, hip_flexion 45°, knee_angle 20°, ankle_angle 5°,
lumbar_extension 5° -- a chosen representative posture, not a measured trial; all values confirmed within
this model's coordinate ranges, no knife-edge/hyperextended values).

**Band tension at this pose: 24.615 N** (Blue TheraBand, ~0.79 m installed length, 100% elongation).

| coordinate | no band (N·m) | with band (N·m) | Δ (N·m) |
|---|---:|---:|---:|
| ankle_angle_r | 10.37 | 11.19 | **+0.82** |
| knee_angle_r | 121.16 | 120.44 | **-0.71** |
| hip_flexion_r | -254.65 | -256.30 | **-1.65** |
| lumbar_extension | -50.11 | -45.54 | **+4.58** |
| ankle_angle_l | 7.41 | 8.26 | +0.85 |
| knee_angle_l | 122.12 | 121.49 | -0.63 |
| hip_flexion_l | -255.80 | -257.53 | -1.73 |

**Reading the table:** the band's direct effect is largest at `lumbar_extension` (+4.58 N·m, ~9% of the
baseline) because the torso IS on the band's direct path (§3); the effect at hip/knee/ankle (0.6-1.7 N·m, ~1%
of baseline) is real but small and INDIRECT -- it is entirely mediated by the extra ground-reaction force
needed to also counteract the band (with-band GRF = 391.8 N vertical + 9.0 N horizontal per foot vs. 383.4 N
purely vertical with no band), not a direct moment arm. This is the expected, physically sensible signature of
"the load travels up through the floor/legs," at a MODEST band tension (24.6 N is light relative to a ~767 N
body weight) -- a heavier band tier (Black/Silver/Gold, all in the same verified table) would proportionally
scale the propagation; this is exactly the hotswappable dial the config file provides.

Baseline magnitudes were sanity-checked against normalized literature ranges for hip-hinge/RDL biomechanics
(order 1.5-3.5 N·m/kg bodyweight for hip extensor moment during a bodyweight hip-hinge is a commonly-cited
ballpark in strength-biomechanics literature): 254.65 N·m / 78.2 kg = 3.26 N·m/kg -- within, not outside, that
range. This is a soft plausibility check, not a tight external anchor (no specific published trial matching
this exact pose/subject was sourced); flagged as an honest gap, not oversold as equivalent to the knee
flagship's Grand-Challenge-anchored evidence class.

**Left/right asymmetry note:** R and L values differ by up to ~25% (e.g. ankle 10.37 vs 7.41 N·m) despite a
nominally symmetric pose. Traced to real anthropometric asymmetry in this subject's own marker set (e.g.
`r_calc` local position (-0.0243, 0.0225, -0.0198) m is not an exact mirror of `L_calc` (-0.0181, 0.0172,
0.0193) m) -- an honest feature of using a real scaled subject, not a bug.

## 6. Pelvis residual QC (honesty gate on the GRF simplification, not hidden)

The symmetric 50/50 GRF split satisfies total FORCE equilibrium exactly by construction, but was never
required to satisfy moment equilibrium about the pelvis independently -- that gap is measured, not assumed
away:

| coordinate | no band | with band | body-weight scale |
|---|---:|---:|---:|
| pelvis_tilt | -605.5 N·m | -603.8 N·m | 766.9 N |
| pelvis_list | -0.03 N·m | -0.03 N·m | 766.9 N |
| pelvis_rotation | 0.01 N·m | 0.01 N·m | 766.9 N |
| pelvis_tx/ty/tz | ~0.0 N | ~0.0 N | 766.9 N |

**Honest read: `pelvis_tilt`'s residual (-605 N·m) is LARGE** -- roughly the same order as the whole system's
characteristic moment scale (body-weight × a representative length, ~613 N·m). This was independently
re-verified via the same virtual-work method (-605.65 N·m, matching ID's -605.52 to <0.03% -- confirming this
large number is a genuine, correctly-computed consequence of the simplified GRF choice, **not** a recurrence
of the §4 bug). It means: applying the GRF at the actual heel-marker location with a plain 50/50 split gets
total force right but leaves a substantial rotational demand unaccounted for at the pelvis -- a real
limitation of this Tier-1 simplification (a true center-of-pressure-based or explicitly moment-balanced GRF
would shrink this). Flagged here rather than hidden; translational residuals (tx/ty/tz) and the two other
rotational residuals (list, rotation) are all near-zero, i.e. well-closed.

## 7. Honest gaps (full list)

1. **Tier-1 linear band law deviates from the measured curve by -45% to +58%** away from the 100%-elongation
   calibration point (§2 table) -- an inherent zero-intercept-vs-concave-curve mismatch, not a coding error.
   Tier 2 (`Blankevoort1991Ligament`) is the documented upgrade path.
2. **Dissipation (0.02 s/m) is assumed, not sourced** -- no published TheraBand damping coefficient found.
   Inactive in this static demo regardless.
3. **Pelvis-tilt residual is large** (§6) -- the symmetric-GRF-split simplification is force-balanced but not
   moment-balanced; a real gap, quantified and disclosed, not resolved here.
4. **Single representative pose, not a measured trial or a sweep** -- this is a Tier-1 static demonstration at
   one chosen hip-hinge posture, not a dynamic/EMG-validated movement.
5. **Baseline magnitude plausibility is a soft (normalized-literature-range) check, not a tight external
   anchor** -- unlike the knee flagship's Grand-Challenge in-vivo anchor, no specific published trial was
   matched to this exact pose/subject.
6. **The OpenSim ID Tool bug (§4) is diagnosed and worked around, not fixed upstream** -- if this model family
   (LaiArnold/Rajagopal-lineage "rolling knee") is used for ID in any FUTURE script, the same workaround (or an
   upstream fix / OpenSim version check) is needed; do not assume a fresh script reading raw ID knee/hip
   columns is safe.
7. **Static optimization / muscle activations were not computed** (the task's "and/or" permitted joint moments
   alone; SO was descoped after the ID-tool bug consumed the debugging budget for this pass) -- honest gap,
   not attempted-and-hidden. The verified virtual-work moments in §5 are a valid SO target (net required
   moment) for a future pass.

## 8. Exact next build step

1. **Static Optimization for muscle activations**, using the §5 verified moments as the SO target (bypassing
   the ID tool's own moment computation, since it's confirmed unreliable for knee/hip in this model) --
   gives the "muscle activations WITH vs WITHOUT band" half of the original ask that this pass didn't reach.
2. **Re-provision Pose2Sim** (per `docs/MECHANISM_MSK_BUILD_PLAN.md` §5 step 3) to extend this same
   force-propagation method to a video-derived trial once a resistance-band clip is available (§6 of the
   build plan: no band-specific clip exists yet in `ig_downloads`).
3. **If this OpenSim ID Tool bug matters beyond this script** (e.g. for the KNEE-CELL flagship's own tooling,
   if it ever runs stock `opensim` ID on a LaiArnold/Rajagopal-lineage model rather than the JAM C++ fork it
   currently uses): flag it upstream / check whether a newer `opensim-core` release fixes it, since this was
   diagnosed only empirically here, not traced to a specific C++ line.

## Files

- `scripts/msk/band_config.json` -- the hotswappable dial (band product/citation/table, target elongation,
  pose, anchor description, GRF assumption). Change any value here to re-parameterize.
- `scripts/msk/attach_band.py` -- the build + verification script (geometry, calibration, bug workaround,
  propagation table, JSON summary).
- `data/msk_smoketest/elastic_band/model_{no_band,with_band}.osim` -- the two built models (band
  present/absent, everything else identical).
- `data/msk_smoketest/elastic_band/coords_{no_band,with_band}.sto`, `id_{no_band,with_band}/` -- the ID Tool
  inputs/outputs (kept for the pelvis-residual QC and as a reproducibility artifact; knee/hip columns therein
  are the KNOWN-CORRUPTED ones, do not read them directly -- use `band_results_summary.json` instead).
- `data/msk_smoketest/elastic_band/band_results_summary.json` -- machine-readable results (geometry,
  calibration, moment arms, Tier-1 sweep check, the verified propagation table, pelvis residuals).
