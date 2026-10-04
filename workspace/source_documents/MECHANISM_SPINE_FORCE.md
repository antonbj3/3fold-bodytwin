# MECHANISM SPINE FORCE — muscle-driven lumbar compressive-force validation vs OrthoLoad VBR (2026-07-21)

Replicates, for the SPINE, the VALIDATED knee force-cert method (`docs/MECHANISM_JOINT_FORCE_VALIDATION.md`
+ `docs/MECHANISM_STATIC_OPT.md`): Static Optimization + a Newton's-law/virtual-work method IMMUNE to the
OpenSim 4.6 `InverseDynamicsTool` bug. Every number below is machine-measured this session
(`scripts/msk/validate_spine_force.py`, exit 0, **bit-for-bit identical across two independent full runs**),
not recalled. Isolation respected: `.venv-msk` only, model/OrthoLoad data read in place, no git commit/push.

## Headline result

| | Tier-1 (pure kinematic, muscle-free) | Tier-2 (muscle-driven, self-computed) |
|---|---:|---:|
| **no_box** (bodyweight-only forward bend) | 41.19 %BW | **296.59 %BW** |
| **with_box** (10 kg two-handed stoop lift) | 52.78 %BW | **389.17 %BW** |
| **IN-VIVO OrthoLoad spine_vbr** (stoop-lift, 10 kg, n=13 trials, 4 subjects) | — | **170.94 %BW** (median) |

**Ratio (Tier-2 with_box / in-vivo median) = 2.28×; (Tier-2 no_box / in-vivo median) = 1.74×.** Both land in
the **pre-registered EXPECTED direction** (twin > VBR in-vivo) — the opposite direction from the knee cert,
and expected for a structural reason stated *before* computing (§6): a vertebral-body-**replacement**
implant typically shares axial load with posterior pedicle-screw/rod instrumentation, so it reads a
*fraction* of true total spinal compression, not the whole thing. Both raw numbers (2274–2985 N) also sit
in the same order of magnitude as classic (recalled, not re-verified this session — see honest gaps)
lifting-biomechanics compression figures for a comparable stooped 10 kg lift.

**A real, machine-confirmed sign bug was found and fixed this session** in the shared crossing-muscle
force-direction pattern this build reuses from `scripts/msk/static_opt_knee.py` — see §5. Fixing it flipped
the Tier-2 result from a non-physical **negative** (tensile) number to the positive (compressive) one above.
The same bug appears to still be live in `static_opt_knee.py` itself (unfixed there, out of this build's
scope) — flagged for the parent/a future editor, not swept under the rug.

## 1. Inputs

- Model: `data/msk_models/LaiArnoldModified2017_erector_spinae_subject2_scaled.osim` (168 muscles = 80
  native lower-limb + 88 erector-spinae; subject2, 78.2 kg; regenerate via
  `scripts/msk/add_erector_spinae.py` if the gitignored file is absent — it was present this session).
- No captured "lifting" mocap trial exists for subject2 (only `walking1-3`, `DJ1-3`/`DJAsym*`, `squats1`,
  `squatsAsym1`, `STS1`, `STSweakLegs1`, `static1` — checked live). Per this cert family's own established
  precedent (`docs/MECHANISM_MSK_ELASTIC_BAND.md`, `docs/MECHANISM_BAND_POSTERIOR_CHAIN.md`,
  `docs/MECHANISM_BAND_STATIC_OPT.md` all use a representative STATIC pose, not a dynamic trial, for spine
  work), this build does the same — disclosed as honest gap #1, not hidden.
- OrthoLoad in-vivo spine anchor: `data/external/orthoload/spine_vbr/database_api/akf/*.akf` (869 files,
  subjects WP1/WP2/WP4/WP5 — Rohlmann-group vertebral-body-replacement instrumented telemetry).

## 2. Why a NEW pose, not the existing hip-hinge cert's pose (forced via OODA, not copied)

`docs/MECHANISM_BAND_STATIC_OPT.md` already forced (3 independent sign checks + an 8-point robustness
sweep) a real finding: at the existing hip-hinge cert pose (`pelvis_tilt=25, hip_flexion=45,
lumbar_extension=+5` — a **neutral spine riding on a tilted pelvis**), gravity's own generalized force at
`lumbar_extension` is **negative** (a flexion-direction restraint) — the 88 erector-spinae muscles (88/88
positive moment arm, true extensors) are structurally incapable of supplying it; Static Optimization there
is genuinely Ipopt-infeasible, not a bug. This task explicitly asks for "a forward-flexed lifting/hinge
posture... where the erector spinae carry it" — i.e. genuine **lumbar flexion**, not a neutral spine on a
tilted pelvis. Re-verified live here, independently, before committing to any pose:

| lumbar_extension (deg) | gravity-only required (N·m) | gravity+box required (N·m) | torso world pitch (deg) |
|---:|---:|---:|---:|
| +5 | −16.25 | +38.20 | −15 |
| −5 | −0.17 | +46.14 | −5 |
| −15 | +15.91 | +52.69 | +5 |
| −25 | +31.51 | +57.63 | +15 |
| **−35 (LOCKED)** | **+46.16** | **+60.82** | **+25** |
| −45 | +59.40 | +62.16 | +35 |

(Fixed at `hip_flexion=45°, pelvis_tilt=10°, knee=15°, ankle=10°` throughout — only `lumbar_extension`
varies, isolating exactly the variable this task is about.) The sign flips between −5° and −15°; the
locked pose (−35°) clears the **pre-registered ≥15 N·m margin** in both the bodyweight-only and
+box conditions (46.16, 60.82 N·m) — not a knife-edge. Both numbers are comfortably inside the erector
group's own force-length-adjusted ceiling at a bent pose (148.6 N·m, `docs/MECHANISM_ERECTOR_SPINAE.md`).

## 3. Box load: a disclosed, deliberate simplification (forced via OODA)

A first attempt routed the 10 kg box through the arm chain (`shoulder_flex`/`elbow_flex`), applying it at
the wrist markers. Diagnosed as unworkable **before** trusting any number from it: with a relaxed
(arm_flex=0) arm, hand height *rose* as trunk flexion increased (backwards — a real person actively
reaches down as they bend over); with an arbitrary reach angle (arm_flex=−20°), the box ended up so far
posterior that adding it *reduced* the required lumbar moment (also backwards). Rather than hand-tune
unverified arm kinematics further, the box is instead attached as a **fixed torso-local offset point**
(anterior +0.35 m, inferior −0.95 m from the R/L-Sternum midpoint) — this sidesteps arm-reach tuning
entirely while giving *direct, disclosed* control over exactly the parameter that governs the lumbar
moment: horizontal load-to-spine distance, the same single parameter the NIOSH lifting equation and
Potvin/McGill-style spine models hinge on. Verified live at the locked pose: the box lands **0.216 m above
the floor, 0.031 m anterior of the ankle** — a physically sane "box on the ground, about to lift" spot, not
asserted blindly.

## 4. Method (two tiers, mirroring the knee cert exactly)

BFS on the model's own joint tree (`scripts/msk/validate_joint_force.py`'s `get_descendant_bodies`, reused
not reimplemented) finds `back`'s descendant ("above-cut") set = `{torso, humerus_r/l, ulna_r/l,
radius_r/l, hand_r/l}` — 9 bodies, 35.54 kg = **45.4% of body mass** (a rough anthropometric sanity
context, not independently literature-checked this session). A live, topological (pose-independent)
path-point scan finds **exactly** 88/88 erector-spinae muscles cross this cut and **0/80** native muscles
do — matching `add_erector_spinae.py`'s own construction (all 80 native LaiArnold muscles are lower-limb
only) — re-verified here as a pre-registered external anchor (PASS), not assumed inherited.

- **Tier 1** (pure kinematic, immune to the ID bug by construction, never calls
  `InverseDynamicsTool`/`Solver`/`JointReaction`): static (a=0) Newton's law on the above-cut free body:
  `R = -(mass_above · g_vec) - F_box`. GRF is **irrelevant** to this specific cut (feet are on the
  pelvis-side, a parallel branch, per `docs/MECHANISM_MSK_ELASTIC_BAND.md`'s own kinematic-tree finding) —
  still built into the model for a physically complete, SO-solvable system, just not part of this formula.
- **Tier 2** (muscle-driven): Static Optimization (quasi-static 10-frame hold, same convention as
  `scripts/msk/band_static_opt.py`) gives each erector-spinae muscle's tension; every crossing muscle's
  live pulling-force vector is identified (never hardcoded) and summed, then **subtracted** from R to
  isolate `F_bone_contact = R − Σ F_muscle_crossing` — same "reaction is a net of contact+muscle" logic as
  `scripts/msk/static_opt_knee.py`'s self-computed method. Independently cross-checked by a
  moment-balance gate (required generalized force via virtual work vs. supplied generalized force from the
  SO solution's own actuator forces + live `computeMomentArm` calls) — **PASS**, near-exact agreement,
  both variants, all 3 lumbar coordinates.
- **Compressive component**: projection of `F_bone_contact` onto the torso body's own local +Y (cranial)
  axis, expressed in ground at the current pose — verified live that torso-local Y = world (0,1,0) at the
  neutral pose (a measured fact, not an assumed convention).

## 5. A real, machine-confirmed sign bug — found, isolated, and fixed (forced via OODA, not waved through)

A first run of this exact pipeline produced a **negative** Tier-2 compressive force (−214.2 %BW no_box,
−283.6 %BW with_box) — non-physical (a disc/facet joint can only push, never pull). Per the "one-shot fail
is not an honest negative" rule, this was forced to a root cause rather than reported as-is:

The crossing-muscle force-direction helper this build reuses the *pattern* of
(`scripts/msk/static_opt_knee.py`'s `knee_crossing_muscles_and_forces`) contains this index assignment:
```python
inside_idx, outside_idx = (k, k + 1) if in_k1 else (k + 1, k)
```
**Decisively shown wrong by a standalone, OpenSim-free unit test** (not trusted from symbolic reasoning
alone): a toy 2-point cable with an "inside" point 5 units cranial of an "outside" point, under tension,
must pull the inside point **caudally** (`normalize(outside − inside)` = `(0,−1,0)`, basic cable mechanics,
model-independent). The code's own index assignment computes `(0,+1,0)` for this exact case — **the exact
negative** of the correct answer, confirmed by `np.allclose(result, -ground_truth) == True`. The fix
(swap the two branches: `(k + 1, k) if in_k1 else (k, k + 1)`) is applied and commented in
`scripts/msk/validate_spine_force.py`'s `crossing_muscles_unit_dirs`. After the fix: Tier-2 compressive
force is positive in both variants (296.59 / 389.17 %BW), Tier-2 ≥ Tier-1 in both variants (co-contraction
adds compression, the expected direction), and the dose-response gate (with_box > no_box) passes.

**This is NOT fixed in `static_opt_knee.py` itself** — out of this build's explicit scope (spine, not a
knee-cert audit). Flagged here as a genuine, machine-verified, actionable finding for the parent / a future
editor: that script's "self_computed" knee contact-force headline (233.20 %BW / 391.11 %BW ratio-vs-JR
comparison in `docs/MECHANISM_STATIC_OPT.md`) uses the same buggy pattern for its own muscle-crossing
subtraction term. **Its official `opensim.JointReaction` cross-check pass is a completely separate OpenSim
code path and is unaffected** — but the self-computed number specifically has not been independently
re-verified against this fix. Whether/how much it would move is unknown without re-running that script with
the same fix — a bounded, well-defined follow-on task, not done here (scope discipline).

## 6. Why the expected direction here is OPPOSITE the knee cert's (pre-registered, not fitted after)

The knee cert expected (and found) the twin's kinematic/muscle-driven estimate to read **lower** than
in-vivo (co-contraction invisible to Newton's law alone). Here, the pre-registered expectation is the twin
reads **higher**: a vertebral-body-**replacement** implant is typically combined with posterior
pedicle-screw/rod instrumentation that **shares/shunts axial load away** from the instrumented body — and
OrthoLoad's own corpus-wide finding (`docs/MECHANISM_ORTHOLOAD_INDEX.md` §7) already flags this
qualitatively ("an instrumented vertebral segment does not carry anything close to full body weight...").
This build's pipeline estimates **total** compressive load at the lumbar joint (contact + all crossing
tissue), not the fraction a load-sharing implant sees. Measured ratio (1.74–2.28×) matches this direction —
the PASS condition here; a ratio < 1 would have been the surprising result needing explanation. Note: the
*exact* fraction a VBR implant captures was not independently, quantitatively pinned down this session
(WebSearch quota was exhausted mid-build — see honest gaps) — this is a qualitative, structurally-reasoned
expectation, not a precisely-calibrated correction factor.

## 7. Machine-checked gates (PASS/FAIL, not eyeballed)

| gate | pre-registered threshold | measured | verdict |
|---|---|---:|---|
| Model topology (external anchor) | 168 muscles, 88 erector, 9 above-cut bodies, 88/0 crossing | exact match | **PASS** |
| Sign/margin (locked pose, both loadings) | required moment > +15 N·m | +46.16 / +60.82 N·m | **PASS** |
| Symmetry (lumbar_bending/rotation, with_box) | \|required\| < 5 N·m | −0.134 / −0.063 N·m | **PASS** |
| Parallel-branch consistency, DIRECT channel (box→hip_flexion) | < 3 N·m | 0.000 N·m exactly | **PASS** |
| SO convergence (both variants) | 10/10 frames, 0 NaN, 0 bound violations, row-spread <1e-4 | 0 NaN, spread ~2e-8 | **PASS** |
| Moment-balance cross-check (required vs SO-supplied, 3 coords × 2 variants) | within max(2 N·m, 15%) | diff ≈0.000 throughout | **PASS** |
| Compressive-force sign (physical, not tensile) | Tier1 & Tier2 both > 0 | all 4 values > 0 | **PASS** |
| Dose-response (with_box > no_box) | strict inequality, both tiers | 52.78>41.19, 389.17>296.59 | **PASS** |
| Tier ordering (Tier2 ≥ Tier1, co-contraction adds compression) | Tier2 ≥ Tier1 | both variants | **PASS** |
| Reproducibility | bit-for-bit across 2 runs | identical | **PASS** |
| **Overall** | all of the above | | **PASS** |

Bonus, non-gated context: the box's effect on `hip_flexion` decomposes 100% into the **indirect**,
GRF-mediated channel (−33.5 N·m; GRF increases from 383→432 N/side to support the extra 10 kg) and exactly
**0%** direct channel — the identical direct/indirect geometric split `docs/MECHANISM_MSK_ELASTIC_BAND.md`
already established for this model family, independently re-confirmed here for a different load case, not
a new phenomenon.

## 8. OrthoLoad anchor — a precise activity+load match, not the pooled bucket (self-contained live re-parse)

`docs/MECHANISM_ORTHOLOAD_INDEX.md`'s own pooled spine "Lifting" bucket (210 trials, 37.3 %BW median) mixes
lying-position limb-lifts, lifting AIDS, and unspecified loads — not a clean match for this build's
scenario. This session's own live re-parse of the raw AKF files (reusing
`scripts/msk/validate_joint_force.py`'s proven `parse_akf`, not re-implemented) finds a literal match:
Comment #1 text `"VBR; Standing; lifting a weight from the ground, knees straight, back bent, weight =
10kg"` — a stoop-lift, matching this build's own pose choice (knees near-straight, spine flexed) **and**
load (10 kg) exactly, not just the activity label.

| condition | n trials | n subjects | median %BW | range |
|---|---:|---:|---:|---|
| **Stoop** (knees straight, back bent) 10 kg, two-handed | 13 | 4 (WP1,WP2,WP4,WP5) | **170.94** | 140.10–193.77 |
| Squat (knees bent, back straight) 10 kg, two-handed [bonus contrast] | 9 | 4 | 166.93 | 138.73–218.52 |

Excludes one-handed ("...with the right hand...") variants and the already-documented WP4
100 N-bodyweight-placeholder anomaly (`docs/MECHANISM_ORTHOLOAD_INDEX.md` §5.1, `bodyweight_N<300`
floor, same fix reused). A representative trial's own raw Fx/Fy/Fz components (checked live,
`wp5_280311_1_132.akf`, peak instant) show the resultant is **99.7% axial/cranial** — i.e. comparing
OrthoLoad's peak resultant force to this twin's *compressive* component specifically is not an
apples-to-oranges move.

## 9. Erector-spinae recruitment (Static Optimization solution, with_box, interior frame)

| subgroup | n muscles | n active (>0.05) | activation range |
|---|---:|---:|---|
| Iliocostalis lumborum | 22 | 22 | 0.182–0.500 |
| Longissimus thoracis pars lumborum | 10 | 8 | 0.015–0.260 |
| Longissimus thoracis pars thoracis | 28 | 28 | 0.103–0.514 |
| Multifidus | 28 | 26 | 0.041–0.433 |

No muscle saturates (peak ≈0.51 of max activation) — a moderate, non-degenerate recruitment level for a
10 kg stoop lift, consistent with the erectors' own 148.6 N·m force-length-adjusted ceiling at a bent pose
comfortably exceeding the 60.8 N·m required here. The pre-existing, still-present, **unbounded** 10 N·m
ideal `lumbar_ext` `CoordinateActuator` (matching this model family's own established convention — bounding
it makes SO infeasible, per `docs/MECHANISM_BAND_STATIC_OPT.md`) supplies ~19.2 N·m (control=1.92, i.e.
1.92× its "nominal" ±1 range) — roughly **31%** of the required 60.8 N·m moment, muscles supplying the
rest. This means the reported Tier-2 number (derived purely from the 88 real muscles' SO tensions) may be
a mild **under**-estimate of true muscle-driven compression, not an over-estimate, if the erectors alone
had to carry 100% of the moment.

## 10. Honest gaps (full list)

1. **Static pose, not a captured dynamic lifting trial** — no such trial exists for subject2 (§1); same
   simplification this whole spine-work sub-family already uses.
2. **Box load is a fixed torso-local offset, not routed through arm/wrist kinematics** (§3) — a deliberate,
   disclosed simplification; does not affect the Tier-1 Newton's-law force balance at all (location-
   independent), only the virtual-work lever-arm term, which is instead directly and transparently
   controlled.
3. **Single lumped lumbar hinge** — this model's own fidelity level (`back`: pelvis→torso, one joint, same
   as `docs/MECHANISM_ERECTOR_SPINAE.md`'s own construction). "Lumbar compressive force" here means the
   resultant at this one lumped joint, not a named vertebral level (L4/L5 vs L5/S1 not distinguished).
4. **VBR load-sharing is a qualitative, not quantitatively calibrated, correction** (§6) — the expected
   *direction* (twin > in-vivo) is pre-registered and confirmed; the exact fraction a VBR implant captures
   was not independently pinned down this session (WebSearch quota exhausted mid-build; a WebFetch did
   confirm the classic NIOSH-lifting-equation paper's DOI resolves live to a real tandfonline.com page, but
   its content itself was blocked, HTTP 403 — consistent with a known publisher bot-block pattern already
   logged elsewhere in this repo; the oft-cited NIOSH 3400 N/6400 N action-limit/MPL compression figures
   are offered only as **recalled**, not independently re-verified, secondary context, clearly weaker
   evidence than the live-reparsed OrthoLoad anchor).
5. **Static Optimization is an effort-minimizing solution, not measured EMG** — same well-documented
   structural limitation as the knee cert's own honest gap #4; this build's erector-spinae-only model
   (no anterior trunk flexors — rectus abdominis/obliques/psoas were explicitly out of scope in
   `add_erector_spinae.py`) cannot express intra-abdominal-pressure-mediated stabilizing co-contraction,
   which in real lifting studies usually *adds* net compression beyond an effort-minimizing estimate.
6. **The unbounded native `lumbar_ext` ideal actuator supplies ~31% of the required moment** (§9) — a
   disclosed, pre-existing model-family limitation, not new to this build.
7. **The crossing-muscle-force sign bug (§5) is fixed here, NOT in `static_opt_knee.py` itself** — a real,
   machine-verified, actionable, but unresolved-in-scope finding for a future editor.
8. **Different populations** — OrthoLoad's WP1/WP2/WP4/WP5 are VBR-implant patients vs subject2, a healthy
   young(ish) OpenCap participant (78.2 kg, 1.96 m); no per-subject anthropometric matching.
9. **Single representative pose per loading condition**, not a distribution across lift styles/speeds/
   subjects — same scope caveat as the rest of this cert family.
10. **Shear component reported (527–689 N, Tier-2) but not independently cross-validated** against a
    literature anchor — context only, not gated.
11. **`mass_above` (45.4% of body mass) is a soft, recalled (not freshly re-verified) anthropometric
    sanity context** — somewhat below commonly-cited generic trunk+head+arms fractions (~55–60%), plausibly
    reflecting this model's own measured/scaled segment masses rather than a generic table; not gated.
12. **GRF Tier-1 simplification** (symmetric 50/50 split, force-only equilibrium) — same disclosed
    simplification as the whole existing cert family.

## 11. Next step

1. **Re-audit `static_opt_knee.py`'s self-computed method against the §5 fix** to determine whether/how
   much the published 233.20 %BW / 391.11 %BW knee numbers shift — a bounded, well-defined follow-on, not
   done here (scope discipline: this build is the spine cert, not a knee-cert re-audit).
2. **If a captured lifting trial becomes available**: re-run with real kinematics + GRF + box-force
   time-histories (reusing `validate_joint_force.py`'s Savitzky-Golay differentiation pipeline) instead of
   the static-pose simplification.
3. **Route the box load through tuned arm/wrist kinematics** as a refinement if a fully anatomically
   literal hand-to-box picture is needed later (not required for the compressive-force conclusion itself).
4. **Extend the erector-spinae donor port to the anterior trunk flexors** (rectus abdominis, obliques,
   psoas — already flagged out-of-scope in `add_erector_spinae.py`) for an IAP-inclusive picture of net
   compression, closer to EMG-informed literature estimates.
5. **A dedicated literature dive** (WebSearch quota permitting) into the Rohlmann VBR device literature
   specifically, to quantitatively calibrate what fraction of true spinal load a VBR implant captures.

## Files

- `scripts/msk/validate_spine_force.py` — the full pipeline (self-contained, re-runnable; imports
  `attach_band.py`/`validate_joint_force.py`/`static_opt_knee.py`/`band_static_opt.py` for proven
  pose-setting/AKF-parsing/BFS/SO-patching/actuator-classification code, not re-implemented).
- `data/msk_smoketest/subject2_spine_stoop_lift/` — built `.osim` model variants (no_box/with_box, +
  reserve-merged), quasi-static `coords_*.sto`, per-variant SO output (`so_no_box/`, `so_with_box/`), and
  `spine_force_validation_results.json` (every number in this document, machine-written).
- `data/msk_models/LaiArnoldModified2017_erector_spinae_subject2_scaled.osim` — the model used
  (gitignored; regenerate via `scripts/msk/add_erector_spinae.py` if absent).
