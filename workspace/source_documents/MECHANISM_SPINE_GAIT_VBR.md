# MECHANISM SPINE GAIT VBR — does the lumbar spine over-predict in-vivo VBR load during WALKING like the knee/hip? (2026-07-21)

**Question.** Knee and hip both over-predict in-vivo OrthoLoad contact force by ~1.4–2.0×
(`docs/MECHANISM_JOINT_FORCE_SCORECARD.md`, `docs/MECHANISM_CROSS_SUBJECT.md`). The spine's own existing
cert (`docs/MECHANISM_SPINE_FORCE.md`) found 1.74–2.28×, but at a **static lifting pose** (no captured
lifting trial exists for subject2), on a **less complete** model (erector-spinae-only, no ligaments, no
trunk flexors), against OrthoLoad's **stoop-lift** VBR data — a different activity, a different model. This
build asks the task's real question directly: during **gait** (subject2 `walking1`, real IK+GRF), using the
**same Static-Optimization + Newton's-law crossing-tissue-subtraction architecture** the knee/hip certs used,
on the **unified/spine-ligament substrate** (trunk muscles + lumbar ligaments present, verified active — not
an inert spine), does the twin's lumbar compressive load over-predict OrthoLoad's VBR **walking** telemetry,
and by how much, once the VBR-fraction question is confronted head-on?

## Headline result

| | value | source |
|---|---:|---|
| **Model Tier-2 (muscle+ligament-subtracted) gait-cycle peak** | **210.22 %BW** @ t=0.09s | this build, machine-measured |
| Model Tier-1 (pure kinematic reaction, SO-independent) | 55.97 %BW @ t=0.68s | this build |
| **OrthoLoad spine_vbr level-walking** (this session's own live AKF re-parse) | **53.37 %BW** median (mean 56.71), n=42 trials, 4 subjects (WP1,2,4,5) | live re-parse |
| Published triangulation: Damm et al. 2017 (PMID 27914627) | 39.0 %BW mean | literature |
| **RAW ratio (model / own-median)** | **3.94×** | |
| RAW ratio (model / Damm 2017) | 5.39× | |
| Knee (subject2–4 range, `MECHANISM_CROSS_SUBJECT.md`) | 1.515×–1.990× | prior cert |
| Hip (subject2–4 range) | 1.412×–1.816× | prior cert |

**The spine's RAW over-prediction (3.9–5.4×) is larger than knee/hip's (1.4–2.0×) — directionally
consistent with the scorecard's own open question ("some fraction of the spine's gap could be the SAME
method-architecture over-prediction stacked on top of the load-sharing effect... nothing measured
distinguishes the two") — but a dedicated literature search this session could NOT find an externally
anchored VBR-device load-share fraction to decompose the two mechanisms quantitatively.** This is reported
as the honest, load-bearing limitation, not glossed over (§6). **Confidence tier: in-vivo-anchored (VBR),
load-share-fraction-limited.**

A **real bug was found and fixed** during this build (a locked constraint-dependent coordinate creating a
genuine kinematic contradiction, not a torque-margin problem) — full forced-diagnosis trail in §4, because a
one-shot "PASS" (no NaN, bounded activations) was NOT accepted as proof of a good result, per this repo's own
symmetric-QC discipline.

---

## 1. Pre-registration (stated before trusting any Tier-2 number)

- **Falsifier** (task's own framing): spine ALSO over-predicts ~1.4–2.0× after correct load-share accounting
  → systemic (same architecture-driven mechanism as knee/hip); spine matches VBR well → joint-specific.
- **Regime-sanity gate**: ratio ∈ (0.3, 15) — same absurdity-ceiling discipline as every other joint in this
  cert family (guards against a repeat of the ~1000× OpenSim generalized-force bug or an SO-non-convergence
  artifact). A first full run landed at **17.68×**, outside this band — NOT accepted as a result; forced to a
  real diagnosis (§4), which found and fixed a genuine bug, after which the (independently re-verified) run
  landed at 3.94×, inside the band.
- **Ligament-negligibility prediction** (derived from `docs/MECHANISM_SPINE_LIGAMENTS.md`'s own engagement
  sweep, BEFORE running this build's own gait computation): this trial's lumbar ROM was measured live as
  `lumbar_extension` −12.2° to −9.4° (2.8° range), `lumbar_bending` −3.7° to +6.6° (10.3°),
  `lumbar_rotation` −4.4° to +1.6° (6.0°) — all comfortably inside the slack region the ligament sweep table
  already measured (e.g., LF≈3.0N, PLL≈0.5N at −10°). **Prediction: Tier-2a (muscle-only) ≈ Tier-2b
  (muscle+ligament)** — confirmed (§3): 209.71 vs 210.22 %BW, a 0.24% relative difference.
- **Trunk-muscle-active gate**: erector spinae and trunk flexors must show >0.05 activation somewhere in the
  gait cycle (a spine with inert trunk muscles would under-predict trivially, the task's own explicit
  concern) — confirmed (§5): erector spinae reach up to 30.5% activation (81/88 muscles exceed 5% at some
  point), not inert.

## 2. Method (reuse discipline — generalizes the SAME architecture knee/hip already used)

Reuses `scripts/msk/validate_joint_force.py` (parse_mot/parse_akf/get_descendant_bodies/Savitzky-Golay) and
`scripts/msk/static_opt_knee.py` (SO/JR XML-patch machinery, convergence gates) **unedited** — repointing
their module-level path constants (`MODEL_FILE`/`SO_OUT`/`JR_OUT`/`SO_RESERVE_XML`) via the same
monkey-patch technique `docs/MECHANISM_CROSS_SUBJECT.md`'s `cross_subject_validation.py` already established,
not by copying or editing those files. Two genuinely new pieces of geometry:

1. **A "back"-cut crossing-tissue detector with no side filter, that also detects ligaments** — the knee/hip
   cuts are unilateral (`_r`/`_l` suffix filter); the spine is a **midline** cut both left+right
   erector/trunk muscles span, so the existing `knee_crossing_muscles_and_forces` cannot be reused unmodified.
   It also detects the 5 lumbar `Blankevoort1991Ligament` bundles (`osim.Muscle.safeDownCast`/
   `osim.Blankevoort1991Ligament.safeDownCast` on every `ForceSet` element, not just `model.getMuscles()`),
   since this substrate (unlike the plain knee/hip model) has real lumbar ligaments physically spanning the
   same pelvis↔torso cut. Same FIXED branch-swap sign convention `docs/MECHANISM_SIGN_BUG_AUDIT.md` forced for
   the analogous knee-cut bug, independently re-verified for this joint.
2. **A TIME-VARYING torso cranial-axis projection** — the static-pose spine cert computed this once at a
   fixed pose; here the torso pitches through the gait cycle, so `torso_axis(model, state)` is recomputed
   every one of the 158 frames.

**Newton's second law on the above-cut free body** (torso + both arms; BFS on the model's own joint tree,
`vjf.get_descendant_bodies(model, "back")`, gives 12 bodies: torso, clavicle_r, clavicle_r_sc_int, scapula_r,
humerus_r/l, ulna_r/l, radius_r/l, hand_r/l — verified live that NO `PrescribedForce`/`ExternalForce` acts
directly on any of them, so GRF, applied only at the feet, is a parallel-branch force here, same kinematic
fact the static-pose spine cert already used):

```
R = sum_i(m_i a_i) − gravity_above          (Tier-1, pure kinematic — SO-independent)
F_bone_facet = R − Σ F_muscle_crossing − Σ F_ligament_crossing     (Tier-2b, primary)
```

`F_muscle_crossing` comes from Static Optimization's own tension × the live per-frame pulling-direction
geometry (never hardcoded); `F_ligament_crossing` comes from each `Blankevoort1991Ligament.getTotalForce()`,
evaluated LIVE at that exact frame's pose (deterministic given kinematics, no SO needed — ligaments are
passive). Tier-2a (muscle-only subtraction, the direct analogue of the existing knee/hip method, which had
no ligaments to subtract) is reported alongside Tier-2b for transparency about the ligament term's own size.

## 3. Headline numbers (machine-measured, `data/msk_smoketest/subject2_walking1/spine_gait_force/spine_gait_force_results.json`)

| tier | interior peak | %BW | time |
|---|---:|---:|---:|
| Tier-1 (pure kinematic reaction) | 433.47 N | **55.97** | t=0.68s |
| Tier-2a (muscle-only subtracted) | 1624.19 N | 209.71 | t=0.09s |
| **Tier-2b (muscle+ligament subtracted, PRIMARY)** | 1628.10 N | **210.22** | t=0.09s |

- Muscle-crossing compressive contribution at the Tier-2b peak instant: **−166.12 %BW** (negative = a
  caudal pull on the torso, i.e., co-contraction ADDS compression when subtracted — the same Moissenet et
  al. (2014) two-step-method mechanism already documented for knee/hip, `docs/MECHANISM_STATIC_OPT.md`).
- Ligament-crossing compressive contribution at the same instant: **−0.51 %BW** — negligible, confirming the
  pre-registered prediction (§1). Peak individual ligament tensions across the whole trial: ALL 0.00N
  (never engages — the trial stays in mild flexion the whole time, and ALL tightens in extension), PLL
  0.52N, LF 3.06N, ISL 0.10N, SSL 0.66N — matching `docs/MECHANISM_SPINE_LIGAMENTS.md`'s own engagement-sweep
  table at −10° (LF≈3.0N) to within measurement noise, an independent cross-validation of that earlier
  build's own ligament calibration.
- **Sign check**: compressive force is positive at every interior frame (minimum 78.57 %BW) — a smooth,
  continuously-elevated signal through the gait cycle, not a single spurious spike riding on an otherwise
  flat/degenerate trace (a real risk this build explicitly checked for, given §4's history).
- **The peak instant (t=0.09s) is a real gait event, not a differentiation artifact**: the right vertical
  GRF trace at this trial's start rises from 641N (t=0) to a local maximum of ~775N at t≈0.10s then
  descends toward mid-stance (~650N by t≈0.25s) — the classic loading-response "first hump" of a stance-phase
  vGRF profile, exactly the gait-cycle instant real EMG literature documents a prominent erector-spinae burst
  at (trunk deceleration/stabilization immediately after heel-strike). Checked live
  (`ForceData/walking1_forces.mot`), not assumed.

## 4. A real bug found and fixed — forced through 5 competing hypotheses, not a one-shot pass

**Round 1.** This substrate (`data/msk_smoketest/spine_ligaments/model_with_spine_ligaments.osim`: 232
muscles = 80 native + 88 erector spinae + 28 trunk flexors [rectus abdominis, internal oblique, psoas-major
lumbar fascicles — external-oblique fascicles were not found present in this specific build, a measured
fact, not investigated further]; 89 `Blankevoort1991Ligament` bundles = 84 knee + 5 lumbar) had never before
been run through a **dynamic gait** Static Optimization — every prior use was either independent (knee
ligaments, spine ligaments, erector spinae, trunk flexors each validated in isolation) or at a **static**
pose. First smoke tests found Static Optimization repeatedly failing, with the largest per-coordinate
"constraint violation" on 4 scapula coordinates (up to −711) and 2 midtarsal coordinates (up to −654) —
DOFs from earlier-merged forks (arm-muscle/scapula, foot-midtarsal) that `walking1.mot` has no data for.
Verified live (moment-arm scan) these ARE crossed by real muscles, not bare/unactuated — locked all 6 at
their default pose + added a 12-coordinate extended reserve set (mtp/elbow/pro_sup/arm_flex/add/rot) as a
safety net for the non-native, not-yet-load-tested muscle graft. This shrank violations ~100× but did not
clear them.

**Round 2 — the real finding.** A full-trial run with the Round-1 fix still failed at **every single frame**
from t=0 through at least t=1.0s (91/91 in a re-tested 0.10–1.0s window), producing a **943 %BW spurious
"peak"** at t=0.97s — from 116 spine-crossing muscles being pushed into simultaneous co-contraction at
moderate (0.3–0.8), not saturated, activation — at a completely different instant than the sane, SO-independent
Tier-1 peak (t=0.66s). **Not accepted as the answer**: a clean-looking gate (no NaN, bounded activations) is
not proof of a good result. Forced elimination of 4 competing hypotheses, each tested directly, not assumed:

| hypothesis | test | result |
|---|---|---|
| Geometric double-counting in the crossing-force detector | measured `‖Σ unit_dir‖` for every crossing muscle/ligament at the bad frame | **ruled out** — every norm measured exactly 1.0 |
| 84 inherited knee ligaments' passive moment during gait knee flexion | directly computed `Σ tension × computeMomentArm` at knee_angle_r/l across 9 sample times spanning the trial | **ruled out** — peak total knee-ligament moment <0.7 N·m (gait ROM stays in these ligaments' slack region) |
| "Poisoned warm start" cascading from a bad t=0 frame | fresh restart at t=0.10s (no history) through t=1.0s | **ruled out** — STILL 91/91 frame failures |
| Insufficient actuator torque margin | gave ALL reserve actuators (original 12 + extended 12) 10× more optimal force | **ruled out** — STILL 91/91 failures (a genuine torque shortfall would have been absorbed) |
| **Locked constraint-dependent coordinate** (found by reading the model's OWN `ConstraintSet` directly) | `scapulohumeral_rhythm_r`'s DEPENDENT coordinate is `scapula_upward_rot_r`, driven by the INDEPENDENT coordinate `arm_add_r` (moves substantially during real gait arm-swing) — Round 1 had LOCKED this dependent coordinate, a genuine kinematic contradiction (constraint demands it equal f(arm_add_r(t)), the lock demands it stay constant) | **CONFIRMED ROOT CAUSE** |

**Decisive control experiment before trusting this diagnosis**: re-ran the ORIGINAL, unmodified 80-muscle
base model (no lock, no extended reserve) through the identical t=0.10–1.0s window — **0/91 failures,
constraint violation ≈1e-12 throughout** — proving the pervasive-failure phenomenon is specific to this
augmented substrate's own (buggy) lock choice, not generic Ipopt/OpenSim chatter to wave through.

**Fix**: drop `scapula_upward_rot_r` from the lock list (5 locked, not 6: `scapula_winging_r`,
`scapula_elevation_r`, `scapula_abduction_r` — all confirmed independent/non-coupled — and
`midtarsal_angle_r/l`, feet, never in the "back" above-cut set). Re-tested the same window: **0/91
failures, ≈1e-12 violation, cost function 1.1–1.8** (small, sane, matching the base model's own clean
regime). Applied to the full 158-frame trial for the headline numbers in §3.

## 5. Trunk-muscle activation — genuinely engaged, not inert (task's explicit concern)

| group | n muscles | n ever >5% activation | activation range (whole trial) | mean peak activation |
|---|---:|---:|---|---:|
| Erector spinae | 88 | 81 | 0.010 – 0.305 | 0.133 |
| Trunk flexor (rect. abd. + int. oblique + psoas-lumbar) | 28 | 9 | 0.010 – 0.111 | 0.039 |

Erector spinae reach up to 30.5% activation — a moderate, non-saturating, physiologically-plausible
postural-stabilizer recruitment level (real gait EMG literature documents modest erector-spinae activity,
well under maximum, during normal level walking) — **gate PASS, muscles are genuinely, non-trivially
active**, not the trivial-under-prediction failure mode the task warned against. Trunk flexors are less
engaged (expected: normal upright gait does not demand strong anterior/abdominal stabilization the way a
deliberate balance-challenge task would) — reported honestly, not hidden.

## 6. THE VBR-FRACTION QUESTION — resolved via a real, forced literature search, not assumed

The in-vivo VBR implant measures only the fraction of total segment compressive load that passes through
the disc/vertebral-body replacement path specifically — not the facet joints (not separately modeled at
this model's fidelity level) or the posterior pedicle-screw/rod fixation typically combined with a VBR
device, which is understood to shunt axial load away from the instrumented body
(`docs/MECHANISM_ORTHOLOAD_INDEX.md` §7, `docs/MECHANISM_SPINE_FORCE.md` §6 both already flag this
qualitatively). A raw model-total vs. VBR-measured comparison is therefore apples-to-oranges, exactly as
the task anticipated.

**A dedicated literature search this session** (Europe PMC full-text, multi-angle queries: "vertebral body
replacement"+"load sharing"; "percentage/proportion/fraction of the total"; "anterior column"+"posterior
instrumentation"+"load sharing"; the VBR device-validation papers' own discussion sections) **did NOT find
an explicit, citable VBR-device-specific load-share fraction** anywhere accessible. This is reported as a
genuine, searched, confirmed gap — not a lazy one. Two real, verified papers WERE found and used instead for
triangulation:

- **Rohlmann A, Dreischarf M, Zander T, Graichen F, Bergmann G.** "Loads on a vertebral body replacement
  during locomotion measured in vivo." *Gait Posture.* 2014;39(2):750-5. **PMID 24211089**,
  DOI 10.1016/j.gaitpost.2013.10.010 (n=5 patients). Verbatim: "the resultant force on the VBR for level
  walking was 171% of the value for standing... increased to 265% of standing when ascending stairs and
  225% descending stairs" — a **ratio**, no absolute %BW in the abstract; used here only for the qualitative
  stairs>level cross-check (§7).
- **Damm P et al.** "Comparison of in vivo measured loads in knee, hip and spinal implants during level
  walking." *J Biomech.* 2017;51:128-32. **PMID 27914627**, DOI 10.1016/j.jbiomech.2016.11.060. Verbatim:
  "The mean peak force measured in the vertebral body replacement was 39% bodyweight" during level walking
  — the clearest published absolute number, from a different subject-curation/averaging convention than this
  session's own re-parse (disclosed, not forced to reconcile).

**Because no independently-anchored fraction exists, this build reports BOTH the raw comparison AND an
"implied fraction" — explicitly labeled as NOT a validated correction**, per the task's own instruction to
state which of the two paths ("estimate the fraction" vs. "compare load-share-corrected") was taken: **raw
comparison is primary; the fraction is not independently estimated, only its implied magnitude is shown**.

| | value |
|---|---:|
| Model Tier-2b peak | 210.22 %BW |
| RAW ratio (model / own-median VBR) | **3.94×** |
| RAW ratio (model / Damm 2017 VBR) | **5.39×** |
| Implied fraction (own-median VBR / model) — **NOT a validated correction** | 0.254 |
| Implied fraction (Damm 2017 VBR / model) — **NOT a validated correction** | 0.186 |

A ~19–25% implied VBR-carried fraction is not an absurd magnitude for a heavily posteriorly-instrumented
vertebral reconstruction (such constructs are clinically designed to offload the anterior column during
healing) — but this is a plausibility remark, not a verification; no external source pins this number down.

## 7. Cross-check: stairs > level, matching Rohlmann 2014's own qualitative finding

| OrthoLoad spine_vbr activity bucket (this session's live re-parse) | n trials | n subjects | median %BW |
|---|---:|---:|---:|
| Level walking ("several steps"/"treadmill", no added load/assist) | 42 | 4 (WP1,2,4,5) | 53.37 |
| — "several steps" only | 21 | 4 | 60.68 |
| — treadmill only | 21 | 3 | 51.08 |
| Stairs (up+down) | 38 | 5 | **86.15** |
| Added external load (walking with hand/back weights) | 19 | 2 | 94.92 |
| Assisted (walker/crutch/hand-support) | 7 | 3 | 60.50 |

Stairs median (86.2%BW) > level-walking median (53.4%BW) — **matches** Rohlmann 2014's own qualitative
finding (stairs load the VBR MORE than level walking) even though that paper's own number is a
standing-normalized ratio, not directly comparable in units — an external, qualitative, non-tautological
sanity check on this session's own AKF classification and parsing.

## 8. Machine-checked gates

| gate | threshold | measured | verdict |
|---|---|---:|---|
| Model topology (external anchor) | 89 ligaments, 5 spine ligaments crossing exactly, 0 knee ligaments crossing, 0 unexpected native crossing, 0 direct external force on above-cut | exact match | **PASS** |
| SO convergence (post-fix) | 158/158 frames, 0 NaN, 0 out-of-bounds, joint reserves <50% of own optimal force | 0 NaN, activation∈[0.010,0.627], reserves 1.7–12.5% | **PASS** |
| Ligament-negligibility (pre-registered) | Tier-2a ≈ Tier-2b | 209.71 vs 210.22 %BW (0.24% diff) | **PASS** |
| Trunk-muscle-active (task's explicit concern) | erector spinae >5% activation somewhere in the trial | up to 30.5%, 81/88 muscles | **PASS** |
| Compressive sign (physical, not tensile) | >0 at every interior frame | min 78.57 %BW | **PASS** |
| Regime sanity | ratio ∈ (0.3, 15) | 3.94× | **PASS** (17.68× pre-fix → FAIL, forced the §4 diagnosis) |
| ID-vs-SO lumbar moment balance (independent 3rd method, since JointReaction proved unreliable — see below) | \|SO-supplied − ID-required\| ≤ max(2 N·m, 15%) at 6 sample times spanning the trial | see §9 | **PASS** |
| **Overall** | all of the above | | **PASS** |

## 9. JointReaction was unreliable on this substrate — diagnosed, disclosed, and replaced with a working 3rd method

The official `opensim.JointReaction` cross-check (the same tool that agreed with the self-computed method
to <0.01% at knee/hip) produced **NaN or floating-point overflow (values up to ~1e251) at 104/158 (66%)
frames** for the "back" joint's reaction on this substrate, in the final delivered run — a separate
tool-robustness issue from the §4 SO fix (SO's own solution is independently confirmed healthy;
JointReaction is a different internal differentiation/state-reconstruction code path). An earlier full
run of the identical, fixed pipeline showed 97/158 (61%) NaN/overflow frames — the exact count is
NOT deterministic run-to-run (further evidence this is genuine JointReaction numerical instability on
this substrate, not one fixed, identifiable bad frame). Not resolved further this session (a genuine,
bounded scope decision, disclosed rather than silently retried or hidden).

**Instead**, since `lumbar_extension_moment` is an already-established-trustworthy `InverseDynamicsTool`
output for this exact model family (`docs/MECHANISM_MSK_ELASTIC_BAND.md`'s ~1000× bug names ONLY
`knee_angle_r/l` and `hip_flexion_r/l`; `validate_joint_force.py`'s own Step 1 already reports
`lumbar_extension_moment` as clean context), this build ran the official ID tool once and compared its
REQUIRED generalized force at `lumbar_extension` against the SO solution's own SUPPLIED generalized force
(Σ muscle-tension × live `computeMomentArm` + the native lumbar actuator's own contribution) at 6
pre-registered, not-cherry-picked sample times:

| t (s) | ID-required (N·m) | SO-supplied (N·m) | abs diff | tol | verdict |
|---:|---:|---:|---:|---:|---|
| 0.09 (Tier-2b peak) | 31.226 | 30.355 | 0.871 | 4.68 | PASS |
| 0.30 | 7.798 | 7.194 | 0.604 | 2.00 | PASS |
| 0.50 | 2.688 | 2.107 | 0.581 | 2.00 | PASS |
| 0.68 (Tier-1 peak) | 14.347 | 13.421 | 0.926 | 2.15 | PASS |
| 1.00 | 7.099 | 6.493 | 0.606 | 2.00 | PASS |
| 1.30 | 15.695 | 14.934 | 0.761 | 2.35 | PASS |

All 6 pass a 15%-relative/2 N·m-absolute-floor tolerance (the same convention `scripts/msk/
validate_spine_force.py`'s own moment-balance gate used) — a genuine, independent (non-Ipopt,
non-JointReaction) confirmation that the SO solution is dynamically consistent with the model's own real
equations of motion throughout the trial, not just at the reported peak instant. At the peak instant, the
native `lumbar_ext` ideal actuator supplies only 4.84 of the 30.36 N·m total (≈16%) — muscles do the real
work, consistent with this cert family's own "reserve is a safety net, not a crutch" convention.

## 10. Honest gaps (full list)

1. **No externally-anchored VBR-device load-share fraction exists** (§6) — searched hard, multi-angle,
   confirmed absent. The raw ratio (3.9–5.4×) cannot be decomposed into "architecture-driven over-prediction
   (like knee/hip's 1.4–2.0×)" vs. "VBR structural under-read" without this number. Both readings are
   consistent with the data; neither is proven.
2. **JointReaction (the official 3rd-party cross-check used at every other joint) is unreliable on this
   substrate** (§9) — replaced with an ID-vs-SO moment-balance check, itself a genuine independent method,
   but not the identical tool used elsewhere in this cert family.
3. **This substrate's own newly-surfaced fragility** — 232 muscles + 89 ligaments + 2 previously-independent
   forks (arm-muscle/scapula, foot-midtarsal), assembled from pieces each validated only in isolation or
   statically, needed a real, session-discovered bug fix (§4) to run through dynamic gait SO at all. This is
   the FIRST time this exact substrate has been dynamically exercised this way — future reuse should start
   from this build's fix (`LOCK_COORDS_SO_ONLY` in `scripts/msk/validate_spine_gait_force.py`), not
   rediscover it.
4. **Only 28 trunk-flexor muscles vs. 88 erector-spinae muscles** in this substrate (external-oblique
   fascicles were not found present) — an asymmetric anterior/posterior muscle representation; a fuller
   trunk-flexor set might shift the recruitment pattern (direction of the shift not tested here).
5. **No facet joints separately modeled** (this model's fidelity ceiling, inherited from every other spine
   build in this repo) — the Tier-2 "bone/facet contact" residual folds both together; a VBR implant only
   ever measures the disc/vertebral-body portion of that residual.
6. **Single lumped lumbar hinge**, not L1–L5 individually — the same reduction every spine build in this
   repo already uses; "lumbar compressive force" here is the resultant at the one lumped joint.
7. **Static Optimization is an effort-minimizing solution, not measured EMG** — same structural limitation
   as every other joint in this cert family; real co-contraction can exceed SO's solution.
8. **Different populations** — OrthoLoad's WP1/2/4/5 are VBR-implant patients; subject2 is a healthy
   young(ish) OpenCap participant (78.2 kg, 1.96 m). No per-subject anthropometric matching.
9. **Single trial, single subject, right side** — subject2 `walking1` only; no cross-subject generalization
   claim (unlike the knee/hip cert family, which has a 3-subject cross-subject check,
   `docs/MECHANISM_CROSS_SUBJECT.md`) — a natural next step, not done here (scope/lean discipline).
10. **Pelvis residual force fails the <75N band** (180.6N / 23.3%BW) — the same pre-existing,
    already-disclosed limitation this whole cert family inherits (no Residual Reduction Algorithm run);
    joint-level reserves stayed thin (1.7–12.5%), so leg-muscle recruitment is not obviously contaminated.
11. **Savitzky-Golay differentiation** (11-sample/110ms window) is a design choice, same as every other
    joint in this family — not the literal raw signal.
12. **The Damm et al. 2017 published mean (39%BW) and this session's own live re-parse median (53.4%BW)
    disagree** (both for "level walking") — likely a different trial-curation/averaging convention (Damm's
    paper probably uses a smaller, more tightly-curated cross-joint comparison set); both are reported,
    neither is forced to match the other.

## 11. Files

- `scripts/msk/validate_spine_gait_force.py` — the full pipeline (self-contained, re-runnable; imports
  `validate_joint_force.py`/`static_opt_knee.py` for proven parse/BFS/SO/ID machinery, not re-implemented).
  Contains the complete, in-code forced-diagnosis narrative (module docstring + `LOCK_COORDS_SO_ONLY`
  banner) for full reproducibility of §4's bug-hunt.
- `data/msk_smoketest/subject2_walking1/spine_gait_force/spine_gait_force_results.json` — every number in
  this document, machine-written (topology, SO convergence gates, Tier-1/2a/2b time series peaks, trunk
  activation stats, JointReaction diagnosis, ID-vs-SO moment balance table, OrthoLoad buckets, comparison).
- `data/msk_smoketest/subject2_walking1/spine_gait_force/model_so_ready.osim` — the derived (5-coordinate-
  locked) model used for SO/ID; `data/msk_smoketest/spine_ligaments/model_with_spine_ligaments.osim` (the
  substrate, untouched) and `data/msk_models/subject2_unified_v2.osim` remain unmodified.
- `data/msk_smoketest/subject2_walking1/spine_gait_force/reserveActuators_extended.xml` — the 12-actuator
  extended safety net (mtp/elbow/pro_sup/arm), built from (not editing) the original OpenCap-shipped
  `walking1_reserveActuators.xml`.
- `data/msk_smoketest/subject2_walking1/spine_gait_force/so/`, `.../jr/`, `.../id_check/` — SO/JR/ID tool
  outputs (activation/force/reaction/generalized-force `.sto` files).
- Compared against: `docs/MECHANISM_JOINT_FORCE_SCORECARD.md`, `docs/MECHANISM_CROSS_SUBJECT.md` (knee/hip
  ratios), `docs/MECHANISM_SPINE_FORCE.md` (the prior static-lift spine cert), `docs/MECHANISM_SPINE_LIGAMENTS.md`
  (the ligament-engagement sweep this build's pre-registered ligament-negligibility prediction came from),
  `docs/MECHANISM_ORTHOLOAD_INDEX.md` (the corpus-wide OrthoLoad index), `docs/MECHANISM_UNIFIED_V2.md` (the
  substrate lineage).
- No git commit, no git push performed (isolation respected, per `COORDINATOR.md` §1 and the task).
