# MECHANISM SHOULDER FORCE — replicating the knee force-cert method for the glenohumeral joint (2026-07-21)

Executes the operator's task: replicate the VALIDATED knee force-cert method
(`docs/MECHANISM_JOINT_FORCE_VALIDATION.md` + `docs/MECHANISM_STATIC_OPT.md`) for the
shoulder. Every number below is machine-measured this session
(`scripts/msk/validate_shoulder_force.py`, exit 0, re-run twice with identical output), not
recalled. Isolation respected: `.venv-msk` only, OrthoLoad + LabValidation data read in
place (never written to), no git commit/push.

## Headline result

| | |
|---|---:|
| **Model has glenohumeral DOFs?** | **YES** — `acromial_r`/`acromial_l`, 3-DOF ball joint, +/-180 deg range each |
| **Model has deltoid/rotator-cuff (or ANY arm) muscles?** | **NO — 0 of 80** (0 of 168 in the erector-spinae variant) |
| **Static-Optimization muscle-driven solve possible?** | **NO** (honest finding, forced via 2 independent enumeration methods — see §1) |
| PREDICTED shoulder reaction force (pure kinematics, synthetic 0→90° abduction, T=2.0s) | **5.31 %BW** |
| IN-VIVO OrthoLoad shoulder (median, n=23 unloaded 90° abduction/elevation trials, 6 subjects) | **72.2 %BW** |
| ratio | **0.073** (predicted / in-vivo) |

**This is an honest-limitation finding, not a failed replication attempt.** The knee method
has two stages: (1) a pure-kinematics Newton's-law reaction force, then (2) Static
Optimization adding real muscle forces to close most of the gap. Stage (1) replicates
cleanly onto the shoulder (below). Stage (2) **cannot run on this exact twin model** — not
because the solve failed, but because the model has no shoulder muscles to solve *for*. A
second, independent honest gap compounds this: no real arm-elevation motion-capture trial
exists anywhere in this twin's own corpus (checked across all 10 available OpenCap
subjects), so the input kinematics used below are an explicitly-synthetic, prescribed
trajectory, not measured mocap. Both gaps are disclosed and load-bearing throughout, per
§7.

## 1. Model audit (the primary, load-bearing finding) — forced via 2 independent methods

`LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim` — **the exact same model file**
the knee cert validated (subject2, 78.2 kg) — has:

- **`acromial_r`/`acromial_l`**: a 3-DOF `CustomJoint` (parent=`torso`, child=`humerus_r`/`l`),
  coordinates `arm_flex_r`, `arm_add_r`, `arm_rot_r` (each range −180°..+180°, default 0°) —
  a standard "combined humerothoracic" ball-joint shoulder representation (the same
  simplification used in e.g. the Rajagopal 2016 full-body model: **no separate
  scapula/clavicle body**, so true scapulohumeral rhythm — which contributes roughly a
  third to half of total arm elevation in real shoulders — is not separately resolved;
  this joint's angle is the *net* humerus-relative-to-thorax angle).
- `elbow_r`/`l` (PinJoint), `radioulnar_r`/`l` (PinJoint, pronation/supination),
  `radius_hand_r`/`l` (**WeldJoint** — the hand is rigidly fused to the forearm, matching
  the model's own "weldHand" name; no wrist DOFs).
- **Zero muscles** attach to any arm body (`humerus_r/l`, `ulna_r/l`, `radius_r/l`,
  `hand_r/l`) — no deltoid (anterior/middle/posterior), no rotator cuff (supraspinatus,
  infraspinatus, subscapularis, teres minor), no pec major, no lat dorsi, no teres major,
  no coracobrachialis, no biceps/triceps long head. All 80 of the model's muscles are
  leg/trunk muscles.
- The entire arm is actuated **only** by 10 ideal `CoordinateActuator`s (one per DOF, both
  sides: `shoulder_flex/add/rot_{r,l}`, `elbow_flex_{r,l}`, `pro_sup_{r,l}`, each
  `optimal_force`=10 N·m) — i.e. **properly determined, not redundant**.

**Why this rules out Static Optimization specifically** (not just "makes it harder"):
Static Optimization exists to resolve *muscle-vs-muscle redundancy* — the same net joint
moment can be produced by infinitely many activation combinations when multiple muscles
cross a joint, and antagonist co-contraction (two muscles pulling opposite ways on the net
moment but the *same* way on joint compression) is exactly the mechanism the knee's SO
step used to close its gap (`docs/MECHANISM_STATIC_OPT.md` §7: rectus femoris +
gastrocnemius co-contraction). With **exactly one ideal actuator per shoulder DOF**, there
is no redundancy to resolve and — more fundamentally — **no mechanism for co-contraction
to exist at all**: a single ideal torque motor cannot compress a joint the way two
antagonist muscles pulling in different net directions can. Running
`opensim.StaticOptimization` on this joint would return one algebraic solution
(`activation = required_torque / optimal_force`, no Ipopt search, no redundancy) — a
**vacuous tautology**, not a muscle-force estimate. It was deliberately **not run** (a
vacuous "pass" is worse than an honest gap).

**This negative was forced to its strongest form before being accepted** (symmetric QC — a
kill claim carries the same burden as a confirmation):
1. Cross-checked across **all 4** LaiArnold model variants present in this repo (base
   80-muscle model, the erector-spinae variant at 168 muscles, the mediapipe-bridge
   variant, the hallux-detail variant) — **0/4 have any arm-crossing muscle**, not a
   single-file coincidence.
2. Re-enumerated muscle-crossing via a **second, independent** OpenSim code path —
   `GeometryPath.getCurrentPath(state)` (resolves wrapping/conditional/moving path points
   at the actual current pose) instead of the static `getPathPointSet()` definition — **0
   mismatches** across all 80 muscles between the two methods.
3. Censused **every** Force element type in the model (not just classes named `*Muscle`):
   only `Millard2012EquilibriumMuscle` (80) and `CoordinateActuator` (13) exist — no
   `Ligament`, `PathSpring`, or other passive-tissue element that could contribute passive
   joint compression is present anywhere in the model.
4. Confirmed the body list has **no scapula/clavicle body at all** (21 bodies total:
   pelvis, 2×(femur/tibia/patella/talus/calcn/toes), torso, 2×(humerus/ulna/radius/hand)) —
   there is no hidden body a muscle could attach to that this check would have missed.

## 2. Corpus check — no real arm-elevation trial exists (checked, not assumed)

Every OpenCap subject's IK trial directory was enumerated live
(`/media/anton/8838D60F38D5FBDE/mechanism_data/LabValidation_withVideos/subject{2,3,4,5,6,7,8,9,10,11}/OpenSimData/Mocap/IK/`):
all trials across all 10 subjects are `DJ*` (drop-jump), `STS*`/`STSweakLegs*`
(sit-to-stand), `squats*`, or `walking*`/`walkingTS*` — **zero** trial names match an
arm/shoulder/elevation/abduction keyword search. This twin's own motion corpus has no
deliberate arm-elevation ROM capture. (Consistent with the corpus's purpose: it is a
gait-lab dataset; arm coordinates exist in the model only for whole-body kinematic/inertial
completeness during walking, e.g. arm-swing counter-rotation.)

## 3. What IS computed: pure-kinematic Newton's-law shoulder reaction force

Since Stage-2 (muscle-driven) is unavailable and Stage-1 needs no real mocap trial to be
meaningful (see geometry below), the same Newton's-law method the knee cert used for its
own pre-SO baseline (`validate_joint_force.get_descendant_bodies` +
`savgol_smooth_and_derivs`, **imported, not reimplemented**) was ported to the shoulder:

```
sum_i(m_i * a_i) = gravity + R          (bodies distal to acromial_r; R = shoulder reaction)
```

**Geometric simplification vs. the knee/hip**: the free body distal to the shoulder
(`humerus_r`+`ulna_r`+`radius_r`+`hand_r`) never touches the ground — unlike the knee/hip
cut, which required the *measured* GRF as an essential external-force input, the shoulder
cut needs **no GRF term at all**. This makes the arm case a *cleaner*, lower-assumption cut
than the leg, even though the input kinematics themselves are synthetic (see §2).

**Input kinematics**: a prescribed 0°→90° abduction of `arm_add_r`, using a quintic
"smootherstep" (zero velocity **and** zero acceleration at both ends — a realistic slow
controlled-test profile, not an impulsive start/stop), padded with 0.3 s of genuinely-static
hold at each end. All other coordinates held at the model's own default value (elbow stays
at 0° = full extension, the standard clinical abduction-test posture). The abduction sign
convention (`arm_add_r` negative = abduction) was **measured live** via forward kinematics
(hand_r moves 0.606 m laterally at `arm_add_r=−30°` vs. 0.122 m at `+30°`), not assumed.

### Duration sensitivity sweep (the dominant free assumption — swept transparently, not tuned)

No real motion exists to time-match, so the assumed elevation *speed* is the single most
consequential free parameter here (unlike the knee's smoothing-window sweep, which is a
numerical nuisance parameter, this one is a real modeling assumption) — reported as a full
sweep, not a cherry-picked value:

| assumed duration (0→90°) | shoulder peak (%BW) | elbow peak (%BW, nesting check) | +2kg hand-load peak (%BW) |
|---:|---:|---:|---:|
| 1.0 s (fast) | 6.469 | 3.311 | 10.669 |
| 1.5 s | 5.607 | 2.705 | 8.895 |
| **2.0 s (primary)** | **5.306** | **2.493** | **8.274** |
| 3.0 s | 5.090 | 2.341 | 7.830 |
| 4.0 s (slow) | 5.015 | 2.288 | 7.675 |

Shoulder ≥ elbow at every duration (mass distal to the shoulder > mass distal to the elbow —
the same "each proximal cut sheds that segment's own weight" nesting signature the knee
cert reported). As assumed speed decreases, the result monotonically converges toward the
closed-form floor below — the expected physical signature (inertial contribution shrinks
as motion slows), not a coincidence.

## 4. Internal-validity gates (machine PASS/FAIL — no real-mocap anchor exists, so these
replace the knee's GRF-sanity/whole-body-residual checks with checks appropriate to a
*synthetic* input)

| gate | pre-registered threshold | measured | verdict |
|---|---|---:|---|
| Endpoint-rest check (accel. in the genuinely-static padding, vs. peak) | ratio < 0.30 | 0.0000 at all 5 durations | **PASS** |
| Closed-form quasi-static floor convergence (slowest duration vs. analytic R=m·g) | diff < 1.5 pct-BW | 5.015 vs. 4.918 %BW (diff 0.097) | **PASS** |
| Smoothing-window stability (windows 70–170ms, at T=2.0s) | range < 1.0 pct-BW | 0.000 (5.306%BW at every window) | **PASS** |
| OrthoLoad ballpark vs. task's own stated range ("tens-to-~100%BW") | 20–200 %BW | median 72.2%BW | **PASS** |
| Regime sanity: predicted/in-vivo ratio | in (0.01, 0.6) — positive, well below the muscle-inclusive anchor | 0.073 | **PASS** |

The endpoint-rest check is a synthetic-control, known-answer validation: the padded
start/end regions are *exactly* constant by construction, so the Savitzky-Golay pipeline
recovering *exactly* zero curvature there (not just "small") is the correct, expected
outcome — confirming the differentiation machinery before trusting its peak-region output.
The closed-form floor is an independent derivation (Newton's law at zero velocity/
acceleration reduces exactly to R=m·g, regardless of arm angle) — a genuine
over-determination check against the numerical pipeline, not a tautology (it uses no
Savitzky-Golay machinery at all). Shoulder-chain mass = 3.846 kg = 4.92 %BW of the model's
78.2 kg total — this is the physical floor any slow-enough elevation must approach.

## 5. OrthoLoad in-vivo shoulder anchor — selection pre-registered from the task's own wording

The AKF parser (`validate_joint_force.parse_akf`, format-verified identical for shoulder:
same `Time/Fx/Fy/Fz/F/Mx/My/Mz/Marker` layout as knee — cross-checked live on a sample file,
`sqrt(Fx²+Fy²+Fz²)` matches the file's own printed `F` to 0.011 N) was re-run directly
against all 140 shoulder AKF files (self-contained — does not depend on the pre-built
`_index/` continuing to exist, though it is corroborated by it). Filter — matching the
task's own stated activity ("arm elevation/abduction to 90°"), decided from the corpus's
own vocabulary before viewing the resulting ratio: Comment #1 contains "Standing Position"
+ ("Abduction"|"Elevation") + "No Weight" + "90 deg".

- **Primary (unloaded, 90°): n=23 trials, 6 distinct subjects** (S1R, S2R, S3L, S4R, S5R,
  S8R) — min 36.3, **median 72.2**, mean 78.2, max 128.1 %BW. Representative
  (closest-to-median) trial: `s4r_140207_1_73`, 72.21 %BW, BW=500N, "Shoulder Joint;
  Standard Test; Standing Position; Elevation; No Weight; 90 deg."
- **Secondary/bonus (2 kg load, 90°): n=3 trials, 3 subjects** — min 130.5, median
  **150.6**, mean 161.3, max 202.7 %BW. (Predicted +2kg case: 8.27 %BW at T=2.0s → ratio
  0.055 — same qualitative story: the in-vivo number jumps by ~2× under a 2kg hand load,
  the muscle-free prediction barely moves, because it has no co-contraction mechanism to
  amplify.)
- Both subsets pass the AKF resultant-column cross-check (recomputed
  `sqrt(Fx²+Fy²+Fz²)` vs. the file's own printed `F`): max abs diff 0.012 N (primary) / 0.012
  N (secondary) across all trials — same tightness as the knee cert's own cross-check.
- The in-vivo median (72.2 %BW) is consistent with the published OrthoLoad shoulder
  literature's general order of magnitude for unloaded arm-abduction tasks (Westerhoff et
  al., Bergmann et al. — tens to ~100 %BW), and with the task's own stated ballpark.

## 6. Ratio interpretation — why 0.073 is expected, not alarming

The predicted number is **~14× smaller** than the knee's own pre-SO ratio (0.396,
`docs/MECHANISM_JOINT_FORCE_VALIDATION.md`), and this is the *correct* qualitative signature
of two compounding, well-understood effects, not a bug:

1. **100% of the muscle/co-contraction contribution is missing here**, vs. the knee's case
   where a nonzero net reaction already existed before SO and muscles amplified/refined it.
   Here there is no muscle path to add *any* co-contraction term back in — the entire
   72.2 %BW in-vivo signal is, by construction of this model's actuation, unreachable.
2. **The arm is a small fraction of body mass** (4.9%) relative to the leg (the knee's
   distal-chain mass is a much larger fraction of body weight), so even the *kinematic*
   floor is proportionally tiny — real in-vivo shoulder loads during arm elevation are
   understood in the literature to be almost entirely deltoid/rotator-cuff-generated, far
   more muscle-dominated than knee loads during stance.

A ratio near or above 1 would have been the actual red flag here (impossible for a
zero-muscle, pure-gravity prediction to legitimately meet or exceed a muscle-inclusive
in-vivo anchor) — 0.073 sitting comfortably inside the pre-registered (0.01, 0.6) sane band,
and tracking the closed-form 4.9% floor closely, is exactly what a correctly-implemented,
honestly-scoped kinematic-only method should produce on this model.

## 7. Honest caveats (full list)

1. **No muscle-driven solve exists for this joint on this model** (§1) — the primary,
   load-bearing caveat. This is not a tuning/convergence failure like the knee's own
   documented SO fragility (`bt_memory/knee-comak-flagship-close-and-fork-build.md`); it is
   a structural absence (0 actuatable redundancy at the shoulder).
2. **No real arm-elevation mocap trial exists in this twin's corpus** (§2) — the input
   kinematics are an explicitly-synthetic, prescribed trajectory, not measured motion. The
   knee cert used 100% real force-plate-driven mocap; this shoulder result does not.
3. **Assumed elevation speed is a free, unmeasured parameter** (§3) — swept transparently
   (1.0–4.0 s) rather than reported as a single unqualified number; the true speed of any
   real subject's controlled abduction is unknown here.
4. **Simplified shoulder joint**: no separate scapula/clavicle body — scapulohumeral rhythm
   is not resolved; `arm_add_r` is the net humerothoracic angle, not pure glenohumeral
   rotation. A more anatomically detailed shoulder model (e.g. one with a
   scapulothoracic/sternoclavicular chain) would decompose this differently.
5. **No hand/wrist DOFs** (`radius_hand_r` is a `WeldJoint`) — irrelevant to this specific
   abduction task (elbow held extended) but limits what other arm activities could be
   analyzed with this model.
6. **Different populations**: OrthoLoad's shoulder subjects are instrumented
   hemiarthroplasty/reconstruction patients; the twin (subject2) is a healthy OpenCap
   participant — the same population caveat the knee cert carried forward, unchanged.
7. **Rigid hand-load approximation** (§3 bonus case): the 2 kg load is modeled as a point
   mass co-moving with `hand_r`'s own COM acceleration — a reasonable but simplified
   rigid-attachment assumption.
8. **Single BFS-derived free-body cut, right arm only** — same scope caveat as the knee
   cert's own single-trial, single-side scope.

## 8. What a proper shoulder solve would need (the honest next step)

To genuinely replicate the knee's Stage-2 (muscle-driven, co-contraction-inclusive) result
for the shoulder, this twin would need a model with:

1. **Deltoid** (anterior/middle/posterior heads) and the **rotator cuff** (supraspinatus,
   infraspinatus, subscapularis, teres minor) at minimum — the primary glenohumeral
   stabilizer/compressor muscles absent here entirely.
2. Ideally also pec major, lat dorsi, teres major, coracobrachialis, and biceps/triceps
   long heads (biarticular, cross both shoulder and elbow) for a physiologically complete
   redundant actuator set Static Optimization could meaningfully resolve.
3. A separate scapula (and ideally clavicle) body with a scapulothoracic joint, so
   scapulohumeral rhythm is resolved rather than folded into a single combined angle —
   standard in richer upper-limb models (e.g. the MoBL-ARMS/Saul 2015 lineage, or
   Holzbaur 2005), none of which are currently vendored in this repo or its mounted
   corpus (checked: `arm26.osim`, present in the OpenSim-JAM build tree, is a 2-muscle
   elbow-only model with no shoulder muscles either — not a substitute).
4. A real arm-elevation motion-capture trial (with or without a hand load) for at least one
   subject, to replace the synthetic prescribed trajectory used here.
Grafting such a model onto this twin's own scaled anthropometry (rather than using a
generic/uncalibrated donor model) is a real, bounded, multi-step follow-on task — not
attempted here, consistent with keeping this session's scope to what the existing twin
model and corpus can honestly support.

## Files

- `scripts/msk/validate_shoulder_force.py` — the full pipeline (self-contained,
  re-runnable; imports `validate_joint_force.py` for the proven parse_akf/savgol/BFS code,
  not reimplemented).
- `data/msk_smoketest/subject2_arm_abduction/shoulder_force_validation/shoulder_force_validation_results.json`
  — every number in this document, machine-written.
- `data/msk_smoketest/subject2_arm_abduction/shoulder_force_validation/synthetic_abduction_T2.0s.mot`
  — the primary-duration synthetic kinematics, written for auditability.
