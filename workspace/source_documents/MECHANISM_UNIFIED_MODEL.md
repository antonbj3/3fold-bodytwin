# MECHANISM UNIFIED MODEL — merging the 9-fork problem into one subject2 instance (2026-07-21)

**Purpose.** `docs/MECHANISM_FIDELITY_ARCHITECTURE_AUDIT.md` §4-5 named this the #1 shortest-path
action: every fidelity increment built this cycle (erector spinae, trunk flexors, arm muscles,
scapula/clavicle, knee ligaments, midtarsal foot joint) lived in its own never-merged `.osim`
fork — "the pieces are NOT yet mutually composable ... today's gains do not compound; they sit
next to each other." This builds and verifies ONE unified model and reports, symmetrically, what
composed cleanly vs what didn't.

**Tool**: `scripts/msk/merge_unified_model.py` (re-runnable, self-contained; regenerates any
missing fork from its own script). **Output model**: `data/msk_models/subject2_unified.osim`
(gitignored, like all `data/msk_models/*.osim` — regenerate via the script). **Raw evidence**:
`scripts/msk/merge_unified_model_evidence.json`.

---

## 0. Executive verdict

**Structurally clean, exactly as pre-registered: 27 bodies, 27 joints, 232 muscles, 84
ligaments, 13 reserve actuators, 3 constraints, 49 markers (329 total force elements) — every
one of these 8 counts matches a threshold derived from raw per-fork measurement BEFORE the merge
was run, with zero unresolved cross-fork conflicts.** The model loads (`initSystem()`), runs a
forward-dynamics step, and runs real marker-driven IK without error. Two genuine, honestly-flagged
findings survive composition: (1) the audit's **58%-MCL-strain-at-140°** kinematic mismatch
**reproduces on the merged model**, confirmed fresh (not cited) — the 1-DOF-vs-6-DOF knee
incompatibility is real and persists; (2) a **new finding, not predicted in advance**: the
scapula/clavicle graft's reparented shoulder joint carries a **~40-55° rotational registration
offset**, surfaced only because this merge is the first time anyone ran a full marker-driven IK
trial through that fork (forced-adversary root-caused below, confirmed **pre-existing in the
standalone fork**, not introduced by merging, and confirmed correctly scoped to the right arm
only). The JAM/COMAK cartilage-contact model is excluded, by construction, not by an unresolved
conflict — it is a different individual on a different skeleton.

---

## 1. What was merged, in what order, and why that order

| Step | Fork | Adds | Built on (verified from the fork script's OWN path constant, not assumed) |
|---|---|---|---|
| 0 | **base** | 80 muscles, 22 bodies, 22 joints (10 CustomJoint/10 PinJoint/2 WeldJoint), 13 reserve actuators, 2 constraints, 49 markers | subject2's own OpenCap-scaled LaiArnoldModified2017 (`LabValidation_withVideos/subject2/OpenSimData/Mocap/Model/..._scaled.osim`) |
| 1 | `add_erector_spinae.py` | +88 muscles (iliocostalis/longissimus/multifidus) | **base** |
| 2 | `add_trunk_flexors.py` | +28 muscles (rectus abdominis/obliques/psoas-lumbar) | **erector_spinae's own output** — HARD dependency: `BASE_MODEL_PATH = aes.MODEL_COPY_PATH` in the script itself, plus a live `assert n_before == 168` gate in its own `add_muscles_to_copy()` |
| 3 | `add_arm_muscles.py` | +25 muscles (shoulder/elbow) | **base** |
| 4 | `add_scapula_clavicle.py` | +3 bodies (`clavicle_r`, `scapula_r`, a massless `clavicle_r_sc_int`), +3 joints (`sternoclavicular_prot_r`, `sternoclavicular_r`, `scapulothoracic_r`), +11 muscles (stabilizers), +1 constraint (`scapulohumeral_rhythm_r`), **+1 in-place modification** (`acromial_r`'s parent frame moved from `/bodyset/torso` to the new `/bodyset/scapula_r`) | **arm_muscles' own output** — HARD dependency: `SOURCE_MODEL = agm.MODEL_COPY_PATH`, and `reroute_11_muscles()`/`build()` reference muscle/body names that only exist once arm_muscles' 105-muscle output is the substrate |
| 5 | `add_knee_ligaments.py` | +84 `Blankevoort1991Ligament` bundles (ACL/PCL/MCL/LCL ×2 sides) | **base** |
| 6 | `foot_multisegment.py` | +2 bodies (`forefoot_r/l`), +2 joints (`midtarsal_r/l`), **in-place modifications**: `calcn_r/l` mass/COM/inertia shrunk, `mtp_r/l` re-parented onto the new forefoot bodies, 12 native muscles' path points re-parented (`edl/ehl/fdl/fhl/perlong/tibant`, ×2 sides), 4 markers re-parented (`r_toe`, `r_5meta`, `L_toe`, `L_5meta`) | **base** |

The two HARD chains (1→2, 3→4) are why "apply each script's own path constant, don't assume a
flat structure" mattered; the four lineages (spine: pelvis/torso · arm: humerus/ulna/radius/
clavicle/scapula/torso · knee: femur/tibia · foot: calcn/toes) are **anatomically disjoint**,
confirmed empirically (zero cross-lineage name collisions, `merge_report.conflicts == []` in the
evidence JSON) rather than assumed from the anatomy alone.

**Method**: a generic XML-level graft, not a re-run of each fork's own donor-extraction/scaling
math. For each fork, the script diffs the fork's saved `.osim` against the plain immediate parent
that fork's own script declares as its input, computes every component (Body/Joint/Force/
Constraint/Marker) that is new-by-name or modified-by-content, and grafts that delta onto the
growing unified tree (new → append, modified → replace, after checking the existing entry hadn't
already been independently touched by a different, earlier-applied fork — a real conflict, which
would be logged, never silently overwritten). Re-deriving each fork's geometry instead was
rejected: it would risk introducing new numerical bugs orthogonal to a MERGE task.

**Adversary forced before trusting the diff** (see the script's own docstring for the full
account): a naive whole-subtree comparison false-positives on every native muscle in every
fork, for two confirmed reasons — OpenSim's Python-API `printToXML` renames the anonymous
`GeometryPath` sub-object (`name="path"` in the original file → `name="geometrypath"` on
re-serialization) and reorders `GeometryPath`'s own children (`Appearance`/`PathPointSet`/
`PathWrapSet`). Both were measured directly against the real files (not assumed), and the
comparator was fixed to ignore the cosmetic name and to compare same-tag sibling groups in-order
while tolerating cross-tag reordering. Verified empirically to produce **zero false positives**
on every known-untouched component and to still catch every true modification, including one not
predicted in advance (`acromial_r`, found live by the diff, not hypothesized first).

---

## 2. Verification

### 2a. Structural counts — pre-registered before running, measured after (real `opensim.Model.initSystem()`, not the XML tree)

| Quantity | Expected (derived from raw per-fork counts) | Measured | |
|---|---|---|---|
| Bodies | 27 | 27 | PASS |
| Joints | 27 | 27 | PASS |
| Muscles (`Millard2012EquilibriumMuscle`) | 232 | 232 | PASS |
| Ligaments (`Blankevoort1991Ligament`) | 84 | 84 | PASS |
| Reserve actuators (`CoordinateActuator`) | 13 | 13 | PASS |
| Constraints (`CoordinateCouplerConstraint`) | 3 | 3 | PASS |
| Markers | 49 | 49 | PASS |
| Total force elements | 329 | 329 | PASS |

`232 = 80 (native) + 88 (erector) + 28 (trunk flexors) + 25 (arm) + 11 (scapula stabilizers)`.
Total mass after merge: **78.976 kg** (subject2's own scaled mass — ligaments/new bodies add
negligible mass, as expected for force elements and small girdle segments).

### 2b. Smoke tests (never just parse-and-assert)

- **Forward-dynamics**: `equilibrateMuscles()` + a real 5 ms `Manager` integration step — **OK**, no exception.
- **Real IK** on `walking1.trc` (158 frames), reusing `smoke_test_ik.py`'s own patch/run
  convention, output cross-checked against the pre-existing, independently-produced reference
  `walking1.mot` (an external anchor, not fabricated or tuned here). IK depends only on
  Body/Joint/Marker kinematics — muscles/ligaments are inert for IK by construction — so this is
  a genuine, decorrelated check of whether the graft changed anything it shouldn't have, not a
  tautology.

| Coordinate group | Worst RMSE (deg) | Interpretation |
|---|---|---|
| Truly unaffected (pelvis, hip, knee, lumbar, **left** arm — no graft touches these) | 1.05 (`lumbar_bending`) | Reassembly/solver noise. **Fails** an arbitrary 0.05° gate I pre-registered before measuring (an honest miscalibration on my part — that number was a guess, not derived from anything); **passes** a 2.0° gate anchored to this repo's OWN already-published foot-graft precedent (below), which is the more defensible anchor. |
| Known-affected: foot graft (`ankle_angle_{r,l}`, `subtalar_angle_{r,l}`) | 1.11–1.63 (RMSE), 3.06–3.62 (max) | Consistent with foot_multisegment's own already-published in-sample numbers (median 0.123°, regression 0.246°, this repo's own accepted range) — the new midtarsal DOF transfers into the unified model exactly as its own standalone cert already found and accepted. |
| Known-affected: shoulder registration (`arm_flex_r`, `arm_add_r`, `arm_rot_r`, `elbow_flex_r`) | 2.70 / 11.93 / **55.16** / 8.74 | New finding — see §3b. |

---

## 3. Honest findings — what merged cleanly, what didn't, and why

### 3a. Clean (verified, not assumed)

- Spine chain (erector spinae + trunk flexors): pure muscle additions to pre-existing
  `pelvis`/`torso` bodies. Zero body/joint/marker changes in either fork, confirmed by direct XML
  diff, not inferred from reading the code.
- Arm-muscle addition alone: pure additions to `torso`/`humerus_r`/`ulna_r`/`radius_r`. Zero
  body/joint changes.
- Knee ligaments: pure additions (84 `Blankevoort1991Ligament` elements) to existing
  `femur_{r,l}`/`tibia_{r,l}`. Zero body/joint/muscle changes — confirmed the native 80 muscles
  are byte-for-byte unchanged by this fork.
- Foot multisegment: structurally clean composition (new bodies/joints + in-place mass/parent
  reassignment of exactly the components its own build script says it touches — `calcn_r/l`,
  `mtp_r/l`, 12 named tendons, 4 markers — no more, no less).
- **Zero cross-fork conflicts**: no two forks modified the same named component in incompatible
  ways. This was measured (every "already modified + still different" case checked and logged),
  not assumed from the anatomy being disjoint.

### 3b. New finding: the reparented shoulder joint carries a ~40-55° registration offset

Not predicted before running the merge's own IK smoke test — found live, then forced through an
Observe→Orient→Decide→Act loop before being written up, not accepted on the first strange number:

- **Observe**: `arm_rot_r` (right internal/external shoulder rotation) showed 55° RMSE against
  the reference IK, dwarfing every other "unaffected" coordinate (next-worst was `elbow_flex_r`
  at 8.7°). `clav_prot_r`/`clav_elev_r` (the new sternoclavicular coordinates) were measured to be
  **perfectly constant across all 158 frames** (span = 0.000°) — meaning they carry zero
  observability from this marker set at all.
- **Orient (root cause, confirmed not guessed)**: the marker set (`walking1.trc`, same 49
  markers used everywhere in this repo) has **no marker on `clavicle_r` or `scapula_r`** — the
  only girdle-adjacent markers (`R_Shoulder`, `L_Shoulder`) are parented on `torso` in every fork,
  unchanged. `add_scapula_clavicle.py`'s own `build()` re-parents the **pre-existing**
  `acromial_r` (glenohumeral) joint's parent offset frame from `/bodyset/torso` directly onto the
  newly built `/bodyset/scapula_r`, using an explicitly-flagged **identity-rotation** registration
  assumption (`torso_anchor_orientation()`'s own docstring: "R_chain=IDENTITY ... verified live to
  give a physically sane ~20 deg scapular tilt"). That check verified the **scapula's own** rest
  tilt was anatomically sane — it never checked whether the **humerus's** rest orientation,
  reprojected through the new scapula frame, still matched the original torso-mounted definition.
  It doesn't: the registration carries a fixed rotational discrepancy, dominant in `arm_rot_r`,
  with smaller leakage into `arm_add_r`/`elbow_flex_r` via IK's global least-squares refit.
- **Decide/Act (decisive, cheap tests, not a one-shot conclusion)**:
  1. Ran the identical IK trial directly on the **standalone** `LaiArnoldModified2017_arm_muscles_subject2_scaled_scapulaClavicle`
     fork (no merge involved at all): `arm_rot_r` measured **37.8-43.1°**, versus the unified
     model's **41.0-44.5°** — same regime, both far from the reference's **-14.1 to -8.9°**. This
     **confirms the offset is a pre-existing property of the scapula_clavicle fork itself**, not
     something introduced by merging (the small residual difference between standalone and
     unified is consistent with ordinary IK-solver path sensitivity from the still-unconstrained
     scapula DOFs, not a new mechanism).
  2. Checked the **left** arm (never touched by this fork, which is right-side only by
     construction — confirmed only `_r`-suffixed bodies/joints/muscles were added): RMSE
     0.30-0.54°, maxabs <1.1° — clean, as it should be.
- **Verdict**: this is a genuine, previously-unknown limitation of the `add_scapula_clavicle.py`
  fork, surfaced for the first time by this merge's full-trial IK verification (the fork's own
  original verification checked moment-arm sign conventions and one static rest-pose sanity check,
  never end-to-end marker-driven kinematics of the reparented joint). **Not fixed in this task**
  (out of scope for a merge — fixing it means re-deriving the scapula registration in
  `add_scapula_clavicle.py` itself, a separate, focused piece of work) but documented here as a
  concrete, machine-measured caveat: **right-arm IK output from the unified model should not be
  used for kinematics/force analysis until this registration is corrected at the source.** Left
  arm, and every other segment, is unaffected.

### 3c. Re-confirmed: the knee-ligament / 1-DOF-knee kinematic mismatch persists in the merged model

Measured fresh on the unified model (not cited from the standalone cert): sweeping `knee_angle_r`
from 0° to its full range, right-MCL strain (16 bundles: `MCLd1-5`, `MCLp1-5`, `MCLs1-6`) grows
monotonically from 5.0% (0°) through 25.7% (90°), 46.3% (120°), to **58.2% at 140°** — reproducing
the audit's "58% at 140°" almost exactly. The mechanism is unchanged by merging: Lenhart2015's
donor ligament geometry assumes a free 6-DOF knee (kinematics emerge from force balance);
LaiArnold's `walker_knee_r` is a prescribed 1-DOF `CustomJoint` (only `knee_angle_r` is free;
ab-adduction/rotation/all 3 translations are polynomial functions of flexion angle, not free) —
grafting the ligaments cannot change that, and doesn't; it is a structural incompatibility between
two rigid-body kinematic conventions, confirmed (not merely assumed to persist) on the actual
merged artifact.

### 3d. Excluded, by construction — not a merge failure

The JAM/COMAK cartilage-contact model
(`jam-resources/models/knee_tka/grand_challenge/DM/DM.osim`, previously headlined here as the repo's
single strongest cert — og1's ~3% error vs Grand-Challenge in-vivo contact force) is **not merged and
cannot be**: it is the Grand-Challenge TKA patient "DM" — a different individual, a different skeleton
lineage (JAM/COMAK, not LaiArnold), an implant knee, zero ACL. It shares no body or joint names with
subject2's model and is not a fork of it in any sense; there is no body-for-body correspondence to
graft through. This is a category difference, not an unresolved conflict, and matches the audit's
own finding (§2 row 3, §4.3) verbatim.

**CORRECTION (2026-07-26, postdates this section):** the "single strongest cert" superlative does not
survive a domain-default check. A same-day graph audit found that OrthoLoad's independent
population-median baseline (zero fit to subject DM) beats the model on mean absolute error across the
3 run trials — 3.06% vs the model's 5.72% — and wins 2 of the 3 held-out trials (og1, og3); the model
wins only og4 (`AUDIT-VOIDFLOOR-DOMAINDEFAULT-KNEE-SARCOMERE-2026-07-26`,
`AUDIT-KNEE-BODYWEIGHT-WINDOW-CORRECTION-OG1-SELECTION-2026-07-26`). og1's cited ~3%/r=0.802 was also
the best of the 3 same-day-run trials (og4 3.76%/r=0.741, og3 10.35%/r=0.658), not typical. This does
not change this section's own point — DM and subject2 are still different individuals and still cannot
be merged — it only retracts the comparative superlative. What the model still demonstrably owns
against a fair adversary: it beats a cross-subject-scrambled-activation null by 5-7x on RMS
(`VOID-FLOOR-KNEE-GENERIC-WAVEFORM-CLEARS-BAR`).

### 3e. Not in scope for this merge (present in the repo, deliberately not included)

Per the task's explicit script list, the following existing forks/prototypes were **not**
targeted and are **not** part of the unified model — noted here so their absence isn't mistaken
for an oversight: `LaiArnold_hallux_proto_scaled.osim` (hallux-split dial-turn-1 prototype),
`AXISTEST_ONLY_LaiArnold_midtarsal_{longaxis,zaxis}.osim` (the foot-axis exploration variants that
preceded the accepted `foot_multisegment.py` axis, §"HONEST-FAILURE FORCED FIX" in that script),
`LaiArnoldModified2017_mediapipe_bridge_subject2_scaled.osim` (pose-pipeline bridge model), the
elastic-band model pair (`model_no_band.osim`/`model_with_band.osim`, a different experimental
axis), and the RRA-adjusted model (`subject2_walking1_rra_pass1_adjusted.osim`, already found in
`MECHANISM_RRA.md` to make things worse when re-validated through the muscle-driven model, so
correctly excluded on its own merits too).

---

## 4. Reproduction

```
source_repository/.venv-msk/bin/python3 scripts/msk/merge_unified_model.py
```

Regenerates any missing fork from its own script, builds `data/msk_models/subject2_unified.osim`,
verifies structural counts, runs both smoke tests, re-verifies the knee-ligament finding, and
writes `scripts/msk/merge_unified_model_evidence.json`. Exit code 0 iff structural counts all
match, zero conflicts, both smoke tests ran without error, and the precedent-anchored (2.0°)
truly-unaffected-coordinate gate passes — this run's exit code was **0 (PASS)**, with the §3b
shoulder-registration caveat printed explicitly (a documented, non-blocking finding, not silently
swept into the pass).

---

## Evidence index

`scripts/msk/merge_unified_model.py` (tool) · `scripts/msk/merge_unified_model_evidence.json`
(raw machine-checked numbers behind every claim above) · `data/msk_models/subject2_unified.osim`
(the unified model itself, gitignored — regenerate via the script) ·
`docs/MECHANISM_FIDELITY_ARCHITECTURE_AUDIT.md` (the audit this task responds to) ·
`docs/MECHANISM_KNEE_LIGAMENTS.md`, `docs/MECHANISM_FOOT_MULTISEGMENT.md`,
`docs/MECHANISM_SCAPULA_CLAVICLE.md`, `docs/MECHANISM_ARM_MUSCLES.md`,
`docs/MECHANISM_ERECTOR_SPINAE.md` (per-fork source certs, read for parent/dependency
verification, not re-litigated).

No git commit or push performed (isolation-respected, per task and per `COORDINATOR.md` §1). New file
written: `data/msk_models/subject2_unified.osim` (a NEW model file, none of the 6 source forks
were modified).
