# MECHANISM UNIFIED MODEL v2 — merging the hand + hip/ankle-ligament forks (2026-07-21)

**Task.** `docs/MECHANISM_LAYER_COMPLETENESS.md` flagged that fragmentation GREW: the unified
model does not contain the hand build or the hip/ankle ligaments, both still living in further,
separate, unmerged forks built the same evening. This closes that gap: `full_hand.py`'s output
(un-welded wrist + index+middle fingers + 2 intrinsics) and `add_hip_ankle_ligaments.py`'s
output (22 Blankevoort bundles) are grafted onto the CURRENT unified model.

## 0. Executive verdict

**PASS, symmetric.** Both forks merge with **zero cross-fork conflicts** and **zero unexpected
component loss**. Structural counts match a pre-registered "sum of parts" exactly (34 bodies /
34 joints / 240 muscles / 111 ligaments / 364 total force elements). The merged model
`initSystem()`s cleanly, runs a 5 ms forward-dynamics step with all-finite state, and runs real
marker IK on `walking1.trc` — every PRIOR certified finding this merge could plausibly disturb
(the knee's 58.2% MCL-strain mismatch, the hip/ankle ligaments' engagement-direction tests, the
hand's textbook lumbrical double-action) **reproduces on the composed artifact**, most to
near-bit-identical precision. One genuinely NEW piece of merge machinery was required (not
present in the original `merge_unified_model.py`) because the hand fork does something none of
the original 6 forks did: it **removes** a component (the old 0-DOF wrist weld). A naive reuse
of the original apply-only merge logic was forced and found to produce a real, non-obvious
defect — documented in detail in SS3b, since it is the one part of this task that was NOT clean
on the first attempt.

## 1. What was merged, and onto which file

**"The current unified model" is NOT `data/msk_models/subject2_unified.osim`.** That file's
mtime (18:08) pre-dates `add_spine_ligaments.py` (19:34) — confirmed via `stat`, not assumed.
Per that script's own isolation note, it read `subject2_unified.osim` READ-ONLY and wrote a NEW
file, so the true current 7-component substrate is:

```
data/msk_smoketest/spine_ligaments/model_with_spine_ligaments.osim
```
measured live (real `initSystem()`) both before writing `merge_unified_v2.py` and again at the
top of its `build_unified_v2()` (this repo is concurrently written by a second mechanism
instance, COORDINATOR.md SS0 — the second measurement guards against drift): **27 bodies / 27 joints
/ 232 muscles / 89 ligaments (84 knee + 5 spine) / 13 reserve actuators / 3 constraints / 49
markers / 334 total force elements**, matching pre-registration exactly on this run (zero
drift detected).

Two forks were grafted onto it, via the SAME generic XML-delta/graft machinery as
`merge_unified_model.py` (imported and reused, not re-derived):

| Fork | Its own declared parent | What it adds |
|---|---|---|
| `full_hand.py` output (`LaiArnoldModified2017_full_hand_subject2_scaled.osim`) | The plain base (`...weldHand_scaled.osim`) — full_hand.py loads `anatomical_hand.py`'s own output directly and extends it, so diffing full_hand's FINAL output against the plain base in one step captures BOTH stages combined, matching the task's own framing of "anatomical_hand.py -> full_hand.py" as one fork | +7 bodies (index_proximal/medial/distal_r, middle_proximal/medial/distal_r, wrist_int_r), +8 joints (wrist_flex_r_joint, wrist_dev_r_joint, mcp2/pip2/dip2/mcp3/pip3/dip3_flex_r_joint) **and REMOVES 1** (`radius_hand_r`, the old 0-DOF weld), +8 muscles (FDP2_r, FDS2_r, ED2_r, FDP3_r, FDS3_r, ED3_r, LUM1_r, DI1_r) |
| `add_hip_ankle_ligaments.py` output (`model_with_hip_ankle_ligaments.osim`) | The plain base (loaded directly — confirmed via its own `MODEL_PATH` constant) | +22 `Blankevoort1991Ligament` force elements only (ILFL_lat/med x2 sides, PFL, ISFL x2 sides = 8 hip; ATFL, CFL, PTFL, DELT_tibcalc/tibnav/tibspring/posttibtal x2 sides = 14 ankle). Zero Body/Joint/Constraint/Marker touches — same pure-ForceSet-addition pattern as the already-merged knee ligaments. |

Output: **`data/msk_models/subject2_unified_v2.osim`** (NEW file — `subject2_unified.osim` and
`model_with_spine_ligaments.osim` are untouched; verified by `stat` mtime/size before and after
this run).

## 2. The CRITICAL check the task asked for: hand-fork vs arm-muscle/scapula-fork collision

The hand fork un-welds the wrist; the arm-muscle/scapula-clavicle forks (already in the unified
base) also modify the arm chain. Measured directly (`compute_delta` on every `SET_NAME`, cross-
referenced against scapula_clavicle's own already-applied delta, read from
`merge_unified_model_evidence.json`, not re-derived):

| SET_NAME | hand fork touches | arm/scapula fork touches (already merged) | overlap |
|---|---|---|---|
| BodySet | index_\*_r (x3), middle_\*_r (x3), wrist_int_r | clavicle_r, clavicle_r_sc_int, scapula_r | **none** |
| JointSet | dip2/dip3/mcp2/mcp3/pip2/pip3_flex_r_joint, wrist_flex_r_joint, wrist_dev_r_joint (+removes radius_hand_r) | scapulothoracic_r, sternoclavicular_prot_r, sternoclavicular_r (+modifies acromial_r) | **none** |
| ConstraintSet | (none) | scapulohumeral_rhythm_r | **none** |
| ForceSet | DI1_r, ED2_r, ED3_r, FDP2_r, FDP3_r, FDS2_r, FDS3_r, LUM1_r | 11 shoulder-girdle muscles | **none** |

**Verdict: no collision.** The two forks touch disjoint ends of the same kinematic chain
(shoulder-girdle vs wrist/hand); neither touches `humerus_r`/`ulna_r`/`radius_r` or the elbow
joint. `merge_unified_v2.py`'s `build_unified_v2()` prints this table live as
`arm_chain_precheck` and it is also in the evidence JSON. This was re-confirmed at the gate that
actually matters (SS4-5 below: `initSystem()`, forward dynamics, and marker IK on the fully
composed model), not trusted from the XML diff alone.

## 3. What did NOT merge cleanly on the first attempt — forced, diagnosed, fixed

**This is the one real finding of this task**, and it is reported here in full per the "kills
are auditable" / symmetric-QC discipline, even though it was FIXED before the final artifact was
produced.

None of the original 6 forks ever removed a named component relative to its own parent (verified:
zero `removed` entries across all 6 forks in `merge_unified_model_evidence.json`). Consequently
`merge_unified_model.py`'s own `apply_delta()` was never exercised on a removal — it only ever
applies `added`+`modified` names; a `removed` name is only ever printed as a warning, never
deleted. The hand fork DOES remove one component (`radius_hand_r`).

**Forced test:** a naive reuse of the original apply-only semantics (no explicit removal
handling) was built and run first. Result: `hand_r` ends up with **two** parent joints — the
stale `radius_hand_r` weld (parent `radius_r`) AND the new `wrist_dev_r_joint` (parent
`wrist_int_r`). `initSystem()` does **not** raise a Python exception. OpenSim's own log instead
prints, and Python-level code sees nothing:
```
[error] Model unable to assemble: ... Optimizer failed: Ipopt: Infeasible problem detected
Assembly error tolerance achieved: 0.34202 required: 1e-10. Model relaxing constraints and trying again.
```
...then **silently returns a "successful" state anyway** via constraint relaxation. This is a
real false-pass that a naive `try/except` around `initSystem()` would never catch — exactly the
kind of adversary the watertight discipline exists to force before trusting a green light.

**Fix:** `merge_unified_v2.py` adds a PRE-REGISTERED `expected_removals` map per fork
(`{"full_hand": {"JointSet": {"radius_hand_r"}}}`) and runs a real deletion pass before the
add/modify pass; any OTHER, non-pre-registered removal is left in place and loudly flagged in
the evidence JSON (`unexpected_removals`), never silently dropped.

**General, machine-checkable gate added** (not log-text eyeballing, per the watertight
discipline's own "never eyeball a figure/log" rule): `check_no_duplicate_parent_bodies(model)`
walks every `Joint`'s child frame to its base (`Frame.findBaseFrame()`) and asserts each body is
claimed by at most one joint — a direct structural/graph-theoretic invariant (in-degree <= 1 per
body in the kinematic tree), not a proxy. Confirmed empirically: flags `{'hand_r': 2}` on the
naive/broken graft, `{}` (clean) after the fix, `{}` on the actual merged v2 model.

## 4. Verification

### 4a. Structural counts — pre-registered "sum of parts", measured on the real `initSystem()`-loaded model

| | measured | expected (current-unified + full_hand delta + hip_ankle delta) | |
|---|---:|---:|---|
| bodies | 34 | 34 | PASS |
| joints | 34 | 34 | PASS |
| muscles | 240 | 240 | PASS |
| ligaments | 111 | 111 | PASS |
| reserve actuators | 13 | 13 | PASS |
| constraints | 3 | 3 | PASS |
| markers | 49 | 49 | PASS |
| **total force elements** | **364** | **364** | **PASS** |

Total mass: 79.0341 kg (vs 78.9759 kg pre-merge). The +0.0582 kg delta is fully and exactly
accounted for by the 7 new hand bodies (0.05816 kg summed directly from the model — index/middle
phalanx cylinders + a near-massless 1e-5 kg `wrist_int_r` dummy frame body) — no unaccounted
mass appeared or disappeared. Structural invariant (SS3): **PASS**, zero duplicate-parented
bodies.

### 4b. Smoke test A — forward dynamics (5 ms)

`equilibrateMuscles()` + `Manager.integrate(0.005)` completes; all 8 new hand/wrist coordinates
finite after the step (e.g. `pip3_flex_r` = 0.096 rad, `mcp2_flex_r` = -0.015 rad — small,
sane, un-actuated drift under gravity/inertia over 5 ms, not exploding). `any_nonfinite=False`.

### 4c. Smoke test B — real marker IK on `walking1.trc` vs the external reference `walking1.mot`

Truly-unaffected coordinates (pelvis/hip/knee/lumbar/left-arm) reproduce the reference closely;
worst is `lumbar_bending` at 1.058 deg RMSE — inside the precedent-anchored 2.0 deg ceiling
(`GENEROUS_REASSEMBLY_TOL_DEG`, same convention as the original merge) and near-identical to the
**original 6-fork merge's own value for the same coordinate (1.050 deg)** — i.e. adding the hand
+ hip/ankle-ligament forks did not perturb this at all beyond ordinary IK re-solve noise.

The two ALREADY-KNOWN, pre-existing issues are carried forward, confirmed **not worsened**:

| coordinate | original 6-fork merge (RMSE / maxabs, deg) | this v2 merge (RMSE / maxabs, deg) |
|---|---|---|
| ankle_angle_r (foot midtarsal graft) | 1.441 / 3.056 | 1.452 / 2.955 |
| subtalar_angle_r | 1.108 / 3.183 | 1.162 / 2.550 |
| arm_rot_r (shoulder registration) | 55.16 / 56.59 | 52.36 / 55.88 |
| arm_add_r | 11.93 / 14.56 | 12.08 / 14.66 |
| elbow_flex_r | 8.74 / 10.27 | 8.19 / 9.45 |

Same regime, small (a few percent) run-to-run differences consistent with ordinary IK-solver
noise from a slightly larger coordinate set — not a new defect introduced by this merge.

**New hand/wrist/finger coordinates — honest disclosure, not a hidden gap:** `walking1.trc` has
no hand/finger markers (the model's marker count is unchanged at 49; neither new fork adds any).
These 8 coordinates are therefore **unweighted** in the IK objective — their solved trajectories
(e.g. `mcp2_flex_r` sits at a stable ~30.1 deg the entire trial, `wrist_flex_r` at exactly 0.0
deg throughout) are an artifact of wherever the underdetermined assembly optimizer's solution
landed, **not physically meaningful hand posture** — confirmed this is not a repeat of the SS3
assembly pathology (zero `[error]`/`Assembly error`/`relaxing constraints` lines anywhere in this
run's full log, checked by direct grep, not assumed). Recovering meaningful hand kinematics from
this walking trial would require actual hand markers, which this mocap dataset does not have.

### 4d. Re-verification: do PRIOR certified findings survive composition?

Reused each build's OWN proven verification code (not reinvented) directly on the composed v2
model/state:

- **Knee** (`reverify_knee_mismatch`, reused unmodified from `merge_unified_model.py`): max MCL
  strain at the knee's range ceiling = **58.2%**, reproducing the already-certified finding
  exactly.
- **Hip/ankle ligaments** (`sweep_1d`/`band_mean` reused from `add_hip_ankle_ligaments.py`, sign
  conventions RE-MEASURED live on the merged model, not trusted from the cached config):
  ILFL tighter-in-extension 300.3N vs 17.1N (PASS, matches the standalone build's own 300N/17N);
  ATFL tighter-in-plantarflexion 10.5N vs 0.0N (PASS); PTFL tighter-in-dorsiflexion 174.8N vs
  88.0N (PASS); tension-only check 0/1397 violations across both sweeps.
- **Hand moment arms** (OpenSim's own `computeMomentArm()`, matched to the standalone build's
  EXACT test pose after a forced correction — see below): LUM1_r at mcp2_flex_r = 8.8911mm
  (anchor 8.8911mm, abs err 2.3e-10mm) and at pip2_flex_r = -2.2570mm (anchor -2.2570mm, abs err
  3.4e-9mm) — the textbook lumbrical double-action (MCP flexion + PIP extension from the same
  muscle) reproduces to near machine precision. DI1_r at mcp2_flex_r = 6.3194mm (anchor
  6.3194mm) also reproduces.

  **Forced OODA, disclosed rather than hidden:** the first attempt at this specific check
  measured moment arms at the all-zero coordinate pose and got 8.146mm vs the 8.891mm anchor (an
  8.4% gap that would have read as a merge defect). Moment arms are pose-dependent by
  construction (via-point tendon geometry), so this was an apples-to-oranges comparison, not a
  finding. Orient: re-read `full_hand.py`'s own `verify()` code and found the anchor was measured
  at a specific flexed test pose (`mcp2/pip2/dip2 = 20/20/10 deg`), not zero. Reproducing that
  exact pose is what produced the near-machine-precision agreement above. Reported here per the
  "honest-negative is not a free pass — force the fix before declaring one" rule: a one-shot
  "FAIL" would have been a false negative from a self-inflicted test-construction bug, not a real
  finding about the merge.

## 5. Forks that remain separate — not touched by this task, and why

Scoped explicitly to the 2 forks the task named. Other standalone `.osim` files in this repo are
**not** anatomical-content forks of subject2's skeleton in the same sense and were left alone:
`subject2_scaled_pcsa_fmax_proxyA_mass.osim` / `...proxyB_geom.osim` (alternate Fmax
*parameterizations* of the same 80 native muscles, both already self-disclosed as
non-improving experiments, not new anatomy); `LaiArnoldModified2017_mediapipe_bridge_...osim`
(a different pose/marker-bridging target, not a muscle/ligament content fork);
`LaiArnold_hallux_proto_scaled.osim` (a separate big-toe DOF prototype, never selected the way
`LaiArnold_midtarsal_proto_scaled.osim` was for the already-merged foot fork); the JAM/COMAK
`DM.osim` implant model and `hip_cartilage_contact/model/synthetic_hip.osim` (both excluded by
construction — different individuals/lineages entirely, not subject2, per the original merge's
own precedent); and the various `band_*`/`elastic_band`/`subject2_spine_stoop_lift` scenario
probes (transient load-scenario scaffolding, not persistent anatomy). None of these were in the
task's scope.

## 6. Reproduction

```
source_repository/.venv-msk/bin/python3 scripts/msk/merge_unified_v2.py
```
Regenerates any missing input from its own script (spine ligaments / full hand / hip-ankle
ligaments), grafts, writes `data/msk_models/subject2_unified_v2.osim`, runs all verification
steps above, and writes `scripts/msk/merge_unified_v2_evidence.json`. No git commit, no git push
(isolation-respected, per `COORDINATOR.md` SS1 and the task).

## Evidence index

`scripts/msk/merge_unified_v2.py` (this merge), `scripts/msk/merge_unified_v2_evidence.json`
(full machine-readable output: counts, arm-chain precheck, conflicts/unexpected-removals lists,
both smoke tests, all three reverification blocks). Cross-checked directly against
`scripts/msk/merge_unified_model_evidence.json` (original 6-fork merge, for the shoulder/foot
IK comparison in SS4c) and `scripts/msk/full_hand_evidence.json` (source of the exact moment-arm
anchors in SS4d). Read directly: `scripts/msk/full_hand.py`, `scripts/msk/anatomical_hand.py`,
`scripts/msk/add_hip_ankle_ligaments.py`, `scripts/msk/add_spine_ligaments.py`,
`scripts/msk/merge_unified_model.py`. IK inputs/outputs:
`data/msk_smoketest/subject2_unified_v2_walking1/`.
