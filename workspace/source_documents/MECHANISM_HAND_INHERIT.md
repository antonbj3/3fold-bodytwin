# MECHANISM HAND INHERIT — CS's HUM-HAND graft assessment (2026-07-21)

Executes the mission: inherit CS's (cad-to-simulation, `upstream/feat/breakthrough-hyperreal`) last-12h
HAND model + video-retarget work into bodytwin, per COORDINATOR.md §2 (bt-copy/re-point, never write back).

**Sharpened evaluation criterion (coordinator reframe, applied throughout this doc):** bodytwin's hand
goal is an **anatomically correct, MUSCLE-DRIVEN hand** — grip must EMERGE from real forearm muscles
(extrinsic flexors/extensors: FDP/FDS/EDC and friends) whose tendons cross the wrist to the finger
joints. That build is owned by a separate agent with a proper anatomical (OpenSim muscle-tendon) donor.
This doc's job is the honest CS-contribution assessment against THAT bar — not to force HUM-HAND into
a role it cannot fill.

## 0. Bottom line

**HUM-HAND is a kinematic/manipulation-robotics hand, not a muscle-driven anatomical hand, and cannot be
the donor for bodytwin's wrist-unwelding/muscle build.** Zero muscle-tendon representation exists
anywhere in its 11-file thread (machine-verified, see §2). It transfers to the anatomical-hand goal in
exactly **one** clean, actuation-agnostic way: its MediaPipe-hand-landmark → joint-angle **video-retarget
pipeline** (a video→kinematics capture layer, independent of what drives the hand). That piece has been
ported into bodytwin, re-verified end-to-end, and run on bodytwin's own corpus footage (§5) — producing a
genuine 31-frame/21-DOF hand trajectory. Everything else (the Newton rigid-body model itself, its PD-servo
grip control, its tool-grasp/manipulation-sequencing chain) is CS-stack-specific and does not drop in.

## 1. What HUM-HAND actually is (format / DoF / topology — measured, not narrated)

Source: `docs/HUM_HAND_V0.md` + `scripts/hand_leg/build_hum_hand_v0.py`, CS commit `973e64836`
(2026-07-18) — **already bt-copied into bodytwin** at that exact commit back on 2026-07-19
(`scripts/hand_leg/build_hum_hand_v0.py`, byte-diff-verified identical to CS's copy, 2026-07-21). This
part of the mission's premise predates this task; the new material below is what CS added in the last
~12h that bodytwin does not yet have.

- **Engine/format:** `newton.ModelBuilder` — a Python-API GPU rigid-body physics engine (NVIDIA
  Newton/Warp, MuJoCo-solver backend), **not OpenSim, not SMPL-X mesh deformation.** Bodies are
  capsule-collision rigid links; joints are `add_joint_revolute` with `target_ke`/`target_kd` (PD
  position-servo control), not muscle-tendon `PathActuator`s.
- **DoF:** 21/hand, the Cobos et al. 2008 ("Efficient Human Hand Kinematics for Manipulation Tasks")
  robotics/prosthetics convention: 4 long fingers × 4 DoF (MCP flex+abd, PIP, DIP) + thumb × 5 DoF
  (CMC flex+abd+opposition, MCP, IP) = 21. This is explicitly **not** SMPL-X's own hand parametrization
  (measured from the real `SMPLX_NEUTRAL_2020.npz`: 15 joints/hand, 45 ball-joint DoF/hand) — HUM-HAND
  *implements* the Cobos convention as revolute joints, merely **anchoring** joint positions and bone
  lengths in real SMPL-X T-pose coordinates (`J_regressor @ v_template`).
  - Declared simplification (CS's own doc, not my inference): a **straight-chain topology per finger**
    with global axis conventions (flexion=local Z, abduction=local Y, thumb opposition=local X) —
    "one simplified direction of motion per joint," not the full 3D bone orientation SMPL-X's own LBS skinning
    would give.
- **No wrist DoF at all.** The wrist is welded to world (`add_joint_fixed`) in every HUM-HAND build —
  the Cobos convention itself excludes the wrist, and CS's build never adds one back.
- **Geometry:** real (SMPL-X-measured, not guessed) per-phalanx segment lengths; capsule collision
  shapes; friction μ=0.62 sourced externally (Han et al., finger-vs-acrylic).
- **What it's used for downstream** (the 8 commits since bodytwin's last sync — full list in §6): video
  hand-pose retargeting, proportion/per-segment calibration of that retarget, tool-grasp pose synthesis
  (screwdriver/spatula), a manipulation sequence planner, contact-sim of tool-on-workpiece, a sim→render
  pose bridge, and a first deformable-material (film) manipulation cell. All of it is **grip-pose /
  manipulation-planning work on a PD-servo-actuated rigid-body hand** — a robotics/animation regime, not
  a physiology regime.

## 2. Machine-verified: zero muscle-tendon representation anywhere in this thread

Per rule 1 (force the adversary, don't assume): grepped `muscle|tendon|ligament|millard|thelen|
coordinateactuator|pathactuator` (case-insensitive) across **all 11** CS `scripts/hand_leg/*.py` files
at `upstream/feat/breakthrough-hyperreal` HEAD — `build_hum_hand_v0.py`, `grasp_impedance_v2.py`,
`grasp_soft_contact_v3.py`, `grasp_synthesis_v1.py`, `hum_contact_sim_cell.py`, `hum_toolgrasp_cell.py`,
`hum_hand_video_retarget_cell.py`, `hum_hand_retarget_calib_cell.py`, `hum_hand_persegment_topology_cell.py`,
`hum_deformable_cell.py`, `hum_render_bridge_cell.py`. **Zero hits, every file** (including the 3
grasp-control experiments — impedance/soft-contact/synthesis — which are the files most likely to smuggle
in a tendon-routing model if one existed; they don't, they're joint-space PD/impedance control on the same
rigid capsules). This is a machine cross-check, not an inference from having read one file.

**Consequence for the anatomical-muscle-driven goal:** there is no FDP/FDS/EDC/EIP/EPL/FPL/APL/lumbrical/
interosseous representation, no tendon path, no wrapping surface, no pulley (A1–A5/palmar aponeurosis),
and — because the wrist itself is welded shut — **no joint for a wrist-crossing tendon to have a moment
arm across in the first place.** HUM-HAND cannot be adapted into a muscle-driven hand by adding muscles
to it; the entire actuation substrate (PD-servo joint targets) is architecturally the thing bodytwin's
goal replaces. The real donor bodytwin needs is a genuine OpenSim musculoskeletal hand/forearm model with
`Millard2012EquilibriumMuscle`/`Thelen2003Muscle` `PathActuator`s for the extrinsic flexors/extensors and
an articulated (non-welded) wrist for their tendons to cross — a different, separately-sourced model
(the kind of donor the Holzbaur/Saul-lineage OpenSim upper-limb models represent), which is out of this
task's scope (owned by the other agent).

## 3. What DOES transfer cleanly — and why it's honest, not forced

**(a) The MediaPipe HandLandmarker video-retarget pipeline — actuation-agnostic, transfers in full.**
CS commit `7a5011a74` (`hum_hand_video_retarget_cell.py` +
`extract_hum_hand_video_retarget_landmarks.py`): 21-landmark MediaPipe hand tracking → algebraic inverse
of the forward-kinematics chain → per-frame joint angles in the same Cobos-anatomical naming (MCP
flex/abd, PIP, DIP, CMC flex/abd/opposition). This is **pure video-landmark geometry** — it has nothing
to do with what actuates the resulting pose. It transfers to *any* hand model sharing the same
anatomical joint taxonomy, muscle-driven or not, because MCP/PIP/DIP/CMC are real anatomical joints, not
HUM-HAND-specific inventions. **Value for the anatomical-muscle build:** a video→target-pose capture
layer, usable as an IK-tracking target or ROM/plausibility cross-check for whichever OpenSim
muscle-driven hand gets built — a validation/reference-data contribution, not a structural one.

**(b) The DoF/ROM taxonomy itself — a minor, optional cross-check, not a structural donor.** The Cobos
21-DoF convention is a legitimate anatomical simplification (real MCP/PIP/DIP/CMC joints, real ROM
limits used in the revolute-joint bounds) — usable only as a naming/ROM sanity cross-check against
whatever coordinate set the anatomical model exposes, nothing more. Its own straight-chain topology is a
declared simplification (CS's own doc: loses true 3D bone orientation / ab-adduction coupling), so it is
not worth encoding as a structural reference beyond "these are the joints and roughly this ROM."

**(c) Not worth porting:** the SMPL-X T-pose segment lengths (a generic neutral-body-mesh estimate;
redundant once a real anatomical OpenSim hand model brings its own validated segment/cadaver-based
geometry) — low marginal value, correctly left out.

## 4. Scope already inherited vs genuinely new (verified before touching anything)

Bodytwin's `scripts/hand_leg/` already contains CS's HUM-HAND v0 + 3 grasp-control experiments
(build_hum_hand_v0.py, grasp_impedance_v2/soft_contact_v3/synthesis_v1, extract_hand_landmarks.py,
hand_ruler_calibration.py, hand_landmarker.task), merged in on 2026-07-19 up through CS commit
`09c3addde` — confirmed both a byte-identical diff on `build_hum_hand_v0.py` and
`git merge-base --is-ancestor 09c3addde HEAD` = true (a real merge, not a loose hand-copy). **The
genuinely new material this task covers is exactly 8 commits since then:**

| commit | cell | in bodytwin before this task? |
|---|---|---|
| `7a5011a74` | HUM-HAND-VIDEO-RETARGET (v0) | no — **ported below** |
| `b7942377e` | HUM-HAND-RETARGET-CALIB (per-finger scale) | no — assessed, not ported (see §7) |
| `6bbe528b5` | HUM-HAND-PERSEGMENT-TOPOLOGY (per-segment scale) | no — assessed, not ported (see §7) |
| `b0b0afac2` | HUM-TOOLGRASP | no — out of scope (tool manipulation, not hand anatomy) |
| `5c41e4a37` | HUM-SEQ-PLANNER | no — out of scope |
| `28a1c17c9` | HUM-CONTACT-SIM | no — out of scope |
| `52fe96671` | HUM-RENDER-BRIDGE | no — out of scope |
| `91fffd87c` | HUM-DEFORMABLE | no — out of scope |

## 5. First concrete step completed: the retarget pipeline runs in bodytwin's stack, on bodytwin's own footage

**Environment fact-check first** (not assumed): `.venv-msk` has OpenSim 4.6 but no MediaPipe;
`.venv-humancap` does not exist in bodytwin; `~/miniconda3/bin/python3` has MediaPipe 0.10.32
(same interpreter CS uses). The retarget cell itself needs only numpy (FingerChain/load_smplx_hand_geometry
import nothing else at module level) — runs under plain system `python3`.

**Clip selection — measured, not eyeballed** (rule 4, SCENE-EYES: MEASURE, never "look" at a frame).
Wrote `scripts/msk/inherited_hand/_sweep_clip_windows.py` and ran a real HandLandmarker detection-rate
sweep, mirroring CS's own field-verification discipline (they field-rejected a gloved watch-teardown
clip at 0/239 before picking their loom-weaving clip): 9 candidate windows across 4
`@samuelwatson__` speed-climbing clips —

| clip | window | detection_rate |
|---|---|---|
| DM0hFD_J4TC | start/mid/finish | 0.0 / 0.0 / 0.0 |
| DNignJfMDF3_4 | start/finish | 0.0 / 0.0 |
| DU83RptkqFq | start/finish | 0.435 / 0.0 |
| **DNrDGcEXGtX** | start | 0.130 |
| **DNrDGcEXGtX** | **finish** | **0.955** |

7/9 windows: 0.0 detection — broadcast-wide speed-climb framing puts the hand too small/fast for
MediaPipe Hands (honest negative, mechanism analogous to CS's own gloved-macro reject — different cause,
same "too hard for this detector" shape). Refined `DNrDGcEXGtX` to 1s resolution over its full 8.45s
duration: `[0.25, 0.0, 0.125, 0.0, 0.625, 1.0, 1.0, 1.0]` — a genuine monotonic rise into a finish/
celebration beat at t=4.0–8.45s, not a single lucky bin. Locked the extraction window there.

**Extraction (real run, not simulated):** `scripts/msk/inherited_hand/extract_hand_video_retarget_landmarks.py`
→ 44 sampled frames (10fps over 4.45s), **38/44 (86.4%) any-hand detection**, cached to
`data/cache/inherited_hand/DNrDGcEXGtX_t4_8p45_hands.json`.

**Retarget cell (real run):** `scripts/msk/inherited_hand/hand_video_retarget_cell.py` → **31 valid
frames × 21-DOF trajectory**, persisted in full to `reports/probes/inherited_hand_video_retarget.json`
(this is the actual, consumable "hand trajectory" the mission asked for — CS's own cell computes this
per-frame array in memory but never writes it to disk; that gap is closed here, declared as bodytwin's
addition in the file's own docstring).

**Every function that does the actual video→angle math (`extract_frame_angles` + its 6 helpers,
`run_synthetic_control`) was ported line-for-line and diffed against the CS original** — a raw `diff`
of the two files' shared function bodies shows **zero code-line differences**, only docstring trims
(machine-checked, not eyeballed; caught and fixed one of my own transcription slips — an unnecessary
"cleanup" rewrite of `_finger_landmarks_forward` — before running, by re-diffing against source).

### Results — honest, side-by-side against CS's own clip (independent generalization test)

| metric | CS (loom-weaving, static hands) | bodytwin (climbing grip, dynamic) |
|---|---|---|
| n_valid_frames | 43 | 31 |
| kinematic_retarget_rms (band <0.05m) | 0.0689m — **misses** | 0.0749m — **misses** |
| per-finger worst | Thumb 0.0970m | Thumb 0.1023m |
| shuffle-null gate (≤1/12 beats) | 0/12 — **PASS** | 0/12 — **PASS** |
| synthetic-control round-trip | exact (~3e-15 rad) | exact (**3.16e-15 rad**, independently re-run) |
| grasp-aperture closing ratio | 9.25× (0.012→0.111m) | 7.27× (0.026→0.189m) |
| raw MediaPipe hand span vs anthro band [0.15,0.21]m | 0.105m — **out of band** | 0.084m — **out of band** |
| verdict | partial | **partial** |

This is a genuine **over-determination / generalization result, not a re-run of the same thing**: the
same qualitative pattern (RMS band-miss from a real cross-model SMPL-X-neutral-vs-filmed-hand proportion
mismatch, alongside a strong real temporal/grasp signal) reproduces on a fully independent clip, athlete,
and content domain (fast dynamic climbing grip vs static/slow loom-weaving manipulation). That
strengthens confidence the RMS-miss is a genuine, general phenomenon (not a fluke of CS's one clip) —
**and confirms the pipeline itself (geometry + temporal-coherence detection) generalizes cleanly, which
is the actually load-bearing claim for reuse as validation/reference data.** It does *not* mean the
absolute-scale retarget problem is solved — it isn't, on either clip.

**One functional caveat worth flagging on its own merits:** the single DoF the pipeline cannot observe
and fixes at 0 — thumb CMC opposition — is arguably the most functionally important DoF for grip itself
(opposition is what lets the thumb meet the fingers at all). The retarget's 20/21-DoF coverage is
real, but the missing 1/21 is not a random omission.

## 6. Files delivered (all new, bodytwin-only, provenance-tagged)

- `scripts/msk/inherited_hand/extract_hand_video_retarget_landmarks.py` — adapted from CS
  `7a5011a74e12badd3fbeca622cd4905290410404:scripts/hand_leg/extract_hum_hand_video_retarget_landmarks.py`
  (repointed at bodytwin's own clip/window; extraction algorithm unchanged).
- `scripts/msk/inherited_hand/hand_video_retarget_cell.py` — adapted from CS
  `7a5011a74e12badd3fbeca622cd4905290410404:scripts/hand_leg/hum_hand_video_retarget_cell.py`
  (core math ported verbatim, diff-verified; adds full per-frame trajectory persistence).
- `scripts/msk/inherited_hand/_sweep_clip_windows.py` — new, bodytwin-only scouting tool (not a CS port).
- `data/cache/inherited_hand/DNrDGcEXGtX_t4_8p45_hands.json` — raw MediaPipe landmark cache (durable, gitignored-scale).
- `reports/probes/inherited_hand_video_retarget.json` — the report + the 31-frame/21-DOF trajectory + ATOMS-style machine-checkable claims.
- This doc.

No file outside these paths was touched. No CS or upstream file was written. No commit was made (per
COORDINATOR.md §1, the coordinator commits).

## 7. Honest graft plan for the welded hand (what's NOT done, and why that's correct)

**Not done, deliberately:** un-welding `radius_hand_r`/`radius_hand_l` in the LaiArnold model, or adding
any finger/wrist DoF to bodytwin's OpenSim tree. Per the coordinator's reframe, that build belongs to a
separate agent working from a genuine anatomical (muscle-tendon) donor — HUM-HAND is not that donor
(§2), so doing this myself risked (a) forcing a robotics DoF convention into an anatomical model where it
doesn't belong, and (b) colliding with the other agent's in-progress work on the same model files (none
of which this task touched — verified via `git status` before and after).

**What the other agent's donor search should look for, based on what's confirmed missing here:**
an OpenSim upper-limb/hand model with real `PathActuator`-class muscles for the extrinsic finger
flexors/extensors (FDS, FDP, EDC at minimum) whose paths cross a **non-welded, articulated** wrist
(flexion/extension + radial/ulnar deviation, i.e. `radius_hand` needs to become a real 2-DoF joint, not
stay a `WeldJoint`), with wrapping surfaces/pulleys at the wrist and MCP level for realistic moment arms.
That is a different model family (e.g. the Holzbaur/Saul-lineage OpenSim dynamic-arm models carry this;
HUM-HAND categorically does not and cannot be extended into it without discarding its entire actuation
layer).

**What this task's output is useful for, once that donor lands:** the 21-DoF trajectory format + the
now-bodytwin-native extraction/retarget pipeline (§5) is ready to (1) re-target onto the anatomical
model's own coordinate names (a joint-name mapping, not a re-derivation, since both use the same
MCP/PIP/DIP/CMC anatomical taxonomy) as IK-tracking targets, and (2) supply a real, non-synthetic ROM/
grasp-aperture cross-check (this run alone shows a real 7.27× aperture-closing event during an actual
climbing grip) for validating that the muscle-driven hand can reach the poses real hands reach.
