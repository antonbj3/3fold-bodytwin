# MECHANISM MSK BUILD PLAN — verified bottom-up plan for the loadable biomechanical twin

Written 2026-07-21 by the plan agent. Every path/version/import claim below was checked LIVE this session
(`ls`/`find`/`grep`/direct Python import/`opensim-cmd --version`/C++ header grep) — nothing here is recalled
from training-data memory of "what OpenSim probably has." Where a claim could NOT be verified this session it is
marked so explicitly in §6. This document does not modify `data/MECHANISM_ANCHOR_GRAPH.json` or any other file.

---

## 1. Capability Inventory (verified)

### 1.1 bodytwin repo itself
- **No venv** anywhere under `~/projects/bodytwin` (checked to depth 4). System `python3` = 3.10.12,
  `numpy` 2.2.6 — this is the shared interpreter; do not `pip install` into it (matches the standing operating
  fact "no bodytwin venv → don't pip-install into the shared env").
- **No `.osim` file** anywhere in the repo (`find -iname '*.osim'` empty).
- `uv` 0.5.9 present at `~/.local/bin/uv`. The **opensim 4.6 cp312 wheel is already cached**:
  `~/.cache/uv/wheels-v3/pypi/opensim/opensim-4.6-cp312-cp312-manylinux_2_27_x86_64.manylinux_2_28_x86_64` — a
  fresh `uv venv --python 3.12 && uv pip install opensim` needs **no network fetch** for opensim itself.

### 1.2 cad-to-simulation-I venvs — verified by actually importing, not just `ls`
| venv | python (verified) | import check (verified live) |
|---|---|---|
| `.venv-humancap` | 3.13.11 | `import mediapipe` → `0.10.32` OK |
| `.venv-newton` | 3.10.12 | `import torch` → `2.11.0+cu128`, `torch.cuda.is_available()==True` |
| `.venv-geoprobe` | — | present, not opensim-related, out of scope here |
| `.venv-lora` | — | present, not opensim-related, out of scope here |

None of the 4 import `opensim` (`ModuleNotFoundError` in all four, verified). MediaPipe model file is exactly
where the handoff doc claims: `/media/anton/183E48713E484A48/cad-to-simulation_video/geoprobe_scratch/weights/pose_landmarker_full.task`
(9.4 MB, present). **The handoff's pose-pass prerequisite claim HOLDS**, verified independently.

### 1.3 OpenSim — THREE independent installations found, none inside a bodytwin/cad-to-simulation-I venv
**A. opensim-core-jam C++ fork build — working, verified live this session.**
`/media/anton/8838D60F38D5FBDE/mechanism_data/opensim_jam_build/opensim-build/opensim-cmd` — ran `--version` live:
`OpenSim version 4.5.2-2026-01-19-a3c872a2b`; `ldd` shows zero missing libraries. Has COMAKTool / ForsimTool /
COMAKInverseKinematicsTool / JointMechanicsTool (JAM's contact-mechanics tools), invoked via
`opensim-cmd run-tool <settings.xml>` — pure C++, no Python bindings, no SWIG. Full source checkout
`opensim-core-jam/` sits alongside (~230 example `.osim` files). Footprint 8.6 GB. **This built the KNEE-CELL
flagship result** (bt_memory-recorded, artifacts still on disk under `run_DM_ngait_og{1,3,3b,4}/results/`):
predicted tibiofemoral contact 2.503 BW vs measured 2.584 BW peak (≈3% error, r=0.802).

**B. video2kin Python pipeline (Pose2Sim + opensim 4.5.2 bindings) — built 2026-07-18, BROKEN IN PLACE (verified + root-caused this session).**
`/media/anton/8838D60F38D5FBDE/mechanism_data/video2kin/{envs/opensim_env, venv_pose2sim, micromamba_bin}`.
Calling `bin/python3` fails `Exec format error`. Diagnosis: `bin/python3` is a 28–40-byte **"IntxLNK" NTFS
reparse-point placeholder** (the marker NTFS uses for a Linux symlink) that the in-kernel `ntfs3` driver mounting
this drive is **not** resolving to the real target — `cat` on the "file" shows literal `IntxLNK` + a UTF-16
target string, not a valid ELF/shebang. Bypassing to the real `bin/python3.10` (verified genuine ELF, 17 MB) gets
past that but then fails on `ImportError: .../libstdc++.so.6: file too short` — the *same* symlink issue
recurring on an internal shared-library symlink inside the conda env. `venv_pose2sim`'s own `bin/python3.10` is
itself a broken placeholder (worse off than the opensim env).
**Consequence: the previously-closed "video→knee-force chain" (mean 8.78% peak error, r=0.94 vs mocap-GT) is
currently non-reproducible in place** — not a science regression, a filesystem one. Fix in §5.
**Generalizable finding (new, worth keeping): neither external drive (`183E48713E484A48` nor `8838D60F38D5FBDE`,
both NTFS) is safe to host a Python venv/conda env** — fine for bulk data (models/csv/video), unsafe for anything
with internal symlinks (venvs, conda envs). Rebuild any such toolchain on ext4.

**C. The cached-wheel route (§1.1)** is the clean way to get opensim 4.6 Python bindings again, on ext4 — this
is the one to use going forward, not (B).

### 1.4 `.osim` files — exhaustive search result
Zero anywhere under `~` (excl. data drives) or `/media/anton/183E48713E484A48`. All real `.osim` files
live on `/media/anton/8838D60F38D5FBDE/mechanism_data/` (df: 931G total, 198G free, ntfs3):
- `opensim_jam_build/` — ~140 example/test `.osim` (incl. `Rajagopal2015_right_leg_9musc.osim`,
  `ThoracoscapularShoulderModel.osim`, `knee_patella_ligament.osim`, `toyLigamentModel.osim`, gait2392/gait2354
  examples) + the flagship's `run_DM_*/results/comak-inverse-kinematics/ik_constrained_model.osim` +
  **`video_bridge/`** (below). 8.6 GB.
- `video2kin/` — `Rajagopal_2015.osim` (bundled opensim-doc example) + Pose2Sim's own
  `Model_Pose2Sim_{simple,muscles_flex}.osim` + 3 already-produced scaled outputs from real prior runs
  (`run_single_*.osim`, `run_subject2_walking1_*.osim`). 4.7 GB.
- **`LabValidation_withVideos/` — 10 subjects (subject2–11, verified count), 23 GB.** Each has
  `LaiArnoldModified2017_poly_withArms_weldHand_{generic,scaled}.osim` from Mocap AND from every
  Video/{HRNet, OpenPose_default, OpenPose_highAccuracy}×{2,3,5-camera} combination (OpenCap dataset).

`opensim_jam_build/video_bridge/`: `DM_video_mapped_full.mot`, `DM_video_mapped_window.mot`,
`comak_settings_video.xml`, `results/comak/` — the **already-proven** coordinate-name/sign bridge from a
Pose2Sim `.mot` into the JAM "DM" knee model (float-exact per `bt_memory/video-to-kinematics-frontend-pose2sim-offline-recipe.md`),
blocked only on a video clip with synchronized GRF for COMAK to converge on.

`knee_gc_measured/`: **47** per-trial measured-knee-force CSVs (Grand Challenge "DM" subject:
`Fx,Fy,Fz,Tx,Ty,Tz,GON,GRFz`), already extracted — the flagship's raw in-vivo anchor, covering
gait/stairs/squat/chair-rise/bouncy/crouch/mtpgait/static/one-leg-stand trials.

### 1.5 Existing MSK-* graph nodes — literature cells, NOT simulateable-model components (important distinction)
`data/MECHANISM_ANCHOR_GRAPH.json` (998 nodes, read-only checked) already has `MSK-JOINT-CONTACT-FORCE-ANCHORS`,
`MSK-ANKLE-TALOCRURAL`, `MSK-TENDON-SPRING-SAFETY-FACTOR`, `MSK-HIP-FEMOROACETABULAR`, `MSK-SHOULDER-GLENOHUMERAL`,
`KNEE-CELL` (all `status=OPEN` except `KNEE-CELL`=`ASSUMED`). Their `claim` fields show these are
**population-epidemiology literature cells** (e.g. `MSK-HIP-FEMOROACETABULAR` = cam-morphology-alpha-angle vs
hip-OA cohort studies; `MSK-SHOULDER-GLENOHUMERAL` = instability-recurrence cohorts) — real and useful, but a
different artifact class than a loadable OpenSim component. Couple the new simulateable build INTO these (e.g. a
built hip model's joint-contact-stress output as a mechanistic leg for `MSK-HIP-FEMOROACETABULAR`), don't confuse
the two. **`MSK-JOINT-CONTACT-FORCE-ANCHORS` is the one directly on point** — it already contains a fully-sourced
per-joint in-vivo-anchor table; §4 below is built on it (see `data/body_twin/agent_outputs/auto__a89669e49d76e0bfb.json`
for the full text), not re-derived.

### 1.6 ig_downloads corpus — verified band-content check
`~/ig_downloads/ig_athletics__semantic.jsonl` (893 rows) — grep count for `"band"` = **0**. The separate 5-clip
`athlete_video_corpus` (strongman deadlifts / UFC weigh-ins / DEXA) is also 0 band-relevant. **No existing corpus
covers the elastic-band use case** — see §6.

---

## 2. Model Selection

**Verified joint list** (grepped directly from the actual `.osim` XML, both files on disk — not from
memory/docs) for `Rajagopal_2015.osim` (`video2kin/.../OpenSenseExample/Rajagopal_2015.osim`) AND
`LaiArnoldModified2017_poly_withArms_weldHand_generic.osim` (`LabValidation_withVideos/subject2/.../Mocap/Model/`):
identical joint sets in both — `ground_pelvis, hip_{r,l}, walker_knee_{r,l}, patellofemoral_{r,l}, ankle_{r,l},
subtalar_{r,l}, mtp_{r,l}, back, acromial_{r,l}, elbow_{r,l}, radioulnar_{r,l}`. LaiArnold additionally has
`radius_hand_{r,l}` as a **`WeldJoint`** (wrist present as a body but fused, 0 DOF).

**Pick: `LaiArnoldModified2017_poly_withArms_weldHand`** (Lai, Arnold & Wu 2017, muscle-driven, "poly" = polynomial
muscle-geometry paths, `withArms` = full upper body), as already staged per-subject (generic + scaled) in
`LabValidation_withVideos/`. Reasons, all verified not assumed:
1. **Already on disk, already scaled to 10 real subjects, already wired into a proven pipeline** (multi-camera
   markerless video → this exact model, mocap-IK ground truth co-registered) — zero acquisition cost, and it is
   the model the "~4.7° RMSE video→kinematics" and the "video→knee-force 8.78% peak error" results were measured
   against, so reusing it keeps the whole existing chain valid.
2. **Multi-segment foot already present, not lumped**: separate `ankle` (talocrural) + `subtalar` + `mtp` joints
   bilaterally — i.e. talus / calcaneus+foot / toes are already 3 distinct bodies connected by 3 distinct 1-DOF
   pin joints. This **already satisfies the operator's literal ask** ("subtalar + MTP, not a lumped foot") out of
   the box. What it does NOT have (see §6): a multi-axis subtalar (real subtalar has an oblique screw axis; here
   it's a single hinge — the standard simplification used by essentially every gait-analysis-grade OpenSim
   model), a separate midtarsal/Chopart joint (talonavicular+calcaneocuboid — foot mid/forefoot moves as one
   rigid piece from calcaneus to the mtp break), or per-toe articulation (mtp is one lumped hinge for all toes).
3. **Full muscle-tendon actuation + all major joints**: the sibling JAM "DM" gait2392-derived model already
   verified to run in COMAK carries 44 Millard2012 Hill-type muscles + 149 Blankevoort1991 ligaments per leg
   (bt_memory-recorded); Rajagopal/LaiArnold are the same lineage (Rajagopal2016 muscle set extended by Lai/Arnold
   with an arm chain) — muscle-driven, not just torque-driven.
4. **Running validated**: the model is the base for OpenCap's LabValidation set which includes drop-jumps and
   dynamic tasks, not just walking (`bt_memory/video-to-injury-load-two-barriers...md` already ran it on
   drop-jump trials).

**Alternatives evaluated (verified where checkable):**
- **Rajagopal2015/2016 "vanilla"** — same joint set (verified identical list above), same muscle lineage, no arms.
  Present on disk too (`video2kin/.../Rajagopal_2015.osim`, `opensim-core-jam/.../model_Rajagopal2015_posed.osim`).
  Functionally a subset of the LaiArnold pick; use it only if the arm chain is unneeded (lighter/faster sim).
- **Hamner2010 running model** — searched exhaustively (`find -iregex '.*hamner.*'` across the whole
  8838D60F38D5FBDE drive) — **not found anywhere on disk**. Would need external acquisition if ever wanted
  specifically (its main advantage over Rajagopal/LaiArnold — validated at running speeds — is now largely
  redundant since LabValidation already includes dynamic/athletic trials on the LaiArnold model). Not recommended
  as the primary pick given the zero-cost, already-integrated alternative.
- **`ThoracoscapularShoulderModel.osim`** (found in `opensim-core-jam/OpenSim/Tests/shared/`) — a specialty
  scapulothoracic-contact shoulder model shipped with core OpenSim. Flag as a literature-known upgrade path if
  the standard `acromial` ball-joint shoulder proves insufficient (real scapular gliding), not verified beyond
  "the file exists" this session (didn't load/inspect its internals).
- **Detailed foot-ankle model beyond subtalar+mtp** (true midtarsal break, multi-axis subtalar, per-toe): a real
  literature category exists, but **I could not verify a specific name/URL/license this session** — WebSearch
  quota was exhausted mid-session (shared budget) and I will not assert a specific SimTK project number from
  unverified recall (that exact failure mode — trusting a recalled identifier — is a known trap; see honest
  gaps). Treat as an explicit acquisition task, live-verified, before building on it.

**How to obtain:** already on disk (primary route, §1.4). Upstream canonical source (WebFetch-confirmed reachable
this session, generic content only): `github.com/opensim-org/opensim-models` hosts the `.osim` files distributed
with OpenSim (this is where the bundled `video2kin` copies trace back to); LabValidation's per-subject scaled
models trace to the OpenCap SimTK project (`simtk.org/projects/opencap`), already confirmed reachable/used per
prior verified session work.

---

## 3. Elastic-Material Approach

Verified directly from the `opensim-core-jam` C++ source (headers + `.cpp` bodies read, not recalled) — **all
three candidate classes live in OpenSim's CORE libraries** (`OSIMSIMULATION_API` / `OSIMACTUATORS_API`), not the
JAM plugin, so **none of this needs the custom JAM plugin build** — a plain `pip install opensim` (or the
already-cached wheel, §1.1) is sufficient for both tiers below.

### Tier 1 — minimal viable: `OpenSim::PathSpring` (verified class, verified tension-only behavior)
Properties (from `PathSpring.h`): `resting_length` (m), `stiffness` (N/m, linear), `dissipation` (s/m), `path`
(an `AbstractGeometryPath` — the same path-point/wrap-object machinery muscles already use to route around body
segments). Force law, read directly from `PathSpring.cpp`:
```
tension = stiffness * stretch * (1 + dissipation * lengthening_speed)
getStretch(s): length > resting_length ? (length - resting_length) : 0.0
```
This is **exactly** real resistance-band behavior — zero force below the slack length, linear tension above it,
velocity damping — confirmed from the actual `getStretch()` body, not inferred. Modeling a band = one `PathSpring`
per band, `resting_length` = its unstretched routed length, `stiffness` = the band's rated N/m, path points
anchored to a fixed ground frame (rack/floor cleat) at one end and to the loaded body segment (foot, waist, wrist)
at the other, with `WrapCylinder`/`WrapSphere` objects on intervening segments if the band wraps over a limb —
zero new code, same pattern as an existing muscle path.

### Tier 2 — nonlinear J-curve, ALSO zero-new-dependency: repurpose `Blankevoort1991Ligament`
Properties (from `Blankevoort1991Ligament.h`, path `OpenSim/Simulation/Model/` — core, confirmed NOT
JAM-plugin-only despite prior session notes suggesting it needed the plugin): `linear_stiffness` (N, slope past
the transition), `transition_strain` (strain where the curve bends from quadratic to linear, default 6%),
`damping_coefficient` (N·s/strain, default 0.003), `slack_length` (m). This is a quadratic-toe → linear
force-STRAIN curve — a materially better fit to a real elastic band's nonlinear stiffening at large stretch
(bands are typically used at 100–300% elongation, well outside a linear regime) than Tier 1, at the cost of
calibrating 2 shape parameters against a manufacturer's published force-vs-%-elongation table for the specific
band (cheap, one fitting pass). Already compiled and exercised this session's ancestor work (the JAM knee model
uses 149 instances of this exact class for real ligaments) — reusing it for a band is a relabeling, not new
engineering.

### Tier 3 — highest fidelity: FEBio soft-body / hyperelastic membrane
For band self-contact, pressure distribution on skin, or non-uniform stretch across a wide flat band — a true
Mooney-Rivlin/Ogden hyperelastic shell, one/two-way coupled to the OpenSim skeleton at contact patches. **Verified
NOT installed anywhere on this box** (`which febio*` empty, `find -iname '*febio*'` empty to depth 4 from root).
Recommend only if Tier 1/2 prove insufficient in practice (e.g. a visible flattening/pressure-mapping use case) —
it is a real new toolchain acquisition + coupling-engineering effort, not a config change.

**Recommendation: build Tier 1 first** (fastest concrete result, exact same effort as adding a muscle), **layer
Tier 2 in in the same pass** if band force-extension calibration data is available (it usually is, from the
product spec sheet) since the API cost is identical. Defer Tier 3 until a concrete fidelity gap demands it.

---

## 4. Per-Joint Data / Anchor Plan

**Kinematics source (all joints, same pipeline):** IG corpus clip (or LabValidation clip) → MediaPipe/Pose2Sim
pose → OpenSim Scale+IK. Already-proven accuracy: ~4.7° lower-limb RMSE / 4.2–5.3° knee (540 comparisons,
matches the ~4.5° OpenCap anchor). Per-joint FORCE anchor differs sharply — table below is built directly on the
already-verified `MSK-JOINT-CONTACT-FORCE-ANCHORS` cell's sourced evidence
(`data/body_twin/agent_outputs/auto__a89669e49d76e0bfb.json`), not re-derived:

| Joint | In-vivo force anchor | Peak force (×BW) | Access | Confidence |
|---|---|---|---|---|
| **Knee** | Grand Challenge (D'Lima/Fregly, simtk.org/projects/kneeloads) **+** OrthoLoad/Kutzner 2010 — **doubly anchored, independent cohorts** | 200–346%BW (Kutzner stair-descend peak 346%; GC gait 180–300%, typ. 200–250%) | GC: fully open, no gate (already on disk, `knee_gc_measured/`, 47 trials). OrthoLoad: request-gated. | **Strongest** — matches flagship pattern exactly |
| **Hip** | OrthoLoad instrumented endoprostheses (Bergmann 2001) | 238–260%BW (walking/stairs) | Request-gated (1 free walking trial instantly; full multi-activity set = named-scientist request, non-commercial license) | Solid — single cohort, large/cited literature, tight range |
| **Shoulder** | OrthoLoad instrumented gleno-humeral (Westerhoff 2009/2011) | 105–283%BW (implant-type-dependent: anatomic vs reverse TSA) | Same OrthoLoad gate | Usable but noisier — pin the implant type per cell |
| **Spine** | OrthoLoad VBR (anterior, Rohlmann 2014) + internal fixator (posterior) | ~39%BW mean peak walking (VBR); absolute 1650N heavy-lift peak (**%BW unreported by source — any xBW figure here is an unverified estimate**) | Same OrthoLoad gate | **Weakest of the 4 anchored joints** — measures implant/segment load, not a ball-socket contact force; different physical quantity than hip/knee/shoulder |
| **Ankle (talocrural+subtalar)** | **NONE PUBLIC** — verified, nothing found in OrthoLoad or elsewhere | n/a | n/a | **No in-vivo anchor exists.** Fall back to a population/epidemiology anchor only (`MSK-ANKLE-TALOCRURAL`'s van Rijn/Valderrabano cohorts — re-sprain rate / OA incidence, not a force number), or an indirect estimate (foot-ground GRF + inverse dynamics moment, not a true joint contact force) |
| **Foot (mtp)** | **NONE PUBLIC**, same search | n/a | n/a | Same honest gap as ankle — no instrumented-implant program exists for the foot |
| **Elbow** | **NONE PUBLIC** — only a 2025 bench/cadaveric prototype instrumented humeral component, never implanted in a patient, no telemetry release | n/a | n/a | No in-vivo anchor; OpenBiomechanics pitching `poi_metrics.csv` has a **computed** `elbow_varus_moment` (large-N real pitchers, but derived via their own pipeline, not instrumented-implant truth) — usable as a *different, weaker class* of anchor, keep the distinction explicit |

**OrthoLoad access, precisely** (verified from the registry + the agent-output JSON, not assumed): free instant
download exists for exactly **1 hip + 1 knee walking trial** (small `.akf` files, email-gated form, no login).
The full multi-patient/multi-activity set for all 5 categories is **individual named-scientist request**, license
= non-commercial/scientific use only, mandatory citation of founder+institution+year+filename+URL. This is a
**human-owned-identity action** (matches the existing pattern for CAMS-Knee/SimTK) — the coordinator should file
it, not a tool.

**Elastic-band posterior-chain loading — no dedicated dataset exists (verified, §1.6).** Closest available proxy
on disk: `weightlift_grf/` (487 MB: `hp_obp.csv` hip-thrust+deadlift GRF/kinematics, `femoral_anchor/subj04.mat`,
`pone0251418/`) — ground-reaction-force during posterior-chain-loading barbell exercises, **not band-specific**
but the right muscle group/loading direction; usable to validate the PathSpring/Ligament force law's *plausibility*
(does simulated posterior-chain moment fall in the same regime as barbell-hip-thrust moment?) while a genuine
band video/force dataset is acquired separately (see §6).

---

## 5. Pipeline + First 3 Build Steps

**Full chain (already proven end-to-end on walking, per bt_memory):** corpus clip (ig_downloads or LabValidation)
→ MediaPipe/Pose2Sim pose (cad-to-simulation-I `.venv-humancap` CPU, verified live) → triangulation/filtering/marker
augmentation (Pose2Sim) → OpenSim Scale+IK → `.mot` → coordinate-bridge into the muscle-driven model (mapping
already proven for the JAM "DM" model, `video_bridge/`) → Inverse Dynamics / COMAK (needs synchronized GRF to
converge — verified hard requirement, IPOPT diverges without it) → validate the resulting joint contact/moment
against the §4 per-joint anchor, restricted to the force-plate-verified time window (the flagship's own hard-won
QC rule — comparing outside that window inflates error against an unmodellable stride).

**First 3 concrete steps for a coder to execute next:**

1. **Re-provision a clean OpenSim Python environment on ext4** (NOT on either NTFS drive — verified unsafe, §1.3):
   ```
   uv venv --python 3.12 source_repository/.venv-opensim
   PYTHONNOUSERSITE=1 uv pip install --python source_repository/.venv-opensim/bin/python opensim numpy
   PYTHONNOUSERSITE=1 source_repository/.venv-opensim/bin/python -c "import opensim; print(opensim.GetVersion()); print(hasattr(opensim,'Blankevoort1991Ligament'), hasattr(opensim,'PathSpring'))"
   ```
   (opensim wheel already cached, §1.1 — should take seconds, not the original ~30 min). This single venv gives
   both the model-selection pick (§2) and the elastic-material classes (§3, both core, no plugin) in one place.
   Separately, re-provision Pose2Sim the same way (`uv venv --python 3.10 ... && uv pip install pose2sim`) —
   also ext4, replacing the broken `video2kin/venv_pose2sim`.

2. **Attach a band `PathSpring` to the existing LaiArnold model and smoke-test in isolation**, before touching
   any real video: load `LabValidation_withVideos/subject2/OpenSimData/Mocap/Model/LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim`
   in the new venv, add one `PathSpring` from a new fixed ground body (`band_anchor_floor`) to an existing
   foot-segment path point, set `resting_length`/`stiffness` from a real product spec sheet, run `opensim-cmd`
   Forward or a short `Manager` integration through a prescribed ankle-dorsiflexion sweep, and confirm (machine
   check, not eyeballing a plot): tension == 0 while path length < resting_length, tension rises linearly (or per
   the Blankevoort curve) once stretched, and the reaction propagates up the kinetic chain in the ID output
   (nonzero hip/knee moment change vs the no-band control). This is the cheap decisive test that the elastic
   tier is wired correctly BEFORE spending effort on real video.

3. **Re-run the proven video→kinematics step on one LabValidation walking trial in the new environment** to
   confirm the rebuilt toolchain reproduces the already-measured ~4.7° RMSE (a regression check against a known
   number, not a fresh unverified claim) — then re-run the existing `video_bridge/comak_settings_video.xml` COMAK
   job against the new venv's opensim/Pose2Sim install to confirm the "video→knee-force 8.78%" result still
   reproduces post-rebuild. Only after that regression passes, extend the same pipeline to an ig_downloads
   plyometrics/squat clip (quality≥4, per the handoff's own priority) for the first athlete-corpus (non-lab)
   run, and to a resistance-band clip once one is acquired (§6).

---

## 6. Honest Gaps

1. **No resistance-band video/force dataset exists anywhere searched** (ig_downloads 893 clips: 0 band hits;
   `athlete_video_corpus`: 0 band-relevant; no dataset directory on the 8838D60F38D5FBDE drive is band-specific).
   The elastic-material engineering (§3) can be built and unit-tested today; validating it against a REAL athlete
   loading a band needs either a targeted new video capture/scrape or an operator-sourced clip — this is a genuine
   acquisition gap, not solved by this plan.
2. **Ankle/foot/elbow have NO in-vivo instrumented-implant force anchor anywhere in the literature** (verified via
   the existing `MSK-JOINT-CONTACT-FORCE-ANCHORS` cell's sourced search, cross-checked, not re-litigated). Any
   ankle/foot/elbow cert can only reach a population-epidemiology anchor (re-sprain/OA incidence) or a
   computed-not-measured proxy (OpenBiomechanics `elbow_varus_moment`), never a true measured contact-force
   number like knee/hip/shoulder. State this explicitly in any closing cell for those joints — do not imply
   parity with the knee flagship's evidence class.
3. **The already-built video2kin Pose2Sim+opensim pipeline is currently non-functional in place** (NTFS3
   symlink-resolution failure, §1.3B) — this plan diagnoses the exact cause and fix but does NOT itself rebuild
   the environment (left for step 1 in §5, a ~minutes-not-hours job given cached wheels).
4. **"Detailed foot beyond subtalar+mtp" (true midtarsal joint, multi-axis subtalar, per-toe) — no specific model
   verified this session.** WebSearch quota was exhausted (shared session budget, 2000/2000 used) before a live
   verification could be completed; I deliberately did not name a specific SimTK project/DOI from unverified
   recall (the same failure class the project's own memory flags — recalled identifiers are frequently wrong).
   If the operator wants beyond-standard foot fidelity, this needs a fresh live search, not this document.
5. **`ClutchedPathSpring` (control-dependent, needs an explicit "clutch" signal) was inspected but NOT
   recommended** — passive resistance bands have no clutch/engagement mechanism, so plain `PathSpring` or
   `Blankevoort1991Ligament` is the correct fit; `ClutchedPathSpring` is left noted only in case a future use case
   (e.g. a band with a quick-release) needs it.
6. **FEBio (Tier 3 elastic option) is verified absent from this machine** — not a blocker for Tiers 1/2, but if
   ever pursued it is a fresh toolchain acquisition, unestimated here.
7. **Spine and shoulder anchors are individually-gated (OrthoLoad) and NOT yet requested** — this plan identifies
   the exact access path and license terms but the request itself is a human-identity action (matches the
   project's own standing pattern for CAMS-Knee/SimTK), not something this plan executes.
8. **I did not independently re-verify that the stock (non-JAM-fork) `opensim` 4.6 pip wheel actually exports
   `Blankevoort1991Ligament` and `PathSpring` in Python** — I verified their C++ source location (core
   `OpenSim/Simulation/Model/` and `OpenSim/Actuators/`, not JAM-plugin-only) and separately verified their
   presence in the JAM fork's compiled tree; step 1 of §5 includes the exact one-line machine check
   (`hasattr(opensim, 'Blankevoort1991Ligament')`) to close this gap before relying on it.
9. **Model wrist articulation gap**: both candidate models' wrist is a `WeldJoint` (LaiArnold) or absent entirely
   from the joint list (Rajagopal) — 0 DOF. The operator's "ALL joints" ask is not fully met for the wrist by
   either on-disk model; flagged here rather than silently assumed solved (a literature/`bt_memory` note claims a
   "carpal/wrist" literature cell was built 2026-07-21 — that is a population/epidemiology graph node, NOT a
   wrist DOF in the simulateable model; do not conflate the two).
