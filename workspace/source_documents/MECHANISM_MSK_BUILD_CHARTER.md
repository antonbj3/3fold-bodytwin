# MECHANISM MSK BUILD CHARTER — the mechanical-twin phase (2026-07-21)

## The pivot (operator, genuine, 2026-07-21)
Shift from the literature cert-graph (**Phase 1, DONE**: 994 nodes, 55% inter-cell connectivity, guard green)
to **BUILDING a loadable biomechanical musculoskeletal twin**. Operator wants:
- **Correctly-modeled muscles + ALL joints INCLUDING feet** — detailed multi-segment foot (talocrural + subtalar
  + MTP), not a lumped 1-segment foot.
- Ability to **handle ELASTIC materials** — e.g. resistance/rubber bands used to **load the posterior chain from
  the ground up**, force propagating through the kinetic chain (gymnasts loading against bands).
- Grounded in a real athletic-movement corpus (gymnasts, foot/ankle drills, plyometrics, weighted lifts) = "a goldmine".

## Mandate (explicit)
"fix it on your own... I can be helpful if it gets difficult." → **AUTONOMOUS build.** Acquire whatever data is
needed (OrthoLoad / in-vivo joint contact forces, OpenSim models, etc.) on my own; figure out the optimal per-joint
approach; build it. (`[[feedback_cad-to-simulation_build_every]]`: build all autonomously when the mandate is clear.)

## FIRST-STEP DOCTRINE (operator, 2026-07-21) — datasets/models are bootstrap, NOT the ceiling
Every dataset/model in this charter (LaiArnold2017's **3-joint** foot, OrthoLoad's **published summary** numbers, the
IG **auto-labeled** corpus) is a **FIRST STEP**, explicitly not the target. The trajectory is progressively higher
fidelity that **we generate**: a per-bone (~26-segment) foot not a 3-joint one; un-welded wrist/hand; subject-specific
geometry from imaging + the filmed corpus; richer elastic-material models; and ultimately the twin generating better
training data than it consumed (the data-generation-for-data-scarcity thesis). **BUILD HOTSWAPPABLE:** every pipeline
stage — skeletal model, per-joint anchor dataset, kinematics source, elastic-material model — must be a swappable
module so a fidelity upgrade is a *dial*, never a rebuild. Never treat a current dataset/model as the answer.
Operator: "you solve it" — solve it, keep cranking the fidelity.

## Approach: BOTTOM-UP
feet/ground-contact → ankle (talocrural + subtalar) → posterior chain (gastroc–soleus–Achilles, hamstrings, glutes,
erectors) → **knee (flagship, already anchored)** → hip → full skeleton + muscle-tendon actuators + elastic-band
actuators. Rationale: matches the operator's force path (band loads the posterior chain from the floor and propagates
up), and the KNEE-CELL flagship already anchors the middle of the chain.

## Data
- **Movement corpus:** `~/ig_downloads/` (893 clips, ~2 GB). Handoff: `~/ig_downloads/BODYTWIN_HANDOFF.md`.
  Semantic rows: `ig_athletics__semantic.jsonl` `{video,category,what,group,gt_channels,quality,...}`. Optics DB
  `ig_capture.db`. Groups: plyometrics_jumps 365+97+26, weighted_lifts/squat 56, None 313 (VLM-review).
  gt_channels: bodyweight-motion 893, slowmo-rate 94, known-load 71. **Labels are AUTO/unverified** → VLM label-review
  before finetune-trust. foot/ankle/mobility = `closeup-low-framemotion-likely-isolation-KEEP` (do NOT discard — the
  high-value foot-model data).
- **Pose pass:** run FROM `~/projects/cad-to-simulation-I` (`.venv-humancap` CPU / `.venv-newton` GPU; **bodytwin lacks
  these venvs**). MediaPipe model: `/media/anton/183E48713E484A48/cad-to-simulation_video/geoprobe_scratch/weights/pose_landmarker_full.task`.
  Priority: `quality>=4` + `plyometrics_jumps/*` + `weighted_lifts/*` + `slow_motion/*`.
- **Joint-force anchors (autonomous acquisition mandated):** OrthoLoad (Bergmann/Charité — in-vivo instrumented-implant
  hip/shoulder/knee/spine; **gated individual-request**, see the OrthoLoad entry in `data/MECHANISM_SOURCE_REGISTRY.json`
  + `[[anchor-graph-tools...]]`) + Grand-Challenge knee (`simtk-knee-grand-challenge` / CAMS-Knee). Decorrelated
  in-vivo force anchor per joint; honest where none exists (detailed foot, elbow — already flagged in the graph).

## Substrate
OpenSim (the KNEE flagship's tooling: opensim 4.6 installed, IK/ID; `simtk-knee-grand-challenge`, opensim-jam,
gait2392/Rajagopal, addbiomechanics). Full-body-muscle + detailed-foot model + elastic path-actuator choice → decided
by the verified capability inventory in **`docs/MECHANISM_MSK_BUILD_PLAN.md`** (written by the plan agent). For rich
elastic-material fidelity, evaluate OpenSim PathSpring/elastic-actuator vs a soft-body/FEBio path — plan decides.

## State / handoff (post-compaction: work Phase 2, not the lit-loop)
- **Phase 1 lit cert-graph: DONE** — 994 nodes, 55% connectivity, guard green. The steady-burn hypothesis-cell loop
  is **PAUSED**: on cron ticks, still `mechanism_autofold_pending.py` + guard + commit to keep the graph durable and
  fold any in-flight cells, but **do NOT fire new `watertight-researcher` hypothesis agents** — the direction shifted.
- **Phase 2 MSK build: ACTIVE.** Concrete next steps live in `docs/MECHANISM_MSK_BUILD_PLAN.md`.
- Couples to existing graph cells: `KNEE-CELL` (flagship), `MSK-JOINT-CONTACT-FORCE-ANCHORS`, `MSK-ANKLE-TALOCRURAL`,
  `MSK-TENDON-SPRING-SAFETY-FACTOR`, `MSK-HIP-FEMOROACETABULAR`, `MSK-SHOULDER-GLENOHUMERAL`, the OrthoLoad thread.
