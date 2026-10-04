## Current coordination after the coordinator limit — 2026-09-23

Anton's latest instruction replaces the model choices in the older queue below: **the build lane medium**, swarm_worker only when quota is available again, the proof lane xhigh only for scoped difficult problems. Global status: `external_research_path`. Previously completed jobs must not be restarted from this historical list.

- The form/function proposal further down has been processed in BT-IM1 and BT-IM2; see RESULTS_INDEX A106/A108. Its older formulation about unacknowledged reading describes the handover time.
- Certificate review and BT-LV1-SUM are delivered and recorded with limitations. No complete guarantee for the entire geometry chain is claimed.
- Next prepared innovation task: [BT-MAT-OBS0 — observable material measurements](SOL6_MATERIAL_OBSERVABILITY_20260923.md). **Started through swarm_worker V4.1 Go around16:11**, own `results/BT-MAT-OBS0/`; contract and nullspace first. Separate review required before approval.
- Ongoing resumptions in sibling tracks: dental GEOM_UNC and DESIGN_crown and the field's isolated C1 aggregation. Resource locks still apply.

## Night 23–24/9 step 1 (Anton 23:20 "Start with 1") — MAP2 overlap map, brief tasks/MAP2_BRIEF.md
- MAP2-A (coordinator-medium, started): knee/video/GRF/foot/hand/fat
- MAP2-B (coordinator-medium, started): geometry/attachments/muscles/contracts/certificates — tests the "new" claims
- MAP2-C (coordinator-medium, started when K5 finished): bone/materials (K2, K3, K5, H3, H3b, X3*, BT-MAT-OBS0/1, BT-MAT-SUM-AUDIT1, BT_PLASTIC_HISTORY_ORDER), soft tissue/tissuetwin, the DS_* results from the sibling session, BT-SUM-GRID1, BT_HISTORY_MINIMAL_STATE_AUDIT, GRAPH_*

# Agent queue — BodyTwin (Anton 22-09: ~5 coordinator → 21:40 "step on the gas" → 7; break down into the smallest components, Wine = addition) — 10:58 Anton: weekly quota 88–89 % → NO new coordinator agents until 22:00; new work as swarm_worker lanes (tasks/lanes/run_lane.sh), coordinator review after 22 — 23:00 Anton: NO lane_runner tonight; no new agents; wait for instructions

## To the coordinator coordinator — handover at Anton's request 2026-09-23

Read [the innovation proposal after RM1: form and function](../notes/INNOVATION_FORM_FUNCTION_2026-09-23.md). Core: models that fit the same landmarks but give different moment arms → choose distinguishing movements/measurements → test joint estimation of joint centres, muscle attachments and muscle paths → validate against the same information and independent ground truth. First compare with already completed H2/H2b and reviews. Proposal left in the queue; no new run started and reading not yet acknowledged.

Every finished agent is replaced immediately. Root reviews output before anything is recorded as "holding".
Shared rules in every brief: PREREG must be hashed (PREREG.sha256) BEFORE the first run — AU1 found that several PREREG lack a pre-run hash; INNOVATE EVERYWHERE (Anton 00:40) — every brief has an innovation part: decomposition → new combination or cross-domain transfer with explicit mathematical mapping → falsifiable test; RESOURCE RULE (proportionate, Anton 22:25: balance safety against performance): only jobs with ≥ 8 GB RAM or GPU go through tasks/heavy_run.sh (bigmem.lock/gpu.lock/resource_gate/MemoryMax); others run directly with thread cap 4 and estimated memory; pilot cases only for unknown large memory; large intermediate results on external_mount<id>/; EVERY number in the final answer must exist in results/<id>/*.json with a source (Anton: lose no data points; root records them in notes/RESULTS_INDEX.md); read-only in source projects; new files only under its own output file / results/<id>/; data to external_mount (/ has ~5 GB); CPU thread cap 4; no `nvidia-smi -q`; `ps -eo pid,args`, not `pgrep -f`; no email/publication/trial/Wine; source vs derivation vs hypothesis separate; UNKNOWN explicit; no own subagents.

## Wave 0 — context (started ~21:15)
- O2 BodyTwin geometry: actual chain published/staging/private, parameterisation, landmarks/attachments, local data → notes/GEOMETRY_CAPABILITY_AUDIT.md (+ reads dental I1 when complete)
- O3 Field/engine results relevant to BodyTwin (A183–A208, U285–U288, SHARED_GEOMETRY_MATH) → notes/FIELD_RESULTS_FOR_BODYTWIN.md with hypothesis candidates
- O4 Researcher scenarios + modalities → notes/RESEARCHER_SCENARIOS.md
- O5 Mechanism transfer → notes/MECHANISM_TRANSFER_INVENTORY.md

## NEXT FREE SLOT
- E1 (interrupted at the crash before output): restart with the same brief when a slot exists
- P2 Keast tibia SSM (keast2023_tibia_ssm) as an independent control of the P1 method

## Wave 1 — candidates (filled from wave 0)
- D1 Demo chain: real public mesh → shape parameter → landmarks/attachments/regions tracked → measurable consequence in a calculation (moment arm / mass/inertia / contact surface) — PREREG + counter-test
- D2 Dataset: public shape model (SSM) / segmented bones with a licence — acquire selectively
- B1 Geometric uncertainty → computational result (sensitivity per shape parameter; what resolution is required)
- B2 Thin regions/material partition through conversion (A194 lesson) with anatomical recipient
- E1 (from O5) Camera uncertainty → segment length → scaled inertia/moment arm; counter-test = coverage (corner view); baseline fixed model + size prior
- B3 Combined interventions — contract + smallest executable example

## From O3 (notes/FIELD_RESULTS_FOR_BODYTWIN.md §4)
- H1 → B1 (started ~21:35)
- H2 identifiability/OED against a biomechanical quantity (baseline: simple anchor at the attachment, A206; global multistart)
- H3 thin layers: volume+first moment per material + β1, laminate stiffness (isotropic-mean baseline; tunnel test)
- H4 moment-arm sensitivity through wrapping events (shape parameter → moment arm) — requires D2 data
- H5 exact multi-material inertia through second moments + affine deformation — the J track, cheap
- H6 fusion refusal with undeclared variants; H7 task metric vs variance ordering in shape space; H8 Q2 trace patch (engine)

## From O4 (notes/RESEARCHER_SCENARIOS.md)
- X1 Attachments+joint centre with stable IDs on bone mesh → moment arm → knee contact force vs Grand Challenge implant force; counter-test: geometry-specific moment arm must beat the difference between two generic scalings (= demo core D1, requires D2)
- X2 AddBiomechanics: contact time+leg stiffness from RAW marker trajectory vs force plate (not dynamically adjusted trajectories — shared source)
- X3 Mechanical load → RANKL/OPG + bone mass (Lemaire); counter-test: interaction ≥10× baseline drift, decoupled arm = exactly 0
- V1 Verify/fix synthetic-flag propagation in bone_remodeling/bmu_turnover_kinetics (+ gates False but exit 0) — in own copy, report for publication decision
- V2 Doc fix: docs/RUNNING 'pytest -q tests' requires PYTHONPATH=src (5 collection errors) — report for publication decision

## From O1/D2 (~22:20)
- P1 started (demo 1: population → individual from a few measures)
- D3' the collaborator's attachment points (box-lifting output, 136 thigh attachments, INTERNAL) → BodyTwin frames → against a real femur (demo 3); compare with TLEM attachments
- D4 Demo 2 "safest": real femur through the BodyTwin chain with surface distance (Hausdorff/RMS), not just volume
- PUBLICATION GATE (does not block internal work, Anton 22:30): TLEMsafe agreement, the collaborator's box-lifting attachments, OpenMandible licence — decided only at sharing/publication. Keast 2023 tibia SSM (CC BY 4.0) may be acquired when needed without asking. the reference model trial/Wine tomorrow.
- M2 Motion-engine adapter in results/ (not the engine's checkout): osim/the reference model → URDF with coupled joints, floating base, external forces; foot-contact model for GRF/events (K1 leaf)
- K1b K1 sensitivity extended to ±40 mm (nonlinear) with p=2, + hip-centre definition as uncertainty (27,8 mm between marker definitions, V0) — V0 used the p=3 file
- V candidates 4–8 from results/V0/candidates.json (segment length→mass/inertia, calibration→uncertainty, DXA prior, bar force, re-ID)
- external_media Volume disconnected after the crash — not inventoried
- C3 GRF and kinematics filtered with the same cutoff (Kristianslund) → peak-moment stability 4–8 Hz; + functional hip-centre calibration as the largest budget term (C2)
- M3 Foot contact new PREREG on another trial/person (OpenCap LabValidation has force plates): floor per sphere + toe point from M2 post hoc; HS delay 43–47 ms is open
- V0b candidate 5: assembly order from image (the NH35 pair), Kendall τ ≥ 0,3; hands/grip/force from video (hand biomechanics for BodyTwin)
- New Volume (EuRoC, 85 GB benchmarks) — inventory when connected again
- H3b curvature-adaptive refinement (radius < cell → refine) in the zone 0–10 mm below the joint surface; target stiffness+subchondral strain ±5 % at lower cost than uniform 0,5 mm
- FV2 (after the field session's audit A225): moving body in the field engine through U318 → possible occlusion/physical plausibility (foot–floor, bone–bone) in IK; FV1 showed that occlusion explains only 1–7 %, so priority on foot–floor (M2/M3/M4: feet 14–33 mm above the plate)
- J2: geometry manager with a population prior (P3) in the collaborator's box lifting instead of raw TPS (J1 failed due to noise amplification) — distinguishable from affine?
- R2b: reclassify tcell/mapk/hpa (and review the other 79 B) according to AU8; add C1b's 9 counter-tests to pytest
