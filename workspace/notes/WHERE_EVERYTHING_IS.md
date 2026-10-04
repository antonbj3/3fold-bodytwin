# Where everything is — orientation map for BodyTwin

Written 2026-10-02 by the coordinator session, in response to Anton’s question: *"Don’t know where anything is, or whether I can trust whether you miss things or not"* and *"will it come along to the new repo?"*.

**Why it exists.** The path has gone cad-to-simulation → `~/projects/bodytwin` → `~/projects/3fold-workspaces/bodytwin`, and each move carried some of the work onward, not all. This file should come along at the next move. Everything here is **measured today** unless stated otherwise, and what I have not checked is marked UNCHECKED.

**Why it is needed beyond a good memory.** I drew two false absence conclusions on the same evening it was written: a search with `grep -rl <ord> *.py` did not recurse and missed 44 files in five subpackages, and a nervous system search used cortex/axon/synapse and missed motor unit/recruitment/EMG, which is where our nervous system work actually is. Both were caught, one by Anton. The map is cheaper than trusting search terms.

---

## 1. Three trees, not one

| where | what | measured |
|---|---|---|
`~/projects/3fold-workspaces/bodytwin/` | **the workspace.** Adapter reading source projects without changing them. `tasks/`, `results/`, `notes/` live here. | 43 cell families, ~15 800 result directories |
`~/projects/bodytwin/` | **the older main project.** The center of gravity is here. | 5 329 scripts in `scripts/physics_exp/`, 164 data directories, 1 309 completed work cells |
`~/projects/bodytwin/scripts/physics_exp/` | experiments whose **results are beside the script** as `*_evidence.json`, not in any `results/` directory | — |

Consequence: *"it does not exist"* must never be said after searching only one tree.

Public: `~/projects/3fold_staging/` has nine repos, including `3fold-bodytwin` with 757 files and **5 model files**, and `3fold-motion-engine` with 1 869 files and 199 model files. The publication check is pinned at commit `b95a8dc`.

## 2. The completed work cells — 1 309 of them

Indexed in **`notes/OLD_CELL_INDEX.json`** (search with `python3 tasks/build_night/index_old_cells.py --quantity <text>`).

- 1 620 cell records, **1 309 distinct cell IDs**, 1 566 with a numerically decisive number
- status in the source: **696 MEASURED, 431 REFUTED**, 27 WITHHELD, 9 VOID, 6 HYPOTHESIS-awaiting-QC, 431 without status
- domains: `bonemet_v1` 1 276, `claims` 308, `glioma_fisher_v1` 28, `hfpef_diastolic_v1` 5, `csf_davson_icp_flux` 1
- **the 431 refuted are as useful as the measured** — a preserved negative prevents a lane from rerunning a question that has already failed
- NOTE: searching for flow, tissue, perfusion, collagen, blood, vessel and pressure gives **zero** hits in the index. The cells belong to another research line (amyloid kinetics, AMYKIN series). Do not count on them for tissue and flow.

## 3. The workspace’s 43 cells

`tasks/free48/sources/` — 35 Q-numbered plus BIORESP, IMMUNITY, MITOSTRESS, SOLBENCH and four SURG_*. Each cell is a first-principles chain with declared parameters, units and dimensional checks (see `Q012_model.py`, which starts "first-principles model of diet x microbiome -> metabolite exposure").

**Numbering is sparse up to Q168 but there is no queue of unbuilt cells** — `notes/RESEARCH_ACTIONS_20260926.json` contains the same 35 Q numbers. Nobody has written down which cells the body *should* have. That is the open map gap.

**Corrected 2/10 23:05:** `SURG_COLLAGEN` and `SURG_HEMOSTASIS` did NOT LACK code — 951 and 1102 lines respectively existed, and my "no cell" was based on the filename pattern `<CELL>_model.py`. Both now have an entry point with the cell contract and run. Remaining without a model file: `BIORESP`, `IMMUNITY`, `SOLBENCH`.

## 4. The surgery substrate (which Anton wants DENSE)

| cell | lines | carries |
|---|---|---|
`SURG_INCISION` | 533 | incision, tool-specific toughness |
`SURG_HEALING` | 638 | healing, regional flow |
`SURG_COLLAGEN` | 951+271 | triple helix Tm, tendon hierarchy, myofascial transmission + the cutting chain molecule→tissue 
`SURG_HEMOSTASIS` | 1102+330 | thrombin generation, platelet plug + incision geometry→bleeding rate with identifiability certificate |
`Q017` | 829 | renal blood flow |
`Q140` | 458 | renal plasma flow |
`Q036` | 638 | normalized regional flow |
`Q107` | 869 | solid reference flow |
`Q021` | 550 | pulmonary blood flow |
`Q156` | 475 | directional pressure drop |
`Q090` / `Q127` | 1 025 / 617 | load-bearing tissue, fiber and satellite cell |

Undrawn edges existing as possibilities: `Q017`↔`Q140` are separated only by **hematocrit**; `Q036`↔`SURG_HEALING` share regional flow; **CORRECTED 2/10 23:30:** `Q156` does NOT couple pressure drop to the flow cells. Its `pressure_delta_mmhg` is an INPUT with default 0, and concerns Darcy flow in brain parenchyma at κ = 1,9e-15 m² — not pressure drop in a vessel segment. Wrong level and wrong quantity. My original row was thus wrong.

## 5. The nervous system — modeled as periphery, not center

Correction of my own incorrect claim the same evening. With the right vocabulary (motor unit, recruitment, EMG, efferent, proprioception, reflex):

**Q043 89 hits** (motor unit recruitment), Q017 38, Q022 25, Q080 24 (receptor cooperativity, desensitization), Q044 21, Q058 14 (coarse/fine muscle model switching), Q031 13, plus Q021, Q088, Q090, Q100.

Old project: `baroreflex`, `cough_reflex`, `emesis_reflex`, eight `msk_*` domains, `prosth_anthropometry_sensitivity`, `HAND_CALIBRATION_V0`, `HAND_MANIPULATION_DEMOS_V0`, plus `amygdala_fear_circuit`, `basal_ganglia_gating`, `blood_brain_barrier`, `dopamine_kinetics`, `serotonin_system`, `synaptic_vesicle_release`.

**The periphery is well covered. The central and cortical are thin.** For surgery, the periphery is the relevant part.

## 6. Mechanism

`~/projects/mechanism`: 7 690 py files, 5,4 GB. Transfer needs are in `notes/MECHANISM_TRANSFER_INVENTORY.md`.
UNCHECKED in detail. What I checked: searching for look/gaze/saccade/fixation gives 37 files, but they are in grasping and optics code, i.e. camera gaze and not eye. **No nervous system material there that the workspace lacks**, as far as I have seen.

## 7. Direction

`tasks/build_night/` — `COMMON.md` (read by each lane: classify the claim, resolution level per quantity, edge at the finest level, key per quantity, search the cell index before derivation), `slots.txt`, `max_lanes`, `steer/<LANE>.md`, `NIGHT_LOG.md` (each finding one row, historical rows are never edited).

Meters: `lane_runner_usage.sh` (window **and credits** — the window is no stop when credits exist), `scan_swarm.py`, `verdict_emission.py`, `index_old_cells.py`, `find_vacuous.py`.

## 8. What is actually broken, measured

- **0 of 202 robot ports** carry a computed quantity; 167 require new computation
- **27 proven signs have 0 consumers**; 25 of the anchors have no decisions at all
- **3 of 903 directory pairs are decidable** — and the metric "0 uncrossed pairs" measured *never mentioned together*, not comparability
- external models are executable but **never run**; the only actual connection was rejected at 0 % valid frames
- `negative_result` exists in 47 of 11 799 results (but 77 % in jobs after the requirement was introduced)

