# DELIVERY — BodyTwin-lanen (22–24 sep 2026)

Status: **ongoing** = running/not committed · **commit** = saved locally · **reviewed** = independent audit (G in `notes/RESULTS_INDEX.md`) · **candidate** = candidate for publication (none yet; requires Anton's decision). A commit is not scientific approval. Nothing is pushed.

## Where the work is

| Del | Plats | Versionshantering |
|---|---|---|
| The lane's code, PREREG, result receipts, notes, briefs | `~/projects/3fold-workspaces/bodytwin/{notes,tasks,results}` (shared workspace, not git) | local repo with separate git directory `~/research/bodytwin_lane_git`, branch `bodytwin-lane`; run `tasks/lanegit.sh <git-kommando>`. Only explicitly added files are tracked; other sessions' files (DS_*, BT-MAT-*, BT_*_2026*, PROOF_LANE_*, GRAPH_*, BT-SUM-GRID1, BT-CERT-AUDIT, BT-LV1-SUM, SOL6_*) are not included |
| Patches to published BodyTwin (R1, BT-R2, BT-R2b) | clone `external_mount` (branches fix/day1-findings, fix/gates-exit), temporary clone `external_media` (branch fix/gates-exit-2, commit cde20a2); format-patch in `results/BT-R2b/patches/` | own git commits in the clones; published b95a8dc untouched |
| Private BodyTwin repo `~/projects/bodytwin` | only read | nothing written from here |
| Cloud infrastructure (lease extension) | `~/research/sol6_recovery_20260923/CLOUD_HUNT_20260923/infra/extend_lease_20260924.sh`, `*_worker.json` (lease_history) | coordinator's directory, not in the lane's repo |
| Large data and run outputs | see `tasks/index/DATASETS.json` (source, licence, path, SHA256) and `tasks/index/index.json`; `external_mount`, `external_media` | not in git (npz/pkl/h5/stl, files > 2 MB, `inputs/` copies) |

## Commits

| Commit | Contents |
|---|---|
| 98ffa97 | anteckningar, planer, briefs, lane-drivare, `cloud_run.sh`, `lanegit.sh` |
| 7973152 | avslutade experiment: kod, PREREG, resultat-JSON/MD, granskningar (2 118 filer totalt i grenen) |

## Leveranslista

Numbers and conditions: `notes/RESULTS_INDEX.md` (row ID in parentheses). Review = AU row.

| Lane/experiment | Commit | Test / outcome | Review | Remaining gaps |
|---|---|---|---|---|
| D1, X1, X1b geometry core (rotation-invariant folding, geomgr) | 7973152 | `cd results/D1 && pytest -q` 11 passed; X1 16 passed; X1b 41 passed (24/9 00:4x) | AU1/AU3 | gonial +60 passes published certificate (MAP2-B) → BT-N5 ongoing |
| D5/D6 demo for the collaborator | 7973152 | `python3 results/D5/demo/run_demo.py --subject z001` 8/8 (22/9); rerun 24/9 waited for load gate, not completed | AU2 | DEMO_CLAIMS_v3 needs updating with the night's results |
| P1, P2, BT-P2b, P3, H7, J2, RM1 (shape model/morph) | 7973152 | script + JSON per folder | AU2, AU3, AU8 | no significant gain against the landmark-RBF baseline → N1 ongoing |
| H2, H2b, AT1, JS1, AU5 (attachment atlas) | 7973152 | script + JSON | AU4, AU5 | attachment regions → BT-N15 ongoing |
| IM1/BT-IM1, IM2/BT-IM2, IM3 | 7973152 | JSON | AU6 (A120): numbers reproduced, conclusions reformulated | IM2 = untested; IM1 percentage driven by injected error |
| K1, K1b, J1, F1, H5, H5b, D3, AB1, AB1b | 7973152 | K1 tests JSON (0,012 %) | AU7/AU8 partial | external solver runtime not run (licence/Wine) |
| C1, BT-C1b, H6, C2, C3, C4 | 7973152 | BT-C1b 9/9 pytest in BT-R2b (A117) | AU8 | graph feedback not registered |
| BT-LV1, BT-LV2 (Lean/interval) | 7973152 | `results/BT-LV1/lean.log` compiles | AU7 (A113) | "error bound governs work" → LV3 ongoing |
| R1, BT-R2, BT-R2b (patches to published code) | clones + 7973152 (receipts) | 16 regression errors before / 16 passed after (A117) | AU8 (BT-R2) | BT-R2b unreviewed; not placed in the shared clone |
| K2, K3, K5, H3, H3b, X3, X3b, X3c (materials/bone) | 7973152 | JSON | AU8 (H3b, X3c) | K5 unreviewed (A116); X3 parallel Pivonka (MAP2-C) |
| Video/GRF: E1, E2, V0–V6, G1, G1b, G2, FV1, FV2, BT-FV3, BT-V7, M1–M4, HD1 | 7973152 | JSON; M2 pytest requires the pinocchio environment (M1's venv) — 1 collection error in default python | AU6 (FV3), AU8 | video parked (Anton 24/9); MAP3 map |
| MAP, MAP2, MAP3 (mapping against graph/code) | 7973152 | TSV | — | not graph-bound |
| BT-B1, BT-B4 (Space The swarm) | 7973152 | BT-B1 `recalculate.py` 2,94 s | not reviewed | — |
| **Ongoing (not committed):** N1, N2, N3, N6, N7a, N7b, N12, N14, LV3 (coordinator); BT-N4, N5, N8, N9, N10, N11, N13, N15 (swarm_worker); BT-B5–B10 (Space The swarm) | — | — | review planned per result | committed when complete |

## Graf/journal

Proposal for graph connections (folds, artifacts, links, 32 gaps) in `tasks/graph_proposals/` (validated, not loaded). Registration through the workspace's `./graph working dispatch/feedback` (`notes/GRAPH_WORKFLOW.md`) remains — the results must go in with status PENDING_INDEPENDENT_REVIEW. Negative results and validity limits are in RESULTS_INDEX and each RESULTS.md.

## Before publication (waiting for Anton)

