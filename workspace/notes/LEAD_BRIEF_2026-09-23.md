# Management brief BodyTwin — mode 2026-09-23 ~23:30

Basis for guiding the work further. Written after Anton pointed out that the agents of the day were not based on the graph and the project's own code. Every claim has a source; "not checked" is written where applicable.

## 1. Where BodyTwin is located

| Plats | Vad | Regel |
|---|---|---|
| `~/projects/bodytwin` | BodyTwin's own repo (formerly "mechanism", fork of cad-to-simulation-I). Code, graph, ledger, project memory `bt_memory/` (2 859 files). Controlled by `COORDINATOR.md` and `docs/MECHANISM_HARDENED_CONVENTIONS.md` | Read from the workspace; will not change from here. Own isolation: memory in repot, neutral dialect in shared runtime state |
| `data/MECHANISM_ANCHOR_GRAPH.json` in the repo | Source graph: 3 950 nodes (3 873 OPEN, 66 ASSUMED, 7 REFUTED, 3 DEFERRED, 1 PROVEN). References `scripts/msk` 3 895 times | Starting point for each task |
| `~/projects/3fold-workspaces/bodytwin` | Work surface (adapter). Loading same graph (`./graph status/show/rank`), writing new work in `notes/`, `tasks/`, `results/` | Generated graphs (`GRAPH.json`, `MERGED_GRAPH.json`, `current/`, `generations/`) never edit |
| `/mnt/games-240/bodytwin_publication_2026_09_17/public` | Publicerad BodyTwin, commit b95a8dc | Patchar is in an isolated clone `/mnt/shared_data/bodytwin_work/R1/bodytwin_fix` (grenar fix/day1-findings, fix/gates-exit, fix/gates-exit-2) |

Code size in the repo (`.py`, machine inventory `results/MAP/private_inventory.tsv`, 15 292 scripts with docstrings): `scripts/physics_exp` 5 512, `scripts/tissuetwin` 3 617, `scripts/msk` 2 151, `scripts/cad` 1 541, `scripts/kernel` 546, `scripts/local` 534. Graph domains: `results/MAP/graph_domain_nodes.tsv` (796 nodes MSK/KNEE/PROSTH/TIS/VIDEO/DERM/BUILD/MODEL/SOLVE/FIX).

## 2. Existing work that the day should have built on

| Area | Graph nodes (status) | Code in repot |
|---|---|---|
| Knee, contact force against implant | `KNEE-CELL` (ASSUMED, n_eff = 1 against Grand Challenge), `KNEE-OG3-FMAX-CORRECTION-RESOLVE`, `KNEE-IDENTITY-BUCKET-DECOMPOSE-DM-JW-SC`, `MSK-KNEE-6DOF-JAM-ACL` (OPEN, proposes Lenhart 2015 6-DOF knee), `MSK-JOINT-CONTACT-FORCE-ANCHORS` | `scripts/msk/knee_jw_transplant/`, `knee_d1_*`, `comak_*`, `contact_waveform_*_cross_subject`, `jam_contact*` |
| Video → kinematics → force | `VIDEO-LINK1-KINEMATICS`, `VIDEO-TO-KNEEFORCE-ENDTOEND` (ASSUMED, with noted contradiction 2026-07-29), `VIDEO-TO-KNEEMOMENT-DROPJUMP-ACL`, `MSK-COM-GRF-VIDEO-ONLY` (8–17 % RMSE according to node), `MSK-KNOWN-LOAD-LIFTING`, `AUTO-MARKERLESS-VIDEO-KNEE-HIP-LUMBAR-CONTACT`; 53 nodes mention OpenCap/LabValidation | `scripts/video_index`, `scripts/mechanism`, `scripts/msk/corpus_ik*`, `body_composition_from_video` |
| Muscles, ligaments, whole body | `MSK-MUSCLE-ARCHITECTURE`, `MSK-DOF-ROM` (osim ranges are clamps, not measurements), `MSK-NEUROMUSCULAR`, `MSK-MUSCLE-SYNERGIES` | `add_arm_muscles`, `add_erector_spinae`, `add_hip_ankle_ligaments`, `add_knee_ligaments`, `add_scapula_clavicle`, `add_spine_ligaments`, `add_trunk_flexors`, `audit_muscles`, `fmax_pcsa*`, `contraction_dynamics*` |
| Cartilage, tendon | `MSK-CONNECTIVE`, `MSK-JOINT-MECHANOSTAT`, `MSK-TENDON-ADAPTATION` | `cartilage_contact`, `cartilage_poroelastic_relaxation`, `cartilage_thickness_material`, `achilles_*` |
| Bone, bone remodeling | `MSK-BONE-ADAPTATION`, `MSK-BONE-REMODELING-MECHANOSTAT`, `MSK-BONE-ENDOCRINE-HOMEOSTASIS`, `MSK-TRABECULAR-ARCH-FACTORIZATION` | `bone_remodeling`, `bone_rankl_opg_lemaire2004`, `bone_wolff_law_mechanostat`, `bmu_*`, `bone_stress*`, `bmd_femur_partial_recovery` |
| Soft tissue, body composition | `MSK-SURFACE-BODYFAT-DEGENBREAK`, `MSK-FAT-FUSION-IDENTIFIABILITY`, `MSK-VISCERAL-FAT-IDENTIFIABILITY`, `MSK-BIOIMPEDANCE-BODYWATER`, `MSK-FASCIA-NETWORK`; TIS domain 255 nodes | `scripts/tissuetwin` (layered shell, envelope, fibers, perfusion — built for another private domain), `calibrate_bodyfat`, `bodyfat_whtr` |
| Hand | `MSK-HAND-DEXTERITY`, PROSTH domain (97 nodes, grips, force-closure, tenodes) | `anatomical_hand`, `scripts/msk/inherited_hand/`, `scripts/hand_leg` |
| Physiology in general | ORG, NEU, ENDO, RENAL, CVD, IMM et al. (hundreds of nodes) | `scripts/msk/*` (cardiac, baroreflex, glucose, calcium_pth …) |

## 3. Dygnets resultat mot grafen

Source graph searched 23:15 for keywords (0 = no nodes).

**Nytt (ingen motsvarighet i grafen):**
- Shape models and morphs: VSD (0), Keast (0), femur-SSM (0), TPS morph (0) → P1, P2, P2b, P3, H7, RM1, J2, X1/X1b geometry core, D1.
- Mounts and Atlas Errors: LHDL (0), Pellikaan (0) → H2, H2b, AT1, JS1, IM3 (star movement against CT + prosthesis).
- Certified calculations: Lean/Range (0) → LV1, LV2 (+ BT-LV1-SUM from another session).
- Buggar i publicerad kod → R1, BT-R2 (patchar i klon).
- Oberoende granskningar AU1–AU8.

**Largely redone (Graph had it):**
- Video/GRF/OpenCap (53 nodes) → E1, E2, V0, V0b, V2–V6, G1, G1b, G2, FV1, BT-V7, V4. The numbers of the day can be booked as independent samples against existing nodes (eg V2 time-GRF 8,3 %BW against the node's 8–17 %).
- Knee/Grand Challenge (98 nodes) → X1b Test 4, IM3 used the same data.
- Bone remodeling → X3, X3b, X3c. Hand → HD1. Body shape/fat → V6.
- Foot contact → M2 (new adapter), M3, M4, FV2, FV3.

**Unchecked:** The C1/C1b contract against the 428 nodes that mention provenance/result_envelope; K1/K1b/F1 (hip reaction, recruitment, exhaustion) against the MSK nodes; K2/K3/H3/K5 (material) against the many density nodes — keyword match was too broad (`results/MAP/overlap_today_vs_graph.tsv`).

## 4. Verified and recorded

`notes/RESULTS_INDEX.md` (A1–A115, B1–B98; status A/R/G), `notes/DAY1_REVIEW.md` (final version + AU6–AU8 addendum), `results/D6/DEMO_CLAIMS_v3.md` (background for the collaborator 25/9), `tasks/graph_proposals/` (65 folds, 32 gaps, validated). Core findings for the collaborator: the gain lies in task-oriented metrics (T5) and averaged landmarks, not in the master model against the landmark-morph baseline (AU2, AU3, RM1); mount dominates the error and must be imaged (H2, H2b, AT1, AU5, JS1, IM3); regular motion carries almost no geometry information (IM1, IM2, IM3).

## 5. Rules learned during the day (applies to all new missions)

1. Each brief starts with searching the source graph (+ detail layer) and `~/projects/bodytwin/scripts/{msk,tissuetwin,physics_exp,cad,kernel}` as well as `bt_memory/`; the brief lists nodes and scripts to build on. Never "missing" without that search.
2. PREREG is hashed before first run; changes are logged last; review (AU) before "holds".
3. Load: `tasks/LOAD_RULE.md` (nice 19, 2 threads, ≤ 2 processes, wait on load/20 > 1,0; ≥ 4 GB or GPU via `tasks/heavy_run.sh`).
4. Disk: / 1,2 GB, shared_data 8,9 GB, sdc1 17 GB (BodyTwin roof 12 GB there). No new downloads without space.
5. No `import *` from other agents' code; own OUTDIR (incident H2b overwrote H2).
6. lane_runner is not used tonight (Anton 23:00). swarm_worker lanes: `--agent build` + mandate required (otherwise plan mode). coordinator weekly quota 5 % (reset 22:00).
7. Inga mail, ingen publicering, ingen trial/Wine; restricted model data and the collaborative's data interna; LHDL icke-kommersiell.

## 6. Ongoing right now

coordinator: AU6 (reviews IM1/IM2/FV3), LV3 (error limit controls resolution), K5 (density–modulus relation). Cron 03314fad monitors (:23/:53), nothing starts.

## 7. Proposal for first step tonight (for Anton's decision)

1. **Overlap map by area** (knee, video/GRF, muscle/ligament, cartilage/tendon, bone, soft tissue, hand, material): read the relevant nodes and scripts, and for each daily result, enter the node/script, matches/contradicts/extends. Result: a table in `notes/`, and graph proposal (`tasks/graph_proposals/`) that connects the folds of the day to the correct existing nodes instead of new ones.
2. **Then build on existing cells**, not in parallel: e.g. the knee cell (n_eff = 1) with the geometry of the day (shape model, mounts, corrected HJC) as a new, independent leg; the video nodes with V3/V5/G1 as samples; the tissue twin shells as soft tissue layers over segments.
3. **the collaborator's track** (new and reviewed) continues from D6/AB1b; external solver runtime is waiting for license/Wine.

## 8. Addition 00:45 — correct input and MAP2 output

**Correct input (BodyTwin's own tools, in `~/projects/bodytwin`; NOT used during the day):**
`python3 scripts/mechanism_preflight.py` (PASS 3 950 nodes; authority `docs/MECHANISM_HARDENED_CONVENTIONS.md`; memory `bt_memory/` via `context_loader.py --memdir bt_memory`) → `python3 scripts/tool_find.py <ord>` (8 429 script; "video" 1 671 hits, "pose" 4 797, "athlete" 14, "opencap" 0 — the word is missing but the work is under LabValidation) → `python3 scripts/capability_index.py --task "<ord>"` (status-stamped) → `anchor_graph_tools.py priority|next-actions`. Split fleet preflight shows 8 dead demons (grade_watcher et al) — fleet thing, not restarted.

**MAP2 (results/MAP2/A,B,C — overlap.tsv 25+46+45 rader, SUMMARY.md):**
- My "missing" was wrong for: form priors (MSKEL-GEOM-OVERDET), contract (result_envelope_v1 in published code), rainflow/Miner (`scripts/fatigue_life.py`), plastic return map, grip (same 45–134 N).
- Extends existing nodes: KNEE-CELL (X1b), MSK-COM-GRF-VIDEO-ONLY (V2/V3), VIDEO-METRIC-SCALE-REFERENCE (G1/G1b), HOLE-BONE-MODULUS-DENSITY-EXPONENT-GEOMETRY (K2/K3/K5), TIS-PERIODISATION-NEEDS-DAMAGE (X3c).
- Never run against each other: existing hip chain 386,77 %BW (1,41× OrthoLoad) against K1/H2/J1.
- **Unmapped:** ecosystem-wide athlete video (twin_capture_lane, sharp_football, pitch_awareness, cad-to-simulation lanes: forceplate_jump_cell, smu_grf_kinematics_xcheck_cell, radar_speed_validation_cell, biomechanics_cert_leg_v0 etc.) and if these pipelines work today. Machine file inventory in progress → `results/MAP/video/video_scripts_all_projects.tsv`.
