# Handoff: field engine and innovation hunt — 2026-09-22, around 20:20 Stockholm

## Anton's latest instruction and responsibilities

Anton wants a coordinator session to take over this session's work on the field engine and ongoing research **before work on the graph merge continues**. No new research jobs or integrations have been started during the handoff. A subagent just assigned to write a fold reader was interrupted; no `fold_reader.py` or associated test files had been created at the check.

The main task is a broad innovation hunt: the field engine, physics/materials/sand, adaptive resolution and DOF, observation/OED and geometry; P1–P4 are frozen published seeds for expansion, also beyond the original domains. Work from your own measured results and errors, first principles, testable hypotheses and primary literature. A literature compilation alone or product cleanup does not meet the goal.

Another session owns the graph engine's research and kernel work. Do not change its checkout or research graph. Dental and BodyTwin must be prepared with full context; do not choose a specific medical application yourself based on Anton's needle or incision example. FPGA implementation is paused.

Anton says the coordinator quota was just reset and that he has about 27 hours to use it, with a new available model. He wants to use coordinator and let swarm_worker rest for a while. Model access/aliases are not verified by this handoff. coordinator Code has been updated and verified to **2.1.280**, executable via `local_config_path/bin/coordinator`. No new coordinator session has been started for him.

## Read first, without rereading all historical material

1. This file and `~/research/field_handover_20260922/LANE_INVENTORY.md`.
2. The research root below: `lanes/INNOVATION_OBJECTIVES_2026-09-22.json`, `lanes/LANE_PROTOCOL.md`, `lanes/SEED_DECOMPOSITION.md`.
3. `HUNT.md`, relevant recent rows before `## Closed`; then the reports for the branch you take forward.
4. `~/HANDOVER_romi_2026-09-21.md` is the older chronological log. Its latest addition was at 18:09 and describes the queue running then. It is **outdated on current process state**; use this new handoff and a live process check.
5. For the graph coupling: `~/projects/graph_workspace/GRAPH_MAP.md` and `lanes/GRAPH_FIELD_BRIDGE_20260922.md`. The actual graph map has now been confirmed by the other session and Anton.

Forskningsrot:

```
../3fold-motion-engine/_private/romi_collab
```

## Verified run state: the entire day queue is delivered

- `build/field_day_delivery_wave_20260922.json`: **57/57 `report_delivered`**, completed 2026-09-22 **19:29:02 +02:00**.
- All 57 `build/<LANE>/RESULTS.md` exist and have the lane name in the first five lines. File hashes and dependencies are in `~/research/field_handover_20260922/LANE_INVENTORY.json`.
- Supervisor PID 538446 has finished; **no own swarm_worker lanes are running** at the process check around 20:15. An OpenCode lane from the other graph session was running; leave it alone.
- Delivery is not the coordinator's scientific approval. Many late reports have not yet been read by root. The automatic queue did not record new results as "holds" and published nothing.
- The first work step is to read the latest audits and raw results and then decide which new hypotheses they actually justify. Do not restart the same finite queue thinking that it keeps thinking or expands itself.

Queue definition file: `lanes/field_day_delivery_wave_20260922.json`.
Historical launch mechanism: `python3 lanes/main_innovation_wave_20260921.py start --wave field_day_delivery`.
This is a reference, not an instruction to restart the completed queue.

If swarm_worker is used again later: `NEW_SESSION=1 lanes/run_lane_v3.sh <LANE> opencode-go/swarm_worker-v4.1-flash`. Read the driver and the brief's limits first. The last permitted shared launch cap was 10 including other sessions. The account key is already configured; never copy or print it. Muse's latest verified attempt gave 403 FreeTierError; reset is unknown.

## Tracks and reports to take over

| Track | Base round | Latest extension |
|---|---|---|
| P1 contact graph, mathematics, sliding regime | U274/A190, U279/A195 | U289/A205; U281/A197 hydrodynamics/transfer |
| P2/P4 events and sensitivity | U275/A191 | U280/A196 |
| P3 exact spatial/temporal queries | U276/A192 | U283/A199 |
| P4 support/load regions | U277/A193 | U284/A200 |
| Local elasticity and actual adaptive DOF | U262/A178, U264/A180, U265/A181, U269/A185 | U285/A201 |
| Digital topology, thin materials and conservation | U266/A182, U278/A194, U272/A188 | U288/A204 |
| Geometric representation and serializer | U267/A183 | Related to U288/A204 |
| Flow/wall geometry | U268/A184 | U282/A198 |
| OED and practical geometry measurement | U256–U260/A172–A176; U270/A186 | U286/A202 |
| Sand/DEM/MPM and numerical position representation | U273/A189 | U287/A203 |
| Graph/field: geometric alternatives and shared evidence | U290/A206, U291/A207 | U292/A208, query-preserving world compression |
| Possible product delivery | U271 | A187 |

`build/NEXT_FIELD_BRIEFS_20260922/` contains four older draft briefs for trace refinement, OED dictionary, local sand coordinates and thin region topology. **They overlap U285–U288 now executed.** Read the outcome first and do not give the same experiments new lane IDs because the drafts remain.

### Earlier root reviews that must retain their scope

- A191: support for event/reset sensitivity in the declared 2D model for initial direction. Not a general material-parameter derivative; final-time kinks and explicit material dependence gave errors.
- A195: standard Gram identity and obstruction in the reviewed translational RPY cases. A new matrix-free advantage is not supported. The Darcy row in U279 had a sign error and is not accepted evidence.
- A188: type identity, conversion and checksum fixes reproduced. Registry/provenance, actual UNKNOWN and preservation of thin regions failed. No approved geometry integration.
- A183: encoding fix and finite enclosure cases reproduced. The exact-expression baseline won on storage. An observed deficit of 4.24e-15 mm is not a universal outward-rounding guarantee.
- A186: U270's negative full-cost result was confirmed; planning cost corrected to 33120. The shape gate was tautological and some precision/coverage numbers were not reproduced.
- A194: total volume can be preserved while the material partition is wrong; crossing-slab gave +3200 mm³ in material 2. Handmade FEM fractions did not substantiate the whole chain.
- U262/A178: actual DOF change and correct global Schur comparison, but the coarse model already met the original error requirement. No supported full-cost gain for that requirement. U264/A180's paid selector was negative.

See the respective audit and HUNT for exact conditions. Do not extend conclusions from these summaries.

## Product and publishing state

P1–P4 and the reproduction packages are already published. Public motion-HEAD at the check: `aae29ed` (clear direct links to four PDFs), previous release `18b91a5`. The old task "finish papers by tomorrow" has passed; do not reopen manuscripts as if they were ongoing referee rounds.

The field engine's three older fixes are also already published:

- staging `../3fold-field-engine`: `597600c3031b1af29e4bdd0351a22ced61baee4d`.
- public `local_path`: `3ffdb0edb119fdfb3484f4d0aa36639412f51508`.
- Publika commits: `2b91296`, `9705033`, `3ffdb0e` (stale owner-celler, massgrid/tre frames, deterministiska owner-ties). Staging var rent vid ny kontroll.

**U271/A187 explicitly reports no further release candidate.** A187's report has been read for this handoff; its tests have not been rerun by root. The chosen payload is empty and `integration_ready=false`. The audit reports 30 regression passes and confirms already published files. This does not mean a new field release is ready.

U273/A189's sand patch is a separate track. Read A189 and later U287/A203 before deciding. U271's local installer accepts `--apply` but lacks an actual write path according to A187; do not use it as a working installer without correction and review.

No new commissions or pushar were made during this handover. Communicate concrete new deliveries with Anton before publishing. No private project names or model attributes in public filer/commits. Staging Author `Anton <[epost borttagen]>`; publikt `Anton Björkegren <antonbj3@users.noreply.github.com>` enligt tidigare instruktion.

Previously authorized follow-up: 40/40 SMTP-accepted, Michael Posa excluded. Old email queues paused. **No further mailing, no arXiv upload.** Anton writes his emails and replies himself.

## Infrastructure still running: the dental move

Isaac Sim assets (~175 GiB) were removed on Anton's explicit instruction. The program and Collected_Robots were preserved; recovery instructions/inventory in `~/research/storage_cleanup_20260922/`.

The dental move's supervisor **PID 654852** is still running. The original **52 005 files, 113 330 792 472 bytes** have finished copying. At the latest check, the step was **`verifying_checksums`**, rsync PID 686331. Live status:

```
external_research_path
```

Source `external_media` → destination `external_mount`.
Script `move_dental_corpus.py` performs a full rsync checksum check, checks that the source is unchanged, replaces the original address with a symlink and deletes the verified old copy **only after checks pass**. A full disk at the original source was the reason for the move.

Let the process finish. No duplicate move or manual deletion based on copying percentage. On error: read supervisor/state/verify logs and preserve source/backup. The script cannot simply be restarted when the destination already exists. Old dataset harvesters and watchdog were stopped because they would otherwise restart when space becomes available; do not restart automatically.

## Graph merge: clearly paused, correct source map now known

Anton wants two dedicated private projects: general BodyTwin, and dental = BodyTwin + manufacturing. The field engine must be a geometry/material foundation for both. The other graph session has provided the correct source map:

**`local_path`**.

- CS-I `data/ANCHOR_GRAPH.json`: **97 nodes, 14 goals, 81 dependencies**. The small node count is correct for the dependency layer; the graph has extensive measured data.
- CS-I's three `CONSTRAINT_NET_{VEHICLE,PROJECTOR,FLEET_FACTORY}.json` and `CONSTRAINT_STRESS_MAP_V1.json` are the variables' quantified relations and stress points.
- `FOLD_LEDGER*.jsonl` in relevant worktrees I/L/F/H/G/B/C/J: about 2770 folds according to the map. Verify the actual disk, deduplicate copies with provenance, bind exact identities; free-text consumers that do not match remain unbound evidence.
- BodyTwin's own biological graph: `source_repository/data/MECHANISM_ANCHOR_GRAPH.json`, **3950 nodes**; 3873 OPEN, 66 ASSUMED, 7 REFUTED, 3 DEFERRED, 1 PROVEN. Source registry alongside. It is not classified as personal data; keep the entire workspace private.
- BodyTwin's inherited CS mirror of 66 nodes is older than CS-I and must not be the manufacturing source.
- **Seed_graph's 43751 nodes are operator text→commits; substrate_graph is machine telemetry. They are not knowledge graphs and must not be added to scientific node counts.** CS-D's older GRAPH_STORE/FEDERATED_GRAPH are superseded generations.

Here root made a context error: first built a merged anchor graph and described readiness too broadly; then moved too quickly to a large seed-index file when Anton questioned the counts. Anton and the other session corrected this. Preserve that correction.

Preliminary workspaces exist in **`local_path,dental}`**. `SETUP_STATUS.json` says `paused_for_session_handoff`, source map now confirmed, full import **not complete**. The original sources are untouched.

Built: lossless namespaced anchor-graph reader, frozen source copies, status preservation, separate coupling notes, registries, profiles through existing `EngineProfile.merge`, version generation and manifest checks. 3950/4047 are the anchor-graph views. Six targeted tests passed. An independent read check confirmed round-trip and no source drift. The new engine is pinned to the graph session's research commit `73e76dd83a601ddf2ccb5bab041512dd426ddb2c`, separate from canonical public `3928720` in `~/projects/graph_workspace/pub`.

Remaining: a shared contract for variables/relation edges/measurements, all measurement logs, stress map, margin/σ/units/validity, explicit source bindings and reviewed cross-domain aliases. Missing σ must not be replaced with a default to make something run; shared evidence must not be counted twice. Contradictions and knowledge gaps are kept open. Original file formats are preserved behind the readers.

Detail: `tools/workspace.py` has late additions for `CONSTRAINT_NETS.json` and an explicitly unresolved relation endpoint not yet in the active generations. No final `refresh` was run after these additions. START/README and active content must be reconciled when merge work resumes. `inherited_evidence_audit_ok=false` is deliberately reported: older BodyTwin has 4793 evidence references that the engine's validator does not recognize and three decided nodes without evidence; not 4793 proven false claims.

## BodyTwin/the collaborator: build on verified context

`~/research/COLLABORATOR_BODYTWIN_CONTEXT_2026-09-22.md` contains read emails, primary sources and actual code review. the external collaborator/AAU has worked for a long time with the reference model and patient-specific geometry/morphing; do not present such capacity as if it were missing in his world. His email asks for a concrete demonstration, physiological geometry/parameterization and a geometry manager.

Current published BodyTwin code: `external_mount`, commit `b95a8dc573d0beaa0d6bc8eca0ed46c0f3ce1658`, compared against remote. Contains strict mesh seam, voxel/tet material boundaries, rigid mm frame, region mass/moment and explicit material-region registration. An actual CPU mesh→SDF bridge exists in `examples/anatomy/compose_mesh_to_field.py`; a kidney mesh's raster-volume error of 0.0975% is not surface precision or clinical validation. No general anatomical shape model/landmark morphing or biomechanical contact chain was found. 110 CPU tests passed, one was skipped in a separate geometry audit.


## Working rules that must not be lost

- Coordinate producer → separate audit → your own source/result reading before "holds". Audit-PASS always means something bounded, sometimes a confirmed negative outcome.
- New experiments in your own build directories. Preserve locked negative results, baselines, costs, coverage and UNKNOWN; do not change error requirements afterward to make a method positive.
- Recursive seed analysis: idea → mechanism → equation → operation → representation → load-bearing assumptions/physical foundation. Every branch ends in a defined building block, verified knowledge or an explicit gap.
- For HUNT bookkeeping: one row before `## Closed`, assert exactly one occurrence. The latest check found one.
- GPU is shared, RTX 5070 12 GB. Read current locks/protocols; no long runs or `nvidia-smi -q`. CPU thread cap 4 in briefs. Modal supports justified experiments; previous balance estimated at ~$10.27 on the active profile, not a current guarantee.
- Use `ps -eo pid,args`, never `pgrep -f`. Code that stops processes is written to a script file first and matches exact process identities.
- Communicate briefly in Swedish. Anton wants broad context before deep choices, clear status reports and no surprise pushes or mailings.
- Save tokens: read relevant reports, not the same images/PDF pages again when unchanged bytes have already been reviewed. An external queue does not wake the coordinator chat automatically.

## Suitable first concrete action for the taking-over session

1. Check the dental move's live state; let completed verification finish.
2. Read U271/A187 to avoid promising a new product release that does not exist.
3. Read the latest extension audits A196–A208 with strong baselines and failures, starting for example with the graph/field bridge A206–A208 and materials/topology A204.
4. Record what you verified yourself and formulate the next bounded, falsifiable experiment from the results. Follow Anton's new choice of worker and responsibilities; do not routinely restart old jobs.
