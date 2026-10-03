# CX-PLANNER — step on the gas: reserve_worker packets (account C reset), more Sol lanes, and optimal tasks for progress in all BodyTwin domains

Anton (24/9 ~15:10): "reserve_worker on C is reset, and we have Space The swarm. Step on the gas properly, 30 agents or so with optimal tasks for progress in all our domains." The coordinator (coordinator) is saving tokens, so you plan and launch.

## Read first (the whole picture)
- `notes/RESULTS_INDEX.md`: A1–A500; especially the recent A306 onward.
- `tasks/NIGHT_PREAMBLE.md` (§0 direction: innovate across the whole field, test against established baselines).
- `results/CX-PARAM/PARAM_MAP.md` (P1–P9).
- `results/CX-KNEENET/KNEE_CONSTRAINT_NET.json`.
- The graph `~/projects/bodytwin/data/MECHANISM_ANCHOR_GRAPH.json`: run `python3 ~/projects/bodytwin/scripts/anchor_graph_tools.py priority` / `next-actions` if the tool exists, otherwise grep status OPEN + high risk/value.
- `tasks/lanes/*.md` (existing lanes; do not redo them).

## Domains to cover (spread effort; not only the knee)
1. Geometry manager/shape (axes 1, 2, 9).
2. Muscle parameters/strength (T1 held A359; T2 running; the patella arm A367).
3. Joint force/null models/constraint nets (A324–A358).
4. Spine/L5-S1 and the upper body (the axis-6 the reference model VBR finding A188; N2b's shoulder band).
5. Speed/scale/emulator/population (A283, A303, A357, A363).
6. Tissue/field/material/injury risk (axes 4, 10, 12, 15).
7. Sparse sensors/video/IMU (axes 7, 11; A321 three sensors −42 %).
8. Validation/certificates/traceability (axis 5).
9. Parametric what-if/surgery planning (CX-WHATIF2 running).
10. Transfer to/from dental, field and graph (read-only in their workspaces; A310 proposals).

## Deliverables
1. **30 reserve_worker packets** `results/BT-D001` … `BT-D030`: harder than the The swarm packets. Each has a bounded but substantial question with a number criterion, a null model/counter-test, and a complete `inputs/` with `DATA_SUFFICIENCY.md` and a load test (same check as CX-SWARMGEN5 — reuse its check_packets.py). Spread at least 2 per domain. Append `C reserve_worker BT-D0xx` to `tasks/lanes/bt_queue.txt` with >> only after the check passes.
2. **6 new Sol lanes** for the highest-value leads that need heavy building (not The swarm-sized). For each, write `tasks/lanes/CX-<NAME>.md` in the same style as e.g. `tasks/lanes/CX-WHATIF2.md`: build on, tasks, PREREG, criteria, counter-test, resources. Launch each with exactly:
   `systemd-run --user --collect --unit=bt-lane_runner-CX-<NAME> --slice=jobs-bodytwin.slice -p MemoryMax=4G -p CPUQuota=200% -p RuntimeMaxSec=28800 --working-directory=. --setenv=PATH=~/.npm-global/bin:$PATH --setenv=HOME=$HOME env MAX_ROUNDS=3 /bin/bash tasks/lanes/run_codex.sh CX-<NAME>`
   Prioritise leads that have NOT been tried. Examples: GC fluoroscopy as the CT↔marker link (A368); L5/S1/spine against VBR with an individual lever arm; IMU→joint force with uncertainty; the certificate stack of the geometry manager as a product; a population product on GPU with the N2b solver.
3. `results/CX-PLANNER/RESULTS.md` starting with `# CX-PLANNER`: a table of domain → reserve_worker packets → Sol lanes → why these are the optimal next steps (value per compute).

## Rules
- lane runner has full permissions, but `~/projects/bodytwin` and other sessions' workspaces are read-only. No emails, pushes or publishing. Do not touch credentials.
- At most 6 of your own Sol lanes. Heavy computation goes to OVH (`tasks/cloud_run.sh`, BodyTwin ≤ 12 vCPU, finish by 22:45) or Modal (public data only).
- Internal restricted model data/the collaborator/LHDL data never goes to Modal.
