# CX-PACKETFACTORY — a deterministic packet factory that keeps 54+ The swarm/reserve_worker slots full without LLM cost per packet

Problem: the local queue (16) and the OVH queue (30) consume ~250 packets in ~15 min. LLM-built packets (CX-SWARMGEN*) cost lane runner budget and do not keep up. Anton: maximise The swarm A/B/C and use up reserve_worker (C) during the eight-hour window (`external_research_path`). The plan wants 25 % independent testing.

## Build `tasks/packetfactory.py` (deterministic, no LLM)
Source: all `results/<id>/` with RESULTS.md and results.json, where the data/code needed for recomputation is INSIDE the result directory or its `inputs/`. Copy exactly what is needed. Size ≤ 50 MB per packet.

Packet families (id prefix = family, e.g. `BT-R-<src>`, `BT-P-<src>`, `BT-X-<src>`, `BT-AG-<n>`):
- **R, reproduction:** recompute the result's key numbers (fields from results.json, max 3) with your OWN code from the raw data in inputs/. Deviation > 1 % relative = flag. The swarm.
- **P, robustness:** perturb one choice that the result depended on (threshold ±20 %, other seed, leave out one person/specimen, other normalisation). Does the conclusion hold? The swarm, or reserve_worker for heavier ones.
- **X, transfer:** the same analysis on another dataset that is already in the result's inputs, or in another result's inputs with the same format (e.g. Imperial ↔ VSD ↔ Keast; JW ↔ DM ↔ SC ↔ PS). reserve_worker C.
- **AG, audit group:** 3 unaudited register rows per packet (from `notes/RESULTS_INDEX.md`, "ej granskad"), with the sources' own files copied. reserve_worker C or The swarm.

Every packet:
- BRIEF.md in the preamble style: inputs/NIGHT_PREAMBLE.md, "Builds on: <src>", the question, a criterion with numbers, the counter-test, and PREREG first;
- `inputs/DATA_SUFFICIENCY.md` plus a load test;
- run `results/CX-SWARMGEN5/check_packets.py`, or the equivalent. Only packets that pass are queued.

Exclude:
- sources whose inputs contain private health data or credentials;
- results that already have an R/P packet;
- CX-* results larger than 50 MB.

Spread profiles A/B/C evenly for The swarm. reserve_worker ALWAYS goes to C.

## Refill: `tasks/packetfactory_loop.sh`
- Every 3 min: count unfinished queue lines. If < 80, generate the next 150 packets (in priority order: AG for unaudited rows → R for the newest results → P → X) and append to `tasks/lanes/bt_queue.txt` with >>.
- Stop when the sources run out or at 23:40.
- Run it as `systemd-run --user --collect --unit=bt-packetfactory --slice=jobs-bodytwin.slice -p MemoryMax=1G ...`.
- Log to `tasks/lanes/packetfactory.log`.

## Deliverables
- The factory, the loop (started), a test that generates 5 packets of each family, and check_packets with 100 % pass.
- `results/CX-PACKETFACTORY/RESULTS.md` starting with `# CX-PACKETFACTORY`: the number of available sources per family, estimated packets in total, and whether the loop is running.

Rules: lane runner has full permissions in the workspace. `~/projects/bodytwin` is read-only (copies are allowed as packet inputs). restricted model data/the collaborator data may go into LOCAL/OVH packets (our machines) but never to Modal/coordinator cloud.
