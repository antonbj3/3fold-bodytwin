# CX-SWARMGEN — build 40 self-contained Space The swarm job packets from open gaps (BodyTwin against published baselines)

Goal: Space The swarm (free, accounts A/B/C) must always have jobs. The queue is empty. You write packets; the queue runs them.

## Sources for gaps (read these first)
- `notes/RESULTS_INDEX.md` (A1–A305): rows with UNKNOWN, "packaging deficiency", "not reviewed", "REFUTED/fails", "next step".
- `tasks/NIGHT_HUNT_BACKLOG.md`.
- `tasks/NIGHT_PREAMBLE.md` (rules for all jobs; copy it into every packet).
- Existing packet format: `results/BT-B143/BRIEF.md` + `results/BT-B143/inputs/`.
- Before every topic: grep `~/projects/bodytwin/data/MECHANISM_ANCHOR_GRAPH.json` and `~/projects/bodytwin/scripts/msk/` so that no packet redoes existing work. Name the node/file in the brief under "Builds on".

## Deliverables
1. 30 research packets `results/BT-B150` … `BT-B179`. Each is a bounded, falsifiable question with a number as the criterion, a null model where relevant (preamble §2: B24 and N1/N1g for joint force), and local computation ≤ 60 s, 1 thread, < 1 GB. Priorities:
   a. Rerun packaging errors with complete data: B125 (Keast meshes), B133 (synced marker/GRF arrays), XF6 (H3 npz + solver state), and others you find.
   b. Follow-up questions to A266–A305 that can move an axis (recruitment vs EMG, knee force vs N1g, geometry manager, speed/emulator, HJC, soft tissue).
   c. New angles against the reference model on axes 1–17 that nobody has taken yet.
2. 10 audit packets `results/BT-AN-G10` … `BT-AN-G19`, 3–4 unaudited rows each from RESULTS_INDEX (marked "not reviewed"). The audit recomputes the key numbers from the result's own files (copy them into inputs/) and gives HOLDS / WITH RESERVATIONS / FAILS per row.
3. For each packet, BRIEF.md in B143's compact format (Swedish, ≤ 25 lines): title line with the id, a line pointing to inputs/NIGHT_PREAMBLE.md, "Builds on", the question, the criterion with numbers, the null model, and deliverables (PREREG.md + PREREG.sha256 first, results.json, RESULTS.md starting with the id).
4. `inputs/` with EVERYTHING the job needs: code, arrays, meshes, the previous job's results.json. The model cannot read outside the packet. At most 200 MB per packet; if data is larger, make a sufficient excerpt and say so in the brief.
5. Validation script `results/CX-SWARMGEN/check_packets.py`: for each packet, every file path mentioned in BRIEF.md exists under the packet, inputs/ is non-empty, there is no absolute path outside the packet in the brief, and the size is ≤ 200 MB. Run it; only packets that pass go on.
6. Only packets that pass: append one line per packet to `tasks/lanes/bt_queue.txt` in the format `PROFIL swarm JOBB_ID`, spreading profiles evenly over A, B and C. Audits get `C swarm`. Append with `>>`; never rewrite the file.
7. `results/CX-SWARMGEN/RESULTS.md` starting with `# CX-SWARMGEN`: table id → question → axis → inputs size → check status.

## Forbidden
- Material from other private domains (tissuetwin, glioma, excluded_category, dental, personal health data). Only musculoskeletal BodyTwin data.
- API keys or credentials in packets. Writing in `~/projects/bodytwin` (read-only). `pgrep -f`, broad kill commands, heavy computation locally (> 60 s).
- Redoing a job that already has RESULTS.md under results/.
