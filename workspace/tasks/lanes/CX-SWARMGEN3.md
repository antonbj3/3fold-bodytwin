# CX-SWARMGEN3 — 40 new self-contained Space The swarm packets (round 3) with a STRICTER data check

Same format and rules as `tasks/lanes/CX-SWARMGEN.md` (read it), with these changes:

1. Lesson from round 1: B150 (Keast meshes without fibula), B157 (VSD meshes missing) and B159 (no personal skin/bone pairs) passed check_packets but lacked the data the QUESTION needs. The check only verified files named in the brief. New check: for each packet, write `inputs/DATA_SUFFICIENCY.md` stating which data the analysis requires, which file carries it, and a small load/shape test (e.g. mesh has ≥ 2 components when fibula is needed; arrays have the right length). Run the test in check_packets. A packet that fails is NOT queued.
2. Topics (priority order):
   a0. NEW LEADS (priority): A359 DXA leg lean beats the strength-scaling baseline laws (NHANES) — small follow-ups on sex/age/BMI residuals and the exponent; A361 TLEM muscles 2.5–4× too strong vs JW isometric strength — which muscle groups, which angles; A358 JW feasible set wrong; A360 muscle-CT pilot. Also follow-up questions on the latest The swarm rows (A328–A356).
   a. Redo B150/B157/B159 with complete data, or give an alternative definition that works with the data we have.
   b. B151: the time axes between markers and GRF appear shifted by ~80 ms. Check the sync across all persons/trials (cross-correlation of the heel marker's vertical velocity against the GRF onset).
   c. Follow-up questions on A324–A327 (N1g, constraint net, early stance, muscle parameters) that are small and bounded. Examples: knee angle as a term per activity; the activity×GRF interaction per person; the sign convention of the ID moment against OpenSim ID for the DM trials; the TLEM slack-length window on straight-line paths per element.
   d. Follow-up questions on the The swarm results booked in `notes/RESULTS_INDEX.md` from A313 onward.
   e. New angles on the board's axes 1–17 that have not been taken, per the preamble §0 direction (innovate across the board, test against established baselines).
   f. 10 audit packets BT-AN-G20…G29 for the unaudited rows from A328 onward (The swarm/reserve_worker rows; each row once) (not the CX-* rows; they are audited separately).
3. IDs: BT-B210…BT-B249 and BT-AN-G24…G33. Build the packets DIRECTLY in `results/<id>/` (not in your own directory), and append the queue lines to `tasks/lanes/bt_queue.txt` with >> once the check passes (spread over A/B/C; audits go to C). If an id directory already exists, skip that id.
4. `results/CX-SWARMGEN3/RESULTS.md` starting with `# CX-SWARMGEN3`: table id → question → axis → size → data check.
