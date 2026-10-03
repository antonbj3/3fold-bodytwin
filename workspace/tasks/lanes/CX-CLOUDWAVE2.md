# CX-CLOUDWAVE2 — collect coordinator cloud wave 1 and keep the pace so that ALL ~250 USD gets used during the window

State: 10 BodyTwin sessions in `tasks/cloud/receipts.jsonl`, balance 225.1 USD (`tasks/cloud/credit_log.jsonl`). Field runs its own (~4–6/h). The pace is too low: the target is ≈ 31–35 USD/h in total until ~23:30. Anton: "use up all 250".
Read: `tasks/lanes/CX-CLOUDLANES.md` (incl. the ADDENDUM: launch from the trusted `lane_repo`, one branch per lane, sequentially; no trust edits), `results/CX-CLOUDLANES/RESULTS.md`, and `tasks/cloud/*.py`.

## Tasks
1. Collection: run collect_cloud.py for all BodyTwin sessions. Git return: `git fetch` in lane_repo and fetch the result branches into `results/CLOUD-<ID>/`; use the event path as a fallback. Record per session: model_id, status, the cost (the balance delta), and whether RESULTS arrived.
2. Evaluation: for each returned session, two lines — the key result, and whether it holds up by our own recomputation where that is cheap. Report plainly if the model has invented data or sources.
3. Next wave: 15 new BodyTwin cloud tasks. Hard, bounded work where coordinator gives the most value. Only public data or our own code, no restricted model data/the collaborator/TLEM. Priorities:
   - follow-ups on what wave 1 found;
   - an independent audit of our strongest claims (A359 T1, A369 what-if curves as a mathematical/numerical question with our own code without internal data, A363 batch solver as an algorithm);
   - derivations that unlock BT-2: the patella mechanism, operating-range identifiability;
   - B1 (the minimal measurement protocol), B5 (disc pressure), and the NATO chain;
   - the hyperrealism roadmap (`results/HYPERREALISM_ROADMAP_20260924/`), with one or two bounded cell/tissue questions that have a public reference.
4. Pace control: after launching, measure the cost per session. Start more at ~30-min intervals, so that the BodyTwin share plus Field's reaches ~0 USD at ~23:30. Write the plan to `tasks/cloud/PACING.md` (balance, sessions/h, the expected end). If the lane runner round ends, leave a script `tasks/cloud/next_wave.sh` that the coordinator can run, which launches the next N prepared lanes.
5. `results/CX-CLOUDWAVE2/RESULTS.md` starting with `# CX-CLOUDWAVE2`.

Rules: never print or copy the OAuth token. Only the promotional credit; stop if paid usage appears. No PRs/merges/public pushes (branches in the private research-cloud-lanes are OK).
