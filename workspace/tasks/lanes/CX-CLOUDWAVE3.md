# CX-CLOUDWAVE3 — own BodyTwin's coordinator cloud pace for the rest of the window (to ~23:15): prepare, launch, collect, evaluate

Predecessor: CX-CLOUDWAVE2 (done). Read its `results/CX-CLOUDWAVE2/RESULTS.md`, `WAVE1_REVIEW.md`, `WAVE2_REVIEW.md`, `tasks/cloud/PACING.md`, `tasks/cloud/next_wave.sh`, and `tasks/lanes/CX-CLOUDLANES.md` (ADDENDUM: the trusted lane_repo, one branch per lane, sequentially, no trust edits).

## Target
Anton: **use ALL ~250 USD**. At ~17:50 local time, ~162 USD remain; the window ends ~23:30. The BodyTwin share is ≈ 10 sessions/h at ~1.7 USD/session. Field runs 6–8/h. Keep the joint pace ≈ 27–30 USD/h. Read the balance ONLY at :00 and :30 (Field reads at :15/:45; 429 otherwise).

## Work loop (repeat until ~23:15, then collect for the last time)
1. Prepare: have ≥ 15 unstarted, preregistered lanes ready at all times (`tasks/cloud/lanes/CLOUD-W4-*` …), hard and bounded, where coordinator gives the most value. Only public data or our own code (no restricted model data/the collaborator/TLEM/LHDL). Priorities:
   - (a) follow-up questions from the returned sessions;
   - (b) an independent audit of our most important claims that has not yet been done: A907 D1 parity as a numerical question, the A363 algorithm, A295 exactness, A303 emulator bounds, A194 moment arms with public geometry;
   - (c) BT-2 knee: patella mechanics/operating range with public models and DOI data;
   - (d) the plan's brainstorm B1–B5 and the hyperrealism roadmap (`results/HYPERREALISM_ROADMAP_20260924/`): 1–2 bounded cell/tissue/organ questions per hour with a public reference;
   - (e) theory questions that unlock things: when is contact from measured kinematics well-posed (coordinate with Field's INVIVO-THEORY, do not duplicate), and identifiability of muscle parameters from strength curves.
2. Launch per 30-min block: ~5 lanes (`next_wave.sh 5`), adjusted to the balance delta.
3. Collect: `collect_cloud.py` (paginated), then fetch Git branches from lane_repo into `results/CLOUD-*`.
4. Evaluate: per returned session, write two lines in `results/CX-CLOUDWAVE3/WAVE_REVIEW.md` — the finding, and your own check where it is cheap. Report plainly any invented data or sources. If a finding corrects our register: add an inline `[CORRECTED …]` in the row (as with A359).
5. `tasks/cloud/PACING.md`: a new section per block.

`results/CX-CLOUDWAVE3/RESULTS.md` starting with `# CX-CLOUDWAVE3`: an overview at the end (sessions, cost, the most important findings, corrections).

Rules: never print or copy the token. Only the promotional credit; stop if paid usage appears. No PRs/merges/public pushes.
