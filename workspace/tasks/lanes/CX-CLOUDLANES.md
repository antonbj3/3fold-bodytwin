# CX-CLOUDLANES — start BodyTwin's share of the coordinator cloud credit (~100 USD) now: launcher, collector, first wave

Read:
- `external_research_path` (Anton: use the whole credit during the eight-hour window);
- `external_research_path` (verified state: coordinator 5.5 works; artifacts retrieved via authenticated session events; Git push BLOCKED until Anton selects the repo in coordinator.ai/code);
- `start_smoke.py`, `inspect_cloud.py`, `launch_receipt.json`, `lane_repo/` in the same directory.

## Return path until Git works
The launch bundles the local git checkout (the input arrives). The result is returned via the session's event API, as the smoke test did. Every cloud task must therefore end by printing its result files verbatim in the final message, between markers:
`=====FILE <path>=====` … `=====END FILE=====`
Include RESULTS.md, results.json, and small code files (≤ 200 KB in total). The collector parses them into `results/CLOUD-<ID>/`. When Anton has granted Git access, switch to the branch return and verify it with a fetch.

## Deliverables
1. `tasks/cloud/launch_cloud_lane.py <lane_dir>`:
   - the lane directory is its own small git repo (git init; BRIEF.md + inputs + code; ≤ 20 MB);
   - start with `coordinator --model coordinator --effort high --cloud "<prompt>"` from the lane directory, in a pseudo-TTY if needed (see how start_smoke.py/launch_no_tty handled it);
   - parse the session URL/ID and write `tasks/cloud/receipts.jsonl` (lane, session_id, time, requested model);
   - refuse duplicate starts per lane (lock).
2. `tasks/cloud/collect_cloud.py`:
   - for each started session, read status plus events via the same endpoints as inspect_cloud.py, using the token from the same auth file (never print the token);
   - extract the FILE blocks to `results/CLOUD-<ID>/`;
   - record the model actually used (model_id from the metadata) and status;
   - read the account balance via `/api/oauth/usage` → `tasks/cloud/credit_log.jsonl` (time, balance).
3. First wave: 8 BodyTwin cloud tasks. **Constraint: only public data or our own code in the bundle. NO restricted model data/TLEM/the collaborator's data/LHDL/internal A matrices.** Pick among:
   - (a) independent reproduction of A359, NHANES strength vs DXA (the cloud downloads from CDC itself; its own code);
   - (b) patellar tendon moment arm vs knee angle from public knee models and literature with DOI → a tracking model with its maximum around 45° (A367);
   - (c) L5/S1/disc: intradiscal pressure against Wilke's in vivo data from lever-arm models (B5);
   - (d) a minimal measurement protocol per person (B1): an OED derivation over DXA/X-ray/landmarks/strength/knee angle with the numbers from our register rows (copy the rows, not internal files);
   - (e) the NATO chain activity → tissue load → ε-N fatigue → injury risk with a public ε-N/bone-fatigue reference and an error budget per link;
   - (f) an independent audit of the constraint-net/N1g reasoning (A324–A358) as a mathematical question: why N_eff = 1, and which information must be individual (only register text + public GC data if it can be bundled legally);
   - (g) a vascular–interstitial transient (Starling/lymph) with a public reference (lane runner plan BT-B);
   - (h) the strength–angle curve as an identifiable individual signature (B4): an identifiability derivation with public curves.
   Each gets BRIEF.md with the strongest baseline, prior counterexamples, criteria and deliverables, following the preamble style (PREREG first).
4. Launch the 8, verify that the first one returns via the collector (a real session, model_id coordinator-5-5), and log the balance before and after.
5. `results/CX-CLOUDLANES/RESULTS.md` starting with `# CX-CLOUDLANES`: sessions, the model per session, balance, return status, and exactly what is required for Git return (Anton's step).

## Rules
- Never print or copy the OAuth token or keys to files, logs or bundles. No PRs, merges or public pushes.
- Only the promotional credit. Do not enable paid extra usage; stop if the balance approaches 0 or if usage appears outside the credit.
- Coordinate with Field: BodyTwin ~100 USD, Field ~100 USD, reserve ~50 USD. Write the reservation in the plan's §14 (`external_research_path`) as a bounded addition in BodyTwin's row.
- lane runner has full permissions; `~/projects/bodytwin` is read-only.

## ADDENDUM 24/9 ~16:00 (coordinator) — READ BEFORE CONTINUING
- The first launch worked (CLOUD-A359, session_01KrvT34D6J9bWS2kY37YKxN). CLOUD-A367 got stuck on coordinator Code's interactive "trust this folder" question, because each lane is a new directory.
- Do NOT solve it by editing ~/.coordinator.json or other trust/permission settings, and do not answer the dialog automatically.
- Instead, launch ALL lanes from the already trusted directory `external_research_path`, which belongs to antonbj3/research-cloud-lanes. For each lane:
  1. `git checkout -B lane/<ID> <baseline>` (the baseline = the repo's main/baseline commit);
  2. copy the lane bundle there (BRIEF.md + code/data, ≤ 20 MB);
  3. `git add` + commit;
  4. launch from that directory;
  5. then move on to the next lane, SEQUENTIALLY (the bundle is taken at launch).
  Afterwards, go back to the baseline branch. Adapt launch_cloud_lane.py to take a lane id + a bundle directory and do exactly this.
- The Git return now works (Anton: the repo selection is done). Tell each session to commit and push on its own branch `coordinator/<ID>` or the session branch, and verify with git fetch in lane_repo. The event return remains a fallback.
- Use ALL ~250 USD (see results/CX-CLOUDLANES/UPDATE_FROM_ANTON.md). Field runs its own waves (4–6/h). BodyTwin should start the remaining 7 now and then a new wave after measuring the cost of the first ones.
