# BodyTwin/Field promotional-credit pacing, 2026-09-24

Target window ends about **23:30 Europe/Stockholm = 21:30 UTC**. Only the 250 USD promotional balance may be used; extra usage is disabled (API `extra_usage.is_enabled=false`, `spend.enabled=false`, paid spend 0 at 14:30 UTC). Stop new starts immediately if paid usage appears or the promotional balance is exhausted. Credit belongs to BodyTwin and Field together. Account balance deltas cannot be assigned to one session while jobs overlap; the session API separately exposes `external_metadata.usage.cost_usd` per session.

## Measured state

- 14:22 UTC: 224.168529 USD remained; 10 BodyTwin sessions already complete from wave 1.
- 14:32 UTC: 212.038828 USD remained; 15 new BodyTwin coordinator sessions launched sequentially from the trusted checkout, all running on first collection. Balance decline over this burst: 12.129701 USD in roughly 10 min, shared with Field. This **is not a sustained cost/session estimate**.
- At 14:32, time to 21:30 is 6 h 58 min. Required mean joint burn: about **30.4 USD/h**. The 10-minute burst rate was ~73 USD/h because 15 sessions were active; do not extrapolate it to the rest of the window.
- All 25 BodyTwin session receipts and account-delta brackets are in `results/CX-CLOUDWAVE2/session_ledger.json`. `individual_cost_usd` now records the session API usage snapshot, with an observation timestamp. Each balance bracket remains a **shared account delta**, and brackets overlap.

## Control at each roughly 30-minute checkpoint

1. Run `python tasks/cloud/collect_cloud.py` once; inspect `balance`, `extra_usage`/paid spend if querying usage, and which sessions are still running. Rebuild the ledger with `python tasks/cloud/summarize_receipts.py`. Do not start more while the 15 wave-2 sessions are actively draining at a rate above the target.
2. Compute `required_rate = promotional_remaining / hours_until_21:30_UTC` and `observed_rate = (previous_balance-current_balance)/elapsed_hours` from successful balance observations only. A 429 leaves the old observation stale; wait and retry later.
3. At the first checkpoint near **15:00 UTC**, if observed joint rate is below required and most wave-2 sessions have returned, start 2–4 prepared wave-3 lanes. If observed is above required, wait another 30 minutes. Use the empirical shared decline after session completion to refine batch size. Field's own 4–6 jobs/hour must be included in the observed rate.
4. For a later checkpoint, a starting heuristic is `N = ceil(max(0, required_rate-observed_rate) * 0.5 / c)`, capped at 4 initially, where `c` is the most recent *blended* USD per extra BodyTwin launch over a completed batch. This is only a control heuristic; use the observed session API cost distribution and adjust if Field changes pace. A running session's usage snapshot is not its final bill.
5. Run `tasks/cloud/next_wave.sh N` to start the next N of 12 preregistered reserve tasks. It verifies promotional balance and that paid usage is disabled once before each small sequential batch; an HTTP 429 stops the batch until a later retry. Twelve reserve lanes cover only the next few checkpoints; prepare further bounded public/own-code work if the balance remains high after those launches. Avoid duplicating completed analyses.
6. Near zero, keep paid usage disabled and verify the actual promotional remainder before any additional start. Stop at 0 or a credit lock; do not enable an overage. Log the final balance and any unspent remainder plainly.

The target is a controlled joint balance of approximately zero by 21:30 UTC. A point forecast is not credible until wave 2 completes and Field's concurrent consumption is measured; the first post-launch balance indicates ample short-term pace, and the next launch decision is intentionally deferred to the 15:00 UTC checkpoint.

## 14:35 UTC update

The second collection returned four wave-2 packages and 11 were still running. Shared promotional balance was **203.719577 USD**. The 14:32–14:35 burst consumed 8.319251 USD, so no reserve-lane start is warranted during this active burst. Remaining-time target to 21:30 is about **29.4 USD/h**. Reassess only after the running wave settles and at the planned 15:00 checkpoint.

## 14:42 UTC update

Twelve of fifteen wave-2 results have arrived; `CLOUD-W2-OPRANGE` and `CLOUD-W2-VASCULAR-ID` are idle while their own background computations continue, and `CLOUD-W2-CELL-GEOM` is still running. The latest promotional balance is **193.294824 USD**. The 14:39–14:42 shared decline is about 2.33 USD over ~2.7 min, still above the ~28.1 USD/h rate needed over the remaining window. No reserve start yet.

## 14:51 UTC update

Fourteen of fifteen wave-2 result packages arrived; operating-range task is idle while its background grid is monitored. Promotional balance **189.743942 USD**. Recent shared decline from 14:48 to 14:51 is 1.061738 USD, roughly 21 USD/h, below the about 27.7 USD/h remaining-time target. If this lower rate persists to the ~14:55–15:00 checkpoint and no jobs resume heavy compute, start 2–3 source-anchored reserve lanes with `next_wave.sh`.

## 14:54 UTC reserve start

The first source-anchored reserve lane (`CLOUD-W3-PATELLA-CONTACT`) launched. A second usage preflight hit HTTP 429, so the planned batch stopped after one launch. `next_wave.sh` now checks usage once per small batch and stops closed on 429. Resume with the next two source-anchored lanes after a successful fresh usage check; do not count the failed preflight as a submitted session.

## 14:57 UTC update

One wave-3 source-anchored patella lane is running; wave-2 operating-range has no returned files yet. Promotional balance **188.082002 USD**. Shared 14:54–14:57 decline is 1.183355 USD, about 24 USD/h, modestly below the ~28 USD/h required rate. The remaining two source-anchored reserve starts can be considered at the ~15:00 checkpoint after a fresh successful usage preflight.

## 15:00 UTC checkpoint

Promotional balance **186.646249 USD**. Hours remaining to 21:30 UTC: 6.5; required joint rate **28.72 USD/h**. The latest 14:57–15:00 shared decline is 1.435753 USD over ~3 min, about **28.7 USD/h**, close to target, while one reserve lane is running. Because the usage endpoint recently returned 429, defer any next batch a few minutes and require a fresh successful promotional/paid-usage preflight. Then start at most two if the rate has slipped below target; otherwise wait for the next 30-minute checkpoint.

## 15:04 UTC reserve batch

A fresh usage preflight succeeded with **185.80 USD** promotional balance and paid extra usage disabled. `next_wave.sh 2` then launched `CLOUD-W3-DISC-IDENT` and `CLOUD-W3-CELL-MESH` sequentially; with the earlier `CLOUD-W3-PATELLA-CONTACT`, three source-anchored reserve lanes are now running. Do not launch another reserve batch until the roughly 15:30 UTC checkpoint unless the measured joint rate falls substantially below target.

## 15:16 UTC per-session cost refinement

The session API's `external_metadata.usage.cost_usd` returned values for all 28 BodyTwin receipts. The sum at this snapshot is **47.855973 USD**; some sessions are still active, so it is a lower/current snapshot, not the final cost. The ten completed first-wave sessions total about **15.95 USD** and the completed wave-2 jobs are mostly **1–3 USD each**. Shared account usage at 15:15 was 69.597022 USD, leaving roughly 21.74 USD of account usage outside these BodyTwin session snapshots, plausibly Field and timing differences. This is an accounting inference, not an attribution audit. For pacing, use a working BodyTwin cost of **~1.8 USD per new session** until newer completed batches refine it. With Field's observed contribution in this window, roughly **3–4 BodyTwin starts per 30 minutes** may be needed when the shared rate falls below the required ~29 USD/h. Nine reserve tasks remain after the three already started; further qualified public/own-code tasks must be prepared for later hours.

## 15:22 UTC update

Promotional balance **179.392623 USD**. All three source-anchored wave-3 packages have returned; operating-range wave-2 remains active after a bounded-stop follow-up. Shared credit declined only 0.213849 USD since 15:18 while completed jobs settled, well below the ~29 USD/h remaining-time target. At the 15:30 checkpoint, launch approximately **4** reserve lanes if a fresh credit/paid-usage preflight succeeds. Nine unstarted reserve lanes remain; this would leave five for the following checkpoint. Prepare more qualified lanes before that buffer runs out.

## 15:30 UTC checkpoint

Seven wave-3 sessions have launched in total; three remain active (NHANES nonlinear, NHANES calibration, KKT degeneration). The promotional balance is **168.940566 USD**. Since 15:26, the shared account declined 6.178989 USD in about four minutes while the batch ramped up; that short burst is well above the approximately **28.2 USD/h** needed to 21:30 UTC. Start **no additional lanes now**. Five prepared reserve lanes remain. Reassess near 16:00 UTC using the completed-batch session costs and Field's contribution.

## 15:40 UTC handoff

All **32** BodyTwin sessions have returned `RESULTS.md` and `results.json` via events, and all show `coordinator-coordinator-5-5` in session metadata. The latest successful shared promotional balance is **162.058221 USD**, so the remaining-time target to 21:30 UTC (5 h 50 min) is **27.78 USD/h**. The shared account declined from 168.940566 USD at 15:30 to 162.058221 USD near 15:40, a short-window rate of roughly **41 USD/h**; **do not start more now**. Check again near **16:00 UTC**. The seven completed wave-3 sessions cost **12.537707 USD** by session API; all 32 BodyTwin sessions cost **56.734401 USD**. Only **five** prepared reserves remain, so prepare further distinct bounded public/own-code work before the next reserve buffer is exhausted.

A final extra-usage query returned HTTP 429. The last successful paid-use preflight at about 15:24 reported `extra_usage.is_enabled=false`, `spend.enabled=false`, paid spend 0. No paid mode was enabled by this work. `next_wave.sh` will refuse a new batch until a fresh promotional/paid-use preflight succeeds. After a 429, wait rather than reusing the last balance.

## 17:56 CEST / 15:56 UTC — CX-CLOUDWAVE3 preparation

The handoff balance at 15:39:53 UTC was 162.058221 USD. No further usage read was made before the :00 checkpoint. Seven W3 lanes were already launched, leaving five W3 reserves. This turn preregistered 22 W4 and 25 W5 distinct public/own-code lanes, so 52 unstarted lanes are available. A Krevolin DOI typo in the unlaunched W4 patella brief was corrected against PubMed (correct DOI 10.1016/j.jbiomech.2003.09.010). The collector now supports `--no-usage` so events can be fetched between shared balance slots. Plan: check the promotional and paid-use state at 16:00 UTC; launch approximately five if preflight passes, then collect/review while preserving 30-minute pacing.

## 18:00 CEST / 16:00 UTC checkpoint

The `next_wave.sh 5` preflight succeeded with promotional balance approximately **153.07 USD** (CLI rounded to cents), paid-use gates clear. This was a 8.99 USD shared decline since 15:39:53 UTC, around 26.8 USD/h; the exact balance was not persisted by the old preflight, so subsequent preflights now append their exact observations to `credit_log.jsonl`. Five remaining W3 lanes launched sequentially and have valid receipts: BATCH-CERT, CELL-COUPLE, NATO-ERROR, OED-CORREL, VASCULAR-SAMPLE. 47 W4/W5 preregistered unstarted lanes remain. No further starts until the 16:30 UTC checkpoint unless an exceptional need becomes evident. Event collection uses `--no-usage` between slots.

## 18:09 CEST / 16:09 UTC collection

Two of the five 18:00 W3 sessions have returned through paginated events. OED-CORREL and NATO-ERROR were independently rerun; both synthetic frozen gates pass, with empirical validity still unknown. A read-only `git fetch origin --prune` in the trusted lane repo again shows only `origin/main`, so events remain the result return path. Three sessions continue. No additional credit read has been made off the :00/:30 slot.

## 18:14 CEST / 16:14 UTC session costs

Four returned 18:00 BodyTwin sessions have current API `usage.cost_usd` of 1.490903, 1.279724, 1.033045 and 0.845281 USD (sum 4.648952 USD, mean 1.162238); CELL-COUPLE is still running and has no settled cost. This is lower than the earlier ~1.7 USD/session planning value. The :30 shared balance must govern the next start count, because Field and settlement timing overlap. Full API snapshot is `results/CX-CLOUDWAVE3/session_costs.json`.

## 18:28 CEST / 16:28 UTC bounded cell follow-up

CELL-COUPLE became idle without result files; its last assistant event said a Picard-iteration full run continued in the background. A follow-up was sent to the same cloud session, asking for one check, at most five further minutes, and an explicit partial/UNKNOWN return if necessary. This avoids silently losing the result or starting a duplicate paid session.

## 18:30 CEST / 16:30 UTC checkpoint

Promotional balance **139.527881 USD**, paid-use preflight clear. Compared with the rounded 16:00 UTC balance (~153.07 USD), the shared account used ~13.54 USD in 30 minutes (~27.1 USD/h), close to the ~27.9 USD/h required over the remaining five hours. Seven W4 lanes launched sequentially: A907-PARITY, A363-ALGO, A295-EXACT, A303-BOUND, A194-ARMS, PATELLA-RANGE and B3-CELL. Forty W4/W5 lanes remain unstarted. One W3 cell result arrived after a same-session bounded follow-up. Next decision at 17:00 UTC based on exact balance delta and settled per-session costs.

## 18:36 CEST / 16:36 UTC closed W3 costs

All five W3 reserve sessions have returned files. Their settled API costs total **7.372789 USD** (mean 1.474558), including 2.723837 USD for the longer CELL-COUPLE follow-ups. The synthetic findings and independent checks are in `results/CX-CLOUDWAVE3/WAVE_REVIEW.md`. Seven new W4 sessions are currently running. No extra shared-balance read was made between checkpoints.

## 18:47 CEST / 16:47 UTC replenishment

Fifteen W6 follow-up lanes are preregistered from the returned W3/W4 findings and the spatial-cell roadmap (H1–H5). This restores the unstarted buffer from 40 to **55** before the 19:00 checkpoint. Priority order places A363 column scaling, A295 dual certificate, A907 switches, A194 source audit, OED contrast and one cell-reduction test early. The A194 source-audit bundle includes pinned public Gait2354 XML from opensim-models commit `d9b05d470b1a481c222372c85b75772faf8f7792`, SHA-256 `b96156dcd721777c4d1ac5fbc2706fac7ee23be64c9c1669dcca39556a416702`, and a note that its credits say CC BY 3.0 but demonstration use, not research validation. The patella data audit has a checked PubMed abstract note. No private geometry was bundled.

## 18:55 CEST / 16:55 UTC review state

Six of the seven 18:30 W4 packages have returned through events; B3-CELL is still computing. W4-A363 and W4-A303 initially omitted runnable code from their final messages; the missing code was fetched by follow-up in each original session. An A194 source-finding error was corrected in its original session after a direct check of licensed OpenSim Gait2354. Another read-only `git fetch origin --prune` still showed only `origin/main`, so no cloud result branches exist to import. The synthetic audits did not verify private A907/A363/A295/A303/A194 numbers. A :00 balance decision is next.

## 18:57 CEST / 16:57 UTC pacing rule for next blocks

`next_wave.sh auto` now makes a single fresh :00/:30 promotional/paid-use preflight, logs the exact balance, and chooses 2–9 starts from a linear drawdown target anchored at 139.527881 USD at 16:30 UTC to zero at 21:30 UTC. It starts with five and changes by one lane per 2.50 USD above/below that target, with a conservative cap from remaining credit (2 USD per proposed session). The decision is logged in `tasks/cloud/auto_decisions.jsonl`. This adapts the BodyTwin share to Field's joint balance without a second rate-limited usage read. Frozen lane criteria and launch sequence remain unchanged.

## 19:00 CEST / 17:00 UTC checkpoint

Promotional balance **116.610279 USD**; paid-use preflight clear. Shared decline since 16:30 UTC is **22.917602 USD in 30 minutes** (~45.84 USD/h), above the roughly 27.9 USD/h straight-line drawdown. The auto rule computed a target of 125.447719 USD, gap −8.837440 USD, and launched **2** W6 follow-ups: A363-COLSCALE and A295-DUAL. Fifty-three preregistered lanes remain unstarted. W4-B3-CELL has run about 30 minutes; a bounded-stop follow-up asked its original session to return completed cases within five further minutes and mark unfinished cases UNKNOWN. No second balance read was made.

## 19:12 CEST / 17:12 UTC returned follow-up costs

Both 19:00 W6 follow-ups returned through events and were rerun locally. API cost: A363-COLSCALE 1.089137 USD, A295-DUAL 1.199681 USD (sum 2.288818 USD). Both pass their separate synthetic frozen gates; limitations are in the two-line review. W4-B3-CELL remains running with no settled per-session cost. This makes the shared balance, not the two cheap BodyTwin receipts, decisive at 19:30.

## 19:17 CEST / 17:17 UTC follow-up input hardening

Seven unstarted W5/W6 follow-up bundles now include the exact returned prior `RESULTS.md`, `results.json`, relevant original code and SHA-256 manifests under `prior/` (A907 switches, A194 public source, OED contrast, vascular rank, patella data, cell reduction and A303 coverage). This lets the new cloud session compare directly rather than reconstructing a predecessor from brief prose. The source files are only public/own synthetic work; no private matrices or credentials were included. Each bundle remains below 20 MB and the frozen PREREG files are unchanged.

## 19:23 CEST / 17:23 UTC W4 closeout

All seven W4 sessions from 18:30 have returned RESULTS and JSON by paginated events. Their settled per-session API cost is **12.699457 USD**, including B3-CELL 1.874270 USD. B3-CELL was bounded after a long sweep: λ=25/50/100 nm completed, λ=200/400 nm UNKNOWN; the returned raw run log and public S3 metadata were checked locally. Two W6 follow-ups launched at 19:00 also returned and cost 2.288818 USD. The seven W4 and two W6 two-line reviews are appended in WAVE_REVIEW.md. Next shared balance read is 17:30 UTC.

## 19:26 CEST / 17:26 UTC public cell geometry follow-up

I independently fetched the public CellMap s2 mitochondrial crop used by W4-B3-CELL; decoded C-order SHA-256 exactly matched the cloud report (`2dd499ca…27ab93063edc7a876fd`). A 12.8 KB compressed copy plus voxel metadata and provenance now sits in the unstarted W6-CELL-REDUCE bundle, together with the earlier crop-closure summary/log. It is explicitly an illustrative cropped geometry, not a measured transport outcome or full cell. Bundle size remains about 124 KB. The W6 frozen synthetic H1 gate was not changed.

## 19:30 CEST / 17:30 UTC checkpoint and pacing correction

Promotional balance **86.881339 USD**, paid-use preflight clear. Shared decline since 17:00 UTC was **29.728940 USD in 30 minutes** (~59.46 USD/h), much faster than the now required ~24.8 USD/h over the last 3.5 hours. The target balance for this slot was 111.480157 USD, a −24.598818 USD gap. The previous auto rule had a floor of two starts, and it launched A907-SWITCH and A194-PUBLIC. I removed that floor immediately after this batch; future `auto` decisions can choose **zero** starts when the shared account is this far ahead of the straight-line drawdown. This correction preserves promotional-only and paid-use preflight on every slot. Fifty-one preregistered lanes remain unstarted.

## 19:40 CEST / 17:40 UTC 19:30 batch closeout

Both W6 follow-ups launched at 19:30 returned via paginated events. A907-SWITCH passed its frozen synthetic switch/regular gates (100/100 each); local rerun reproduced them, while the prior W4 detector tied on the gate and private A907 remains UNKNOWN. A194-PUBLIC found a licensed but demo-only Gait2354 model and no verified matched numeric measured curve or research-model license, so the frozen data-readiness verdict is UNKNOWN. I independently checked the bundled XML and found public candidate repositories, but no validated model-versus-measurement comparison. The A194 session initially omitted its cited audit code/raw file; both were fetched in the same session and a local rerun reproduced the parsed raw data and UNKNOWN gate. Settled API costs are 2.106238 and 1.726213 USD (3.832451 USD together). All 16 sessions newly launched by CX-CLOUDWAVE3 have now returned; 51 preregistered lanes remain unstarted. No off-slot shared-balance read was made.

## 20:00 CEST / 18:00 UTC value-qualified block

The promotional preflight returned **53.087151 USD remaining** and no paid usage enabled/detected; the shared account declined 33.794188 USD since 17:30 UTC. Under Anton's updated value rule, four task bundles were chosen for real public data or named method blockers rather than by the old synthetic reserve order: CLOUD-W6-STRENGTH-ID (public SC Biodex curves, held-out speed and explicit gauge), CLOUD-W6-KNEE-MEASURE (calibrated force endpoint and null-model falsifier), CLOUD-W6-FRAME-GAUGE (absolute CT/fluoro frame acquisition and rank), and CLOUD-W6-PATELLA-DATA (six measured curves retrieval). All four got session receipts at 18:00–18:02 UTC; they were running at the first no-usage collection. The older auto pace rule is superseded by the value criterion. Next shared balance read at :30, with paid-use preflight before any launch.

## 20:07 CEST / 18:07 UTC launcher guard

`next_wave.sh` now uses only the explicit `priority_order.txt` allowlist. The former fallback over every W3–W6 bundle could launch self-generated synthetic exercises, contrary to Anton's 19:50–20:00 value update. The two unstarted allowlisted reserves are the public NHANES cohort-transport check and the public CellMap geometry reduction. Additional launches require a fresh value-qualified entry and the scheduled promotional/paid-use preflight.

## 20:14 CEST / 18:14 UTC partial patella recovery

Three of the four 20:00 sessions returned complete FILE blocks and have local independent checks in WAVE_REVIEW.md. CLOUD-W6-PATELLA-DATA reached `requires_action` on a Bash step to clean scratch links and make a local commit; no push was present. The account rejected CLI attach to the existing cloud session. Its authenticated `Write` events contained a complete draft RESULTS.md and two draft code files, copied with an explicit RECOVERY.md; final `results.json`/raw curves are missing, so the model numbers are not independently reproduced. The preregistered empirical gate is UNKNOWN (0/6 measured curves), regardless of the blocked finalization. No duplicate session was started.

## 20:26 CEST / 18:26 UTC patella numerical recovery

The blocked PATELLA-DATA session still has no final FILE-block return. I recovered its producer-written report and code edits from authenticated events, then ran the reconstructed final code locally with OpenSim 4.6 and the pinned public opensim-models commit. The reported model peaks reproduced (53.562/49.792 mm FSA at 0°); `results.json` and 1° curves are labelled local reconstructions in RECOVERY.md. This does not alter the 0/6 empirical UNKNOWN gate. No duplicate cloud launch or paid usage was initiated for recovery.

## 20:30 CEST / 18:30 UTC public-data block

Promotional preflight returned **24.076935 USD remaining**, no paid use. Shared decline since 18:00 UTC was **29.010216 USD** in 30 minutes. Two explicit value-qualified reserves launched sequentially: CLOUD-W6-CELL-REDUCE (public CellMap crop plus output-aware operator comparison) and CLOUD-W4-A359-FAIR (public CDC NHANES XPT cohort transport against the adjusted baseline). Both have preregistered gates and separate source manifests. No old synthetic-only fallback was started. With shared decline near 29 USD per half hour, the promotional balance may be exhausted before 21:00 CEST; check only at the next :00 and stop on any paid-use indication.

## 20:40 CEST / 18:40 UTC public-data return

CLOUD-W4-A359-FAIR returned full FILE blocks. Its preregistered sex×age transport gate failed on public NHANES: the primary DXA-only feature gave −11.46% held-out RMSE gain versus adjusted anthropometry; a local independent XPT rerun reproduced n=1494/1485 and 84.431/94.110 N. The settled session API cost is 1.552663 USD. CLOUD-W6-CELL-REDUCE remains running. No off-slot account balance read or new launch occurred.

## 20:50 CEST / 18:50 UTC 20:30 batch closeout

CLOUD-W6-CELL-REDUCE returned report and JSON; its code was omitted from FILE blocks and reconstructed from authenticated session Write/Bash events, with provenance in its RECOVERY.md. The frozen synthetic numerical gate passed 24/24; a local nine-rate rerun for held-out seed 2000 reproduced 1.047% maximum error and its geometry counts. Its settled session cost is 2.377845 USD. CLOUD-W4-A359-FAIR and CLOUD-W6-CELL-REDUCE are both idle; there are no other active BodyTwin jobs in this lane. Shared balance remains last measured at 18:30 UTC (24.076935 USD); next account read is 19:00 UTC. No off-slot balance read or new launch occurred.

## 21:00 CEST / 19:00 UTC checkpoint and next qualified reserve

Scheduled promotional-balance read returned **13.796661 USD** (236.203339/250 used), locked_reason null; extra paid usage disabled, paid spend disabled, paid amount 0. Shared decline since 18:30 UTC was 10.280274 USD in 30 minutes. The earlier 29 USD/30-min rate did not persist, so assuming early exhaustion would have been wrong. No launch was possible after the balance slot without an additional off-slot preflight. A new value-qualified `CLOUD-W6-MUSCLE-EMG3` bundle is being prepared from 12 hash-checked public Grand Challenge Biodex/EMG files for DM/JW/SC and the previous failed SC benchmark. Its frozen test is held-out 90°/s prediction with a ≥10% gain on both directions in at least two persons, plus a separate muscle-parameter gauge check. This is a distinct data-driven follow-up to the three-person measurement packages, not a repeat of the SC-only synthetic fit or Field BODYGRAPH. Launch only at a subsequent :00/:30 slot after a fresh promotional/paid-use preflight and only if the shared credit remains.

## 21:30 CEST / 19:30 UTC qualified EMG launch

Promotional preflight returned **8.904557 USD** remaining, with extra paid use disabled and paid spend 0. Shared decline since 19:00 UTC was 4.892104 USD in 30 minutes, again slower. `next_wave.sh 1` launched the value-qualified public three-person EMG/strength transfer task `CLOUD-W6-MUSCLE-EMG3`, session `session_012bdbP9Sxa6hAk8YkyfnHur`. The bundle was 10.1 MB, hashes verified, with no implant sensor calibration assumed. No other lane started. Next shared balance read is 20:00 UTC; collect session events without balance in the meantime.

## 21:46 CEST / 19:46 UTC three-person EMG return

CLOUD-W6-MUSCLE-EMG3 returned all promised files and settled session API cost **3.697110 USD**. Frozen held-out 90°/s torque gate passed for DM and SC, failed for JW; individual muscle parameters remained rank 3/6. A local full rerun from the bundled public CSVs and returned code reproduced all six relative gains exactly to four decimals and the result JSON verdict. Event-order audit verified the method/prediction freeze preceded the first 90°/s torque parser use; sensitivity to a 1.25× between-trial EMG gain can reverse the gate. This is not an implant-force or muscle-parameter validation. Shared account balance remains last measured at 19:30 UTC; next allowed read is 20:00 UTC.
