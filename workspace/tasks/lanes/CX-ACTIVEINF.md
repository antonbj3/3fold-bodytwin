# CX-ACTIVEINF — does co-contraction follow uncertainty (active inference) rather than mechanics? A test on the Grand Challenge

Read `results/CX-ACTIVEINF/FROZEN_HYPOTHESIS.md` (AI-1..AI-4 are frozen; `FROZEN.sha256` is created by the coordinator; do not change it).

Data:
- `results/CX-INVERSEOC/locations.json` (lo/hi/meas per frame) and `epsilon_results.json`, plus the cache in external_media*.npz (ε₂ per frame);
- `results/L1/prep/*.npz` (emg 15 channels: semimem, bifem, vasmed, vaslat, rf, medgas, latgas, tfl, tibant, peronl, soleus, addmagnus, gmax, gmed, sartorius; b = moments; grf);
- the original GC files for the recording order (see `results/L1/code/n12_io.py`, `l1_prep.py` for paths: file timestamps, trial numbers, the session description/README in the GC docs).

## Tasks (PREREG.md + sha256 BEFORE the test)
1. **Establish the CHRONOLOGICAL order per session.** Use the GC documentation, the timestamps in the data files, or the file metadata. Report how reliable it is. If only the trial number exists: use it, and say so.
2. **AI-1:** a within-series slope of median (meas−lo), and of ε₂, against the order, per series. Use a mixed model or the pooled Spearman.
   - Compare altered with normal gait.
   - Permutation test of the order within the series.
   - Controls: walking speed (the pelvis marker speed proxy), GRF peak, knee moment peak. A drift in speed or load must not explain the effect.
3. **AI-2:** altered vs normal gait within a person, matched on GRF level.
4. **AI-3:** the phase profile of meas−lo and of the hamstring EMG, after controlling for knee moment.
5. **AI-4:** the hamstring and quadriceps EMG slopes across repetitions.
Also report JW4's non-gait activities if they are repeated.

Every outcome is a node that expands. If AI-1 does not hold, report what the excess follows instead (mechanics, fatigue, speed). Deliver RESULTS.md starting with `# CX-ACTIVEINF`, results.json, the script, and pytest. Light computation: run locally with 2 threads. Internal data stays local. lane_runner has full permissions in the workspace; `~/projects/bodytwin` is read-only.

## ADDENDUM (the coordinator, before the lane's PREREG, 23:10): measure the RIGHT thing in the RIGHT place. This takes precedence over the tasks above
Anton asked whether we measure the right thing. We do not fully. meas−lo mixes co-contraction with MODEL ERROR: lo carries a known geometric offset (A1425, C0 ~−300 N; 22 % of frames below lo). If the kinematics drift across repetitions, lo moves, and that can look like "learning".
1. **PRIMARY measures (model-free), in AI-1..AI-4:**
   - (a) the measured hamstring EMG (semimem+bifem) at matched knee moment. Per person, fit a regression EMG ~ |M_knee| + GRF on ALL trials, and take the residual. Report the quadriceps and gastroc as controls.
   - (b) the measured implant force's residual after a per-person regression on ‖GRF‖ and |M_knee| (+ |M_ankle|) in the same phase bin, with NO LP model.
   meas−lo and ε₂ are SECONDARY. Report how much of any trend in meas−lo comes from a change in lo versus a change in meas.
2. **THE RIGHT PLACE:** the primary window is heel strike to 20 % of stance (loading response), defined from the GRF threshold. Report mid/late stance separately. Do not use medians over the whole trial window.
3. **Selection:** ε₂ exists only on LP-inside frames. Report the share of inside frames per trial, and whether it changes with repetition.
4. **Labels:** do NOT use L1's normal/altered for the test. Classify from the GC documentation: habitual/normal gait (ngait, ngait_og, ngait_og_ss, smooth? state what the docs say) vs modified gait modes (mtgait, tsgait/trunksway, medthrust, mtpgait, wpgait, bouncy, crouch). Report the classification and its source before the test.
