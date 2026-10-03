# CX-NONGAITOPS — build L1 operators for DM/SC/PS non-gait trials (squat, chair rise, calf raise, one-leg stand, crouch, bouncy)

Background: L1 operators exist only for gait (113 trials) and JW4 non-gait (11 trials, `results/CX-INVERSEOC/non_gait.py` → /media/anton/sdc1-tmp/bodytwin/CX-INVERSEOC/non_gait/). BT-LG (A1419) shows that the law and N1g diverge the MOST in non-gait. Field U472 is testing whether heavy tasks are more informative (JW4). DM/SC/PS need the same chain.

## Tasks
1. Copy the non_gait.py pattern: static calibration per session (DM6, SC3, PS5; statics in n12_io SESSIONS), marker QC, foot reconstruction, and plate choice by nearest implant-side foot (implant side DM R, SC L, PS L).
   Build A, b, F0, cj, C0, Mkx, mx, d, Lmt, emg (if it exists), meas_tot, grf, angle, and law for every non-gait trial in `results/BT-LG-*` (list them; exclude jw_lungef1).
2. **Cross-checks before any use:**
   - the law via this ID must reproduce CX-LGSYNTH's law RMSE per trial (tolerance 0.02 BW; report deviations and why: the ID method differs from the BT-LG proxy);
   - the LP lo/hi via HiGHS; report the feasible share per trial.
3. Save to /media/anton/sdc1-tmp/bodytwin/CX-NONGAITOPS/<SESS>__<trial>.npz with a json per trial (QC). Write RESULTS.md starting with `# CX-NONGAITOPS`, with a table per trial of QC, feasible share, lo/meas/hi medians, and N1g/law RMSE.
No hypothesis test here: this is infrastructure for the heavy-task node. Run under bigmem.lock with 2 threads. Internal data stays local. lane runner has full permissions in the workspace; `~/projects/bodytwin` is read-only.
