# CX-HYBRIDOOS — the frozen hybrid on data it has NEVER seen: non-gait activities (squat, chair rise, stairs, calf raise, crouch, bouncy)

Background: A1775 CX-HYBRID. Predictor (c) = w·N1g + (1−w)·(lo_q + c(phase)), with the quadriceps-arm geometry. It beats N1g for DM/SC/PS in gait (LOPO). JW has a known broken L1 arm (A1450/A1452/A1777). All learning (w, c, N1g's k) was done on GAIT. The non-gait operators have never been used for choices in the hybrid, so they are real out-of-sample data:
- JW4: external_media*.npz (11 trials);
- DM/SC/PS: external_media*.npz (30 trials; note the law cross-check deviation 29/30).

## Freeze in PREREG.md (+sha256) BEFORE any non-gait scoring
- Use EXACTLY CX-HYBRID's learned parameters (w, c per phase bin; LOPO per held-out person as in CX-HYBRID; no new learning).
- N1g = L1's frozen gait coefficient per person/group, transferred unchanged.
- Phase in non-gait: define it from the GRF load window (0–100 % of the loaded interval) and freeze that.
- Criterion: the person-median RMSE of (c) < N1g's on non-gait, AND better for ≥ 3/4 persons.
- Report per activity group as well, including JW separately, and with a result excluding JW labelled as a SECONDARY analysis. JW's exclusion is justified by A1777 (the L1 arm is broken), but it must not be the primary result.
- Control: lo_q time-shifted by half a load window must lose the gain.

Reuse results/CX-HYBRID (score.py), results/CX-QUADARM/geometry.py, and bodytwin_core/contact_band_batch.py. Runtime ≤ 60 min, 2 threads under bigmem.lock, at most 200 MB of intermediate files. Deliver RESULTS.md starting with `# CX-HYBRIDOOS` and results.json. Internal data stays local; no jw_lungef1. lane_runner has full permissions in the workspace; `~/projects/bodytwin` is read-only.
