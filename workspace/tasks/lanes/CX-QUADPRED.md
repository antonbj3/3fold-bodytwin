# CX-QUADPRED — 30-MINUTE HUNT: does the corrected quadriceps arm make a MODEL predictor beat N1g?

Background:
- A1449/A1450/A1702: swapping only the quadriceps arm (Rajagopal; `results/CX-QUADARM/geometry.py` `knee_geometry(model="l1_quadarm_rajagopal", trial=key)`) lowers the frames under lo from 22.55 % to 9.71 %, with 60 better / 3 worse per trial.
- A1410 CX-INVERSEOC: with L1's original geometry, no muscle objective beat N1g (stress-2 0.742, metabolic proxy 0.573, vs N1g 0.396 BW person median). CX-QUADARM tested only lo+c (still worse).

## The task (HARD time box: done and written within 30 minutes; take a subset if needed and say which)
Freeze in PREREG.md (+sha256) before scoring:
- the predictors: stress-2 min-E, metabolic proxy (INVERSEOC's definitions in `results/CX-INVERSEOC/analyze.py`/`score.py`), and the ε=10 % band midpoint (`bodytwin_core/contact_band_batch.py`);
- ALL with the quadriceps-arm geometry, against N1g (L1 LOPO) on INVERSEOC's mask.
Criterion: person-median RMSE < 0.396 BW, or better than N1g for ≥ 3/4 persons.
Reuse caches; run the QP batched (contact_band_batch / analyze.qp) with 2 threads under bigmem.lock. If the full mask does not fit in time: 3 trials per person, stratified.
Deliver `results/CX-QUADPRED/RESULTS.md` starting with `# CX-QUADPRED`, results.json, and the script. lane runner has full permissions in the workspace; `~/projects/bodytwin` is read-only. Internal data stays local. NOTE the disk: external_media has ~2 GB free; write at most 200 MB of intermediate files.
