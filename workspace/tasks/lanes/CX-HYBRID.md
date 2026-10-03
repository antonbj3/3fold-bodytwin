# CX-HYBRID — does the MODEL add information BEYOND N1g? (the model's shape + N1g's level, LOPO)

Background:
- N1g (k·‖GRF‖) wins against every pure model prediction (A1772 QUADPRED: stress-2 0.718, metabolic 0.538, ε-midpoint 0.839 vs 0.396 BW).
- With the corrected quadriceps arm, the model's lower bound lo_q is more physical (under lo 22.55 % → 9.71 %).
- meas−lo is phase-dependent (A1418).
- The real question is not "model OR N1g" but whether the model's time-resolved SHAPE adds anything to N1g.
Read: `results/CX-QUADPRED/` (score.py, per-frame caches), `results/CX-QUADARM/geometry.py`, `results/CX-MINCONTACT/analyze.py` (LOPO c), and `results/CX-INVERSEOC/locations.json`.

## Predictors (PREREG.md + sha256 BEFORE scoring; everything LOPO over persons; the held-out person's implant force NEVER used)
- (a) lo_q + c(phase): c per 5 phase bins, learned on 3 persons.
- (b) N1g · g(lo_q/‖GRF‖): a monotone 1-D correction of N1g by the model's ratio (isotonic or a 2-parameter linear fit, learned on 3 persons).
- (c) the blend w·N1g + (1−w)·(lo_q + c(phase)), with w learned on 3 persons.
- (d) the same as (c), but with metabolic-proxy_q instead of lo_q.
- A control: the same learning with lo_q time-shifted by half a trial (it must lose), and with L1's original lo (how much does the arm correction add?).
Criterion (frozen): the person-median RMSE is ≥ 5 % below N1g (0.396 → ≤ 0.376 BW) AND better than N1g for ≥ 3/4 held-out persons, AND the time-shift control loses that gain.
Report the peak and early stance too.

Deliver RESULTS.md starting with `# CX-HYBRID`, results.json, the script, and pytest. Runtime ≤ 60 min; reuse the caches; 2 threads under bigmem.lock. At most 200 MB of intermediate files (the disk is almost full; write to /media/anton/sdc1-tmp). Internal data stays local. lane runner has full permissions in the workspace; `~/projects/bodytwin` is read-only.
