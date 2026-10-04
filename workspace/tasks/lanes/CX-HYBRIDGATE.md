# CX-HYBRIDGATE — use the model's shape ONLY where the geometry is sound, detected WITHOUT the implant

Background: A1775 CX-HYBRID. The blend (c) of N1g + (lo_q + phase c) beats N1g for DM/SC/PS (PS −27 %). JW gets worse, and JW has a known broken patella geometry: the patellar tendon arm is 3–6 mm at 40–60° and ρ is clipped (A1450/A1452/A1703/A1774). Excluding JW afterwards is not allowed. A validity gate computed from the MODEL/MARKERS only (NO implant force) is allowed if it is frozen in advance.

## Freeze in PREREG.md (+sha256) BEFORE scoring
Gate per frame (or per trial; choose one and freeze it), implant-free: the model is "sound" if
- the derived patellar tendon arm A[3]/ρ for the quadriceps elements is ≥ 15 mm (median over the elements),
- AND ρ is not at 0.4/1.3,
- AND the LP is feasible.
Predictor: (c) from CX-HYBRID (the same LOPO-learned w and c; NOT relearned on the gate) where the gate is sound; N1g otherwise.
Criterion: person-median RMSE ≤ 0.376 BW AND better than N1g for ≥ 3/4 held-out persons.
Also report:
- the share of frames that are sound, per person;
- the result with the gate applied during the learning too (LOPO, as a secondary variant);
- a control: a random gate with the same share per person must NOT give the same gain.

Reuse `results/CX-HYBRID/` (score.py, caches), `results/CX-QUADARM/mechanism.py` (ρ, arm per frame) and `results/CX-PATELLAJW/`. Runtime ≤ 45 min, 2 threads under bigmem.lock, at most 200 MB of intermediate files on external_media Deliver RESULTS.md starting with `# CX-HYBRIDGATE` and results.json. lane_runner has full permissions in the workspace; `~/projects/bodytwin` is read-only. Internal data stays local.
