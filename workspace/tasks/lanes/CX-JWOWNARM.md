# CX-JWOWNARM — JW's quadriceps arm from JW's OWN OpenSim model (in place of the DM transplant/Rajagopal)

Background:
- A1703 CX-JWPATGEOM: the archive has JW's individual OpenSim model with three patellar ligaments (lat/center/med; rest lengths 34/26/30 mm) and CT bones/implant. A frame transform to L1's marker frame is missing.
- The moment arm as a function of knee angle is, however, intrinsic to the model and needs NO frame transform. Only the knee angle needs to be matched.
- A1449/A1702: Rajagopal's arm (another person's) already lowers JW gait under lo from 48 % to 21 %.
- A1775/A1776: JW is what drags the hybrid down.

## Tasks (PREREG.md + sha256 FIRST)
1. Load JW's OpenSim model (OpenSim 4.6, as in CX-L1ARMS; paths in results/CX-JWPATGEOM/results.json). Compute the knee extension arm for the quadriceps (via the patella/patellar ligament mechanism in that model) and cj over 0–90°. Compare against L1's, Rajagopal's and Im et al.'s.
2. Swap ONLY the quadriceps arm for JW (JW1 + JW4 gait, JW4 non-gait), as in CX-QUADARM/geometry.py. Report the LP share under lo.
   Frozen criteria: JW gait under lo ≤ 20 % (better than or equal to Rajagopal's 20.8 %), AND JW4 non-gait a reduction ≥ 25 %.
3. Run CX-HYBRID's predictor (c) with JW's own arm (the same LOPO-learned weights as in CX-HYBRID, NOT relearned). Report the person-median RMSE against N1g, and whether JW then beats N1g.
Deliver RESULTS.md starting with `# CX-JWOWNARM` and results.json. Runtime ≤ 60 min, 2 threads under bigmem.lock, at most 200 MB of intermediate files on external_media Internal data stays local; no jw_lungef1. lane_runner has full permissions in the workspace; `~/projects/bodytwin` is read-only.
