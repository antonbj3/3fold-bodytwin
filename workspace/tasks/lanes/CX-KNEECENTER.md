# CX-KNEECENTER — functional knee centre/axis (SCoRE/SARA) instead of the marker midpoint: how much of the law's, the LP set's and C0's error is actually knee-centre error?

## Background
- A1431 CX-RQANGLE: JW bouncy RMSE ranges from 0.705 to 1.647 BW depending ONLY on which marker defines the knee centre. The lever-arm curves are secondary.
- A1425 U447: the frames below Fmin are repaired by lowering the contact offset C0 by a median of −312 N, or by 3.9 % in contact arm. That is a geometric error of the same kind.
- A1413 Field FB_HUNT_KNEE_CENTER_SENS: in gait2392, SCoRE/SARA (±2.5 mm) gives 3.8 pp RMS against 6.1 pp for a 7.5 mm error.
- Code: `results/CX-LGSYNTH/aggregate.py`, `results/CX-RQANGLE/analyze.py`, `results/L1/code/` (how L1 defines the knee centre and C0), `results/CX-ANKLETERM`/`CX-LAW2`.

## Tasks (PREREG.md + sha256 BEFORE scoring against the implant force)
1. For each person/session (JW1, JW4, DM2, DM6, SC3, PS5), estimate the functional knee centre and axis with SCoRE (centre) and SARA (axis). Use the trials with the largest knee range of motion (gait, crouch, squat, calf raise; NOT jw_lungef1). Express the result in the thigh and shank marker clusters. Report:
   - the offset from the marker midpoint and from L1's centre, in mm;
   - the stability across trials (bootstrap over trials).
   The implant force is NOT used here.
2. Recompute the knee moment with the functional centre/axis and score:
   - the parameter-free law (constant arm + HR×Grood) on the 50 BT-LG and 108 L1 trials, against N1g;
   - L1's LP bounds lo/hi with the new geometry, on the INVERSEOC mask. How many of the 2,647 frames below Fmin remain? Does the change in C0 match U447's −312 N?
3. Frozen criteria: (a) the law's person-median RMSE improves by ≥10 % for ≥3/4 persons; (b) the share of frames below Fmin falls by ≥50 %.
4. Counter-test: a random centre shift of the same size as the SCoRE offset, in random directions, must not give the same improvement.

Deliver RESULTS.md starting with `# CX-KNEECENTER`, results.json, `functional_knee_center(session)` in bodytwin_core-compatible form, and pytest. The LP part goes under bigmem.lock with 2 threads, or on OVH via `tasks/cloud_run.sh` (finish by 07:30). Internal data stays local. lane_runner has full permissions in the workspace; `~/projects/bodytwin` is read-only. Every outcome is a node that expands.

## Addendum 23:20 (A1435, CX-EPSBAND)
L1's C0 = −F_knee · y_tibia (n12_model.null0). The ID knee reaction projected on the tracked tibial axis. The reaction force itself barely depends on where the centre is (force balance), but y_tibia does: the tibial long axis runs from the knee centre to the ankle centre.
So also recompute y_tibia with the functional knee centre (and a functional ankle centre/malleoli midpoint), then C0 and the LP bounds. Report:
- the change in C0 in N per phase, against U447's diagnostic −240/−404/−304 N (early/mid/late);
- how many of the 2,647 frames below lo remain.
A tilt of a few degrees in y_tibia under a large F_knee shear component can give hundreds of N.
