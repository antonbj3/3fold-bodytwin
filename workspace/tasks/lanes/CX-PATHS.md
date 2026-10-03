# CX-PATHS — are straight-line muscle paths the root cause? Wrapped, individual paths → moment arms → strength curves → knee force

Pattern from today: every lever on the muscle model falls in the same way. None of the following gives consistent gains:
- recruitment (A312);
- strength factor (A361/A362);
- Fmax/Lopt/slack from strength curves (A364: the model cannot reproduce measured strength–angle curves even with 3 free parameters per group; parameters end up at bounds);
- the physics bound for JW (A358).

Common factor: L1's muscle paths are straight lines without wrapping, on TLEM geometry (for JW, a DM transplant). Hypothesis: the moment-arm curves as a function of knee angle are wrong. That distorts strength curves, fibre lengths (1.79 L0 at 90°) and contact force.

## Build on (read first)
- `results/L1` (A matrices = moment arms, prep/*.npz with muscle–tendon lengths), `results/CX-SLACK`, `results/CX-GC-STRENGTH2`, `results/CX-JWGEOM` (if done).
- `results/N12b` (the fixwrap variants), `results/N6` (wrapping/contact fraction, A269), `results/N7b` (h: wrap cylinders fitted to the individual's bone), `results/N7c/geomgr`, `results/BT-XF9` (SDF band).
- `~/projects/bodytwin/docs/MECHANISM_MOMENT_ARM_VALIDATION.md` and `MECHANISM_BIARTICULAR_MUSCLE.md`; `~/projects/bodytwin/scripts/msk` (grep moment_arm, wrap, patella).
- Literature ground truth for knee moment arms as a function of angle, e.g. the patellar tendon moment arm and hamstrings. Only cite sources with a DOI; otherwise UNKNOWN.
- The field lane's SDF paths/geodesics, if they exist (grep the field engine for geodesic/wrap). Write "Builds on".

## Tasks (PREREG.md + sha256 before the first run)
1. Moment arm vs knee angle (0–100°) for quadriceps (via patella), hamstrings and gastrocnemius under three path variants:
   (a) L1 straight lines;
   (b) wrapped (cylinders/via points per N12b/N7b);
   (c) individually fitted wrapping on the person's bone (geomgr).
   Compare against the literature curves: RMSE in mm and shape.
2. With (b)/(c): redo CX-SLACK's strength–angle fit. Does the rank rise, do the parameters leave the bounds, and does the held-out 90°/s curve fall below 20 % error?
3. Redo the V0/S1b gait and C1 (the certified feasible set, CX-BEATN1G2) with the best path variant. Does JW's F_min become correct? Does the model approach N1g?
4. Counter-tests: permuted wrap parameters and another person's bone.
5. Criterion for "the root cause is the paths": (b) or (c) lowers the moment-arm RMSE vs the literature by ≥ 30 % AND the held-out strength curve by ≥ 30 %. Symmetric skepticism.

Resources: lane runner has full permissions. Internal TLEM data stays local or on OVH (`tasks/cloud_run.sh`, BodyTwin ≤ 12 vCPU, finish by 22:45). Locally: nice, 2 threads. Write in `results/CX-PATHS/`. `RESULTS.md` starting with `# CX-PATHS`, plus results.json and code with pytest.
