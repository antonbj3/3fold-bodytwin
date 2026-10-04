# CX-LUNGEID — external knee load per fluoro frame for JW's lunge (inverse dynamics) + total contact estimate, for Field F-8 (U391→ force-controlled contact)

Consumer: the Field lane. The displacement-controlled contact was ill-posed: the pose error is ~1 mm against the true deformation of ~0.1 mm. Field now solves a small pose correction so that the contact force balances a load. The comparison against eTibia then tests the medial/lateral distribution and the pressure field.

## Data (read-only; facit sealed)
- `external_media`: jw_lungef1_trajectories.csv (markers), jw_lungef1_grf.csv, the fluoro segments.
- `results/CX-LUNGESYNC/lunge_frames.csv` + SYNC.md: 84 frames, time base.
- FORBIDDEN: `jw_lungef1_knee_forces.csv` and `results/CX-LUNGESYNC/facit_etibia.csv` (eTibia = Field's facit).
- JW's OpenSim model and static trial are in the Grand Challenge 4th-competition zip (`external_media`); grep `~/projects/bodytwin/scripts/msk` for knee_jw / grand_challenge / inverse dynamics pipelines, and `results/L1`, `results/N12b` (the knee chain) for reuse. The OpenSim venv is `source_repository/.venv-msk/bin/python` (read-only use).

## Deliverables in `results/CX-LUNGEID/` (PREREG.md + sha256 first)
1. `knee_load_tibia_frame.csv`, one row per fluoro frame (84):
   - the net intersegmental knee force (3) and moment (3) in the TIBIA frame (state the axes and whether the definition matches the fluoro file's joint coordinates: a/p, s/i, m/l);
   - N and N·m;
   - method: ID from markers + GRF (OpenSim ID with JW's model, or a segment-by-segment Newton–Euler; state which);
   - uncertainty: marker filter ±, segment inertia ±.
2. `total_contact_estimate.csv`: the total contact force per frame from our models, stated separately:
   (a) N1g (already in lunge_frames.csv — refer to it);
   (b) the net ID load + quadriceps force from the knee extension moment / the patellar tendon moment arm (a simple muscle estimate, lever arm per knee angle from L1/TLEM or literature; state the source);
   (c) the medial/lateral split from the external adduction moment (a moment-arm/condyle-distance model, e.g. ~40–50 mm between condyles from JW's implant geometry).
   These are MODEL values for Field to balance against, not facit.
3. `README_FOR_FIELD.md` (≤ 25 lines): frames, signs, what the total contact force contains (external load + muscles), and why balancing against the net ID force alone would give too low a total.
4. `RESULTS.md` starting with `# CX-LUNGEID`: plausibility checks — the GRF moment arm, and the knee flexion moment's sign against the knee angle in the lunge.

Resources: local, nice, 2 threads, ≤ 60 s per step. lane runner has full permissions in the workspace. When done, write `results/CX-LUNGEID/READY_FOR_FIELD.md`.
