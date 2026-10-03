# CX-ROLLBACK2 — unlock the extensor moment arm: use Field's geometry-calibrated fluoro frame (U391) with an explicit mm uncertainty

Blocker in A912 (`results/CX-ROLLBACK/RESULTS.md`): "no verified rigid transform between the fluoro 6-DOF and the patellar tendon line/knee axis". The same blocker stopped A366 and A368.
Unlock: Field's U391 (F-8) searched 576 frame conventions (axes × Euler order × flexion sign × origin) against GEOMETRY + GRF only, with no facit. The best convention gives max 3 mm, p95 2.1 mm penetration. That is a frame with mm-level uncertainty, not a verified transform, but it is sufficient for a moment arm of 40–50 mm (error ~2–6 %).

## Sources
- U391: `the public staging tree/3fold-motion-engine/_private/romi_collab/build/U391/` (RESULTS.md, the chosen convention, code). Read only; use the convention with a source reference. If U391 is not found there, grep romi_collab/build for U391.
- `results/CX-ROLLBACK` (A912: AP rollback curve over 692 fluoro frames, code), `results/CX-JWGEOM` (JW's implant/bone geometry; the STL files are in METRES), `results/CX-SLACK` (strength curves), `results/L1` (the knee chain).
- FORBIDDEN: eTibia/knee_forces (Field's sealed F-8 facit).

## Tasks (PREREG.md + sha256 first)
1. Apply U391's convention to all fluoro trials (lunge, stairs, step-up, openfe, twist) → the femoral and tibial component poses in a common frame.
2. Extensor moment arm vs flexion: the patellar tendon line (from the patella/tibial tuberosity in JW's geometry) against the instantaneous finite helical axis between frames. Propagate the frame uncertainty (sample convention variants whose penetration is ≤ the best + 1 mm) → a band on the arm per angle.
3. Where is the peak? Compare with Krevolin 2004 (peak ~45°, 4–6 cm) and with L1's arm (peak at 0°).
4. If the arm is identifiable within ±10 %, replace L1's patella arm and rerun CX-SLACK's strength fit for JW → does the held-out 90°/s curve fall from 35 % to ≤ 20 %? Then V0/S1b gait against N1g.
5. Counter-tests: the second-best convention; a permuted frame order; openfe (unloaded) vs lunge (loaded).

Resources: lane runner has full permissions in the workspace. Internal TLEM data stays local or on OVH. Locally: nice, 2 threads. Write in `results/CX-ROLLBACK2/`. `RESULTS.md` starting with `# CX-ROLLBACK2`, plus results.json and pytest.
