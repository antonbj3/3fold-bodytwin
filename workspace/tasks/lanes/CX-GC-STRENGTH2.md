# CX-GC-STRENGTH2 — T3 for more persons: model capacity at 90° (DM, SC, JW gc4) and a combination with individual geometry

Read first:
- `results/CX-GC-STRENGTH/RESULTS.md` (A361): JW's measured isometric knee strength at 30° shows TLEM's muscles are 2.5–4× too strong. Calibration lowered V0 1.060 → 0.704 BW. DM/SC were measured at 90°, outside the gait angles, so they were not used.
- `results/CX-JWGEOM` (running or done: individual knee geometry).

## Tasks (new PREREG.md + sha256 before the first run)
1. Compute the model's knee extension/flexion capacity at 90° (and at 30°) from a static pose, instead of gait frames. Use L1's/TLEM's geometry (results/L1/code, tlem_arch.json) with the knee set to 90°: moment arms from the muscle paths at that angle, F0 × force–length (L1's S1a curve) × moment arm.
2. Calibrate the quad and flexor group factors against the measured strength for DM (gc6 isomet90), SC (gc3 isomet90), and JW (gc4 isomet90 plus gc1 30° as a consistency check: does the same factor come out at both angles?).
3. Rerun V0/S1b gait with calibrated Fmax for DM, SC and JW. LOPO criterion: V0 median ≥ 10 % lower in ≥ 2/3 persons. Counter-test: another person's factor.
4. If CX-JWGEOM has delivered individual geometry, combine the two (geometry + strength) for JW and report both separately and together.

