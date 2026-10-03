# CX-ROLLBACK — measured femoral rollback (fluoroscopy) → the extensor moment arm vs knee angle → quadriceps strength curve → the knee chain (BT-2)

Lead:
- Cloud W4-PATELLA-RANGE (`results/CLOUD-W4-PATELLA-RANGE`, review in `results/CX-CLOUDWAVE3/WAVE_REVIEW.md`): in a planar patella model, whether the moment-arm peak lands at 0° or around 45° is governed by ROLLBACK (Krevolin 2004: 4–6 cm peak near 45°, a 3D finite helical axis).
- A1270: JW's moving patella path has no valid candidate at 0/7 angles.
- A367: L1's arm peaks at 0°.
- A364/A362: the quadriceps strength curve cannot be reproduced.
- UNIQUE DATA: JW's fluoroscopy gives the measured tibiofemoral pose with "a(+)/p(-) translation" per frame over large flexion — measured rollback.

## Data
- Grand Challenge 4th competition `/media/anton/sdc1-tmp/3.DataforFourthCompetition-latest.zip` → fluoroscopy trials:
  - lunge: jw_lungef1, already extracted in `/media/anton/sdc1-tmp/bodytwin/CX-LUNGESYNC/`;
  - other fluoro trials in the same archive: stairs, step-up, openfe (open-chain leg extension, UNLOADED — important for separating load from geometry), twist.
- JW's implant geometry: `results/CX-JWGEOM` (the STL files are in METRES).
- Strength curves: `results/CX-SLACK` (isometric 0/30/60/90° + isokinetic).
- FORBIDDEN: eTibia/knee_forces files (Field's F-8 facit).

## Tasks (PREREG.md + sha256 first)
1. Rollback curve: femoral AP translation relative to the tibia vs flexion, per fluoro trial. Compare loaded (lunge, stairs) with unloaded (openfe). Watch the flexion sign: it is ≈ −100° in the lunge although the label says flexion(+). Determine the convention from openfe, where the knee angle is known from the markers.
2. Extensor moment arm from measured kinematics: the patellar tendon line from JW's implant/bone geometry, plus the instantaneous helical axis from the fluoro poses (a finite helical axis between neighbouring frames) → moment arm vs angle. Is the peak ~45° (Krevolin)?
3. Replace L1's patella arm with the measured one. Rerun CX-SLACK's strength fit (reuse its code) → does the held-out 90°/s extension curve fall from 35 % to ≤ 20 %?
4. If the strength curve improves: rerun V0/S1b gait for JW (L1 chain) and C1 (CX-BEATN1G2) → against N1g.
5. Counter-tests: a shifted/permuted rollback curve; DM's/SC's geometry with JW's rollback.

Resources: lane runner has full permissions in the workspace. Internal TLEM data stays local or on OVH. Locally: nice, 2 threads. Write in `results/CX-ROLLBACK/`. `RESULTS.md` starting with `# CX-ROLLBACK`, plus results.json and code with pytest.
