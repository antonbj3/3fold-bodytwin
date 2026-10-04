# CX-CT2MARKER — place each Grand Challenge person's individual bone geometry correctly in the gait frame (CT → marker transform), and fix the patella mechanics

Leads:
- A366 (`results/CX-JWGEOM`): individual femur/tibia/patella registered to 2–4 mm for JW/DM/SC/PS. But there is no verified transform from the CT frame to the gait markers' segment frames, so the geometry is individual but wrongly placed against the motion. The physics bound got worse.
- A367 (`results/CX-PATHS`): L1's patellar tendon moment arm peaks at 0° (51 mm), while Krevolin 2004 has its peak around 45°. The patella mechanics are suspect.

## Tasks (PREREG.md + sha256 before the first run)
1. For each person, look in the Grand Challenge archives (`external_media*Competition-latest.zip`; README/Geometry Data; the static trial; marker positions in CT if any; implant coordinate systems; the fluoroscopy/kinematics files in later competitions) for the data that link CT and markers. Inventory it per person.
2. Build the transform per segment:
   - registration of marker-based anatomical landmarks (epicondyles, malleoli, ASIS via the static trial) against the same landmarks on the individual bone;
   - estimated skin offset (A280/N51: ~20 mm AP);
   - a residual per landmark.
   Criterion: residual ≤ 10 mm median. Report which segments are identifiable.
3. With the transform: rebuild L1's A/b with the individual geometry (reuse CX-JWGEOM's registrations and point transfer) and rerun C1 (CX-BEATN1G2). Report JW's corrected frames specifically.
4. Patella: compare the model's patellar tendon moment arm against Krevolin 2004 and other DOI sources over 0–100°. If the shape is wrong (peak at 0°), test a patella tracking that follows the femoral groove from the individual geometry. Report the effect on the quadriceps strength curve (CX-SLACK's held-out 90°/s) and on knee force.
5. Counter-tests: another person's transform, and permuted landmarks.

Resources: lane runner has full permissions. TLEM/restricted model data stays local or on OVH. Locally: nice, 2 threads. Write in `results/CX-CT2MARKER/`. `RESULTS.md` starting with `# CX-CT2MARKER`, plus results.json and code with pytest.
