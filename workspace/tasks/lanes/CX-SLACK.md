# CX-SLACK — identify the muscle's operating range (optimal fibre length/tendon slack) from measured strength at several knee angles

Lead: A361 + A362 (`results/CX-GC-STRENGTH`, `results/CX-GC-STRENGTH2`).
- JW's strength factor is 0.71 at 30° but 9.37 at 90°.
- At 90° the model's quadriceps fibres sit at 1.79 L0 (active force–length factor 0.034), giving 10 Nm capacity against a measured 96 Nm.
- So there is no single strength factor. The model's force–length operating range is wrong. That fits L1 (A312): TLEM's slack lengths do not fit our straight-line paths (soleus ~2× L0, 27–56/166 elements out of range).
- the reference model calibrates tendon slack in its own calibration studies (HipMuscleCal, according to CX-PARAM). Our variant: identify it from the person's measured strength–angle curve.

## Data
Grand Challenge:
- JW isometric at 30° (gc1) and 90° (gc4); isokinetic 60/90°/s files exist in the same archives, which give strength over the whole angle range;
- DM gc6 at 90°;
- SC gc3 at 90°.
Extract only what is needed (`/media/anton/sdc1-tmp/bodytwin/CX-GC-STRENGTH/` already has the synced archives). Isokinetic data needs a force–velocity correction; use L1's S1b curve and state it.

## Tasks (PREREG.md + sha256 before the first run)
1. Measured moment–angle curve per person (isometric points + isokinetic curves, corrected for velocity, gravity/limb weight).
2. Parameters to identify per muscle group (quadriceps, hamstrings, gastrocnemius): Fmax factor, L0 factor, and tendon slack factor. That is 3 per group; take care with identifiability. Report the Fisher information/singular values and which combinations are identifiable from 1 vs 2 vs many angles.
3. Fit on the calibration curves; hold out one angle/velocity as a test. Also compare against the reference model-style slack calibration (slack set so that fibres are at L0 in a reference posture) as a baseline.
4. Rerun L1's V0/S1b gait with the identified parameters for JW (and DM/SC if the curve allows). Does knee force approach implant force/N1g?
5. Counter-tests: another person's parameters, and a permuted angle order.

Criteria:
- the held-out angle is predicted within 20 % of the measured moment;
- gait V0 median ≥ 10 % lower than the base for the persons with curves;
- the counter-test gives a smaller gain.

Resources: lane runner has full permissions. Internal TLEM data stays local or on OVH. Locally: nice, 2 threads. Write in `results/CX-SLACK/`. `RESULTS.md` starting with `# CX-SLACK`, plus results.json and code with pytest.
