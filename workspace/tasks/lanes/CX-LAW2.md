# CX-LAW2 — expand the residual nodes of the parameter-free knee law (A1391): dynamics + person, and new activities

Read: `results/CX-ANKLETERM/RESULTS.md` (A1391: the ankle law gives median −2.1 % against N1g's peaks without fitting; the time series is weaker, within 10 % in only 19/108; strong person dependence, with JW worst), `results/CX-BTHYPO/FRAMING_NOTE.md` (the branch rule), and `results/CX-BTHYPO/COORDINATOR_PREDICTION.md` (do not change it).
## Tasks (a NEW PREREG.md + sha256 before running; a new version of the law with its own frozen additions)
1. Dynamic terms: add the segment's inertia/acceleration term to the net moment (full ID instead of quasi-static, if the law used quasi-static) and a time shift for muscle force (electromechanical delay 30–80 ms from the literature). Does the time series move closer to N1g?
2. The person node: is JW's deviation explained by his DM-transplant geometry? Test with individual lever arms from CX-JWGEOM's bone registration (A366) where it is valid, and report the effect per person.
3. New activities without training: compute ID (CX-LUNGEID's pipeline) for Grand Challenge stairs, step-up, chair rise and squat (NOT jw_lungef1), and predict k per activity with the law. Compare against the measured implant force (lbf → N). This tests the frozen point 5 (k_stairs/k_gait 1.3–1.6).
4. Report per point: the numbers against N1g/B24, and the residual decomposition (which node expands).
Write in `results/CX-LAW2/`. `RESULTS.md` starting with `# CX-LAW2`, plus results.json and pytest.
