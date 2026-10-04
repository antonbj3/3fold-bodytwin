# CX-GEOMPRIOR — bring the The swarm analyses' geometry priors into the geometry manager (integration, after verification)

Source: BT-DA-027…038 (read RESULTS.md). The positive findings to verify and integrate:
- DA-034: height+weight+sex → femur length (VSD, person-held-out);
- DA-035: one femur's head radius → the contralateral one;
- DA-030: tibial length + plateau width → femoral head radius;
- DA-033: the Imperial Gaussian prior transfers to VSD for 3 scalar parameters;
- DA-029: the left-right difference as a location prior (the interval does not transfer).
Also note the "no" results (DA-028 anteversion transfer, DA-031, DA-032, DA-036, DA-037) as priors that are NOT to be used.
## Tasks
1. Recompute each positive finding with your OWN code from the source tables (independent verification; The swarm < Sol). Verdict: HOLDS / DEVIATES / A MATTER OF DEFINITION.
2. Integrate the ones that hold as conditional priors in `bodytwin_core/geometry` (e.g. `prior_from_anthropometry(height, weight, sex)`, `contralateral(...)`), with prediction intervals calibrated per source, and a pytest that reproduces the key number.
3. Report in the README which priors exist, their held-out error, and their validity range.
Write in `results/CX-GEOMPRIOR/`. `RESULTS.md` starting with `# CX-GEOMPRIOR`.
