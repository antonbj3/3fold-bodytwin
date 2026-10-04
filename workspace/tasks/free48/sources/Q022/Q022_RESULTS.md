BT-HX-Q022

Builds on: `BRIEF.md`, `inputs/NIGHT_PREAMBLE.md`, `inputs/QUESTION.md`; the referenced BodyTwin repository/anchors were not mounted, so no node was invented. Public primary anchors are Fabiato 1981 (DOI 10.1085/jgp.78.5.457, Fig. 10: pCa 5.30–5.40 and full activation near 4.90), Schouten 1992 (DOI 10.1113/jphysiol.1992.sp019180: 74 +/- 20% peak-force change, 0→70 cmH2O), Schulz 1991, and Goto 1991.

`model.py` is a one-thread normalized lumped chain: electrical activation gate → SR/cytosolic Ca mass balance → Hill/length/velocity active force → strain oscillator; perfusion pressure gives saturating flow, oxygen supply, calcium reload, and a Gregg-like gain. The electrical-off placebo, causal ordering, unit ledger, and all frozen checks pass. Results are model-generated, not measured data; details and provenance are in `results.json`.

At the synthetic pulse, p=0 gives peak active force 0.769273 and p=70 gives 1.334063 (dynamic ratio 1.734187); the static gain ratio is 1.74. Electrical→activation lag is 0.006 s, activation→force lag 0.033 s, and peak cytosolic Ca is 4.536488e-6 M. Flow is 0, 0.5, 1.0 at 0, 35, 70 cmH2O. The ±50% sensitivity covers `k_on`, `ca_ec50`, and `gregg_gain`; the latter spans -13.854% to +13.750% peak force and fails the static-anchor check when perturbed, as expected.

What fell: no paired human/internal traces, absolute force/flow calibration, or distributed-heart validation. Next: synchronized electrical/ECG, strain/force, and perfusion pressure/flow traces; estimate parameters and test the frozen lag, saturation, and placebo criteria.
