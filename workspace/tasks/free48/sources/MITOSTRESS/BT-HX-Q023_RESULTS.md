BT-HX-Q023

Built: a reduced HPA delay-feedback model in `model.py` (CRH-like drive → ACTH → free cortisol → slower total pool, with delayed cortisol negative feedback), plus `test_model.py`, `PREREG.md`, `PREREG.sha256`, and `results.json`. `UNKNOWN`: no internal hormone trace or individual CRH/ACTH/cortisol trajectories were supplied; every trajectory value is a model prediction, not a measurement.

Source anchors in `results.json`: ACTH-to-cortisol response lag 3.0 ± 0.92 min, free-cortisol half-life 2.2 ± 1.3 min, total-cortisol half-life 49 ± 2.4 min, and ACTH half-life 20 ± 1.3 min (Dorin 2012, DOI 10.1210/jc.2011-2227; Keenan 2004, DOI 10.1152/ajpendo.00167.2004).

Default prediction: ACTH onset 4.6 min, free-cortisol onset 10.9 min (onset lag 6.3 min); peak lag 5.4 min; free and total peaks occur at 22.3 and 33.5 min. Feedback suppression contrast is 0.08999 versus 0 without feedback; trajectories remain finite and free-cortisol peaks before total cortisol.

The frozen 1–5 min onset-lag criterion failed at the default (6.3 min). The ±50% runs show `delay_ac` onset lags 4.8/6.3/7.8 min, `delay_fa` suppression contrasts 0.03610/0.08999/0.13359, and gain contrasts 0.01740/0.08999/0.16685. In the explicit feedback-delay scan, 0–1 min did not meet the transient criterion, 2 min produced only 0.00888 suppression, and 4 min was the first scanned value above the frozen threshold (0.03610). Thus this model needs a measurable delayed feedback path, not instantaneous feedback, but these are model-defined delays, not biological estimates.

Verification: `python3 -m unittest -v test_model.py` and compilation passed. Next step is timestamped simultaneous CRH/ACTH/free-cortisol/total-cortisol measurements around a known intervention, followed by compartment-aware fitting.
