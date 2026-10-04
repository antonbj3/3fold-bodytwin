BT-HX-Q107

## Outcome
The runnable model in `model.py` couples the fundus and antrum to wall pressure, a 3/min slow-wave, a time-dependent pyloric opening and separate mobilisation/size selection of solid particles. The run used a 1 s internal step and a 5 s saved time series in `results.json`.

- **Liquid:** full geometry gives `T50 = 16.38 min`, `30.89 mL` remaining at 120 min and `7.72 %` residual. Reference: `20 min`; criterion `[10, 40] min`: **PASS**.
- **Solid:** full geometry gives `T50 = 79.16 min`, `36.61 g` remaining at 120 min and `36.61 %` residual. Reference: `90 min`; criterion `[45, 180] min`: **PASS**.
- **Conservation:** maximum relative volume/mass deviation `2.57e-15` (`≤1e-8`): **PASS**. All six tests in `test_model.py` pass, including the orifice-law limiting case and the analytical exponential null model.

## Geometry and sensitivity
The scalar null model gives `T50 = 23.10 min` (liquid) and `92.42 min` (solid). Full geometry therefore does not improve the T50 error against this single literature anchor; the improvement requirement on held-out curves is therefore **UNKNOWN**, not passed. Exploratorily, the solid residual at 120 min is `36.61 %` versus `40.66 %` for the null model, but this is no measurement validation.

The largest ±50 % effect on T50 is the pyloric radius (`2.88` maximum relative change), the pyloric opening range (`0.840`) and particle mobilisation (`0.284`). Geometry can thus improve liquid curves when it changes the time-dependent opening; solid food requires measurement of retained particles. No such data are available here.

## Provenance and next steps
**Source:** Collins, P. J. (1983), *Gut*, `T50 ≈ 20/90 min`, **UNVERIFIED, from memory**. **Derivation:** SI equations, unit checks and experimental assumptions in `model.py`; all numbers in `results.json`. **Hypothesis:** geometry helps when pressure–gating and particle retention govern; not yet tested on independent data.

The next step is synchronised measurements of volume, pressure, flow and retention curves for both phases, while geometry, particle size distribution and the tissue's viscoelastic properties are measured and calibrated on held-out people.
