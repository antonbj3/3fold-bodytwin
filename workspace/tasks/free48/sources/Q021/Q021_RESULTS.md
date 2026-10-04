BT-HX-Q021

## Results

There is no previous model or result file in the working directory; only `BRIEF.md`, `inputs/` and `agent.log` existed. `PREREG.md` was created before the first run and frozen with SHA-256 `607e6318dc411d055dd00c875389acbc0e1c9ed0062bbfbd33567c4ed1a2f9e4` in `PREREG.sha256`.

### Built model

`model.py` has nine regions and couples at a fixed point:

- perfusion: `Q_i = (P_art,i - P_ven,i) / R_v,i`, with hydrostatic pressure gradient and HPV as increased hypoxic resistance;
- ventilation: `V_i` proportional to `C_i * driving_pressure / R_aw,i`;
- O2 content: Hill saturation curve with Hb binding and dissolved O2;
- regional O2 mass balance: `Q_i(C_a,i-C_v) = (1000/760) V_A,i(P_IO2-P_AO2,i)`;
- `P_AO2_global` is ventilation-weighted regional mass-balance PO2. The independent alveolar gas equation is reported as a consistency check.

The parameter table with units and source/assumption is in `model.py:57`; the unit check is in `model.py:405`. Analytical limiting cases and placebo control are in `test_model.py:24`, `test_model.py:40` and `test_model.py:54`.

### Model output

All model numbers below come from `results.json` and are derived, not measurement data.

| Quantity | Regional rest | Regional heavy load | Homogeneous rest | Homogeneous heavy load |
|---|---:|---:|---:|---:|
| A-aDO2 (mmHg) | 6.667433 | 2.441848 | ~0 | ~0 |
| Regional minus homogeneous A-aDO2 (mmHg) | 6.667433 | 2.441848 | — | — |
| Nondependent/dependent perfusion | 0.321475 | 0.540052 | 1.0 | 1.0 |
| Mean V/Q: nondependent / middle / dependent | 2.081978 / 0.954189 / 0.817655 | 1.282036 / 0.965413 / 0.894357 | 1 / 1 / 1 | 1 / 1 / 1 |
| O2 uptake (ml/min) | 363.253339 | 2904.002057 | 372.049305 | 2924.531350 |

The model thus shows a clear regional effect at rest (6.667433 mmHg difference against the homogeneous model), but a smaller regional effect after the assumed load (2.441848 mmHg). Recruitment and bronchodilation flatten the gradients in this parameter case. It is a result of the model’s assumptions, not a claim about normal exercise physiology.

### Preregistered references

- Hall et al., DOI `10.1152/japplphysiol.00659.2013`, Table 4: perfusion `2.9 ± 1.7` (nondependent), `4.5 ± 1.7` (middle), `4.2 ± 1.7` (dependent) at rest; `4.0 ± 1.5`, `5.7 ± 2.1`, `5.1 ± 1.6` during exercise, all in `ml min^-1 ml^-1`. Derived nondependent/dependent ratios are `0.690` and `0.784`.
- Tedjasaputra et al., DOI `10.1152/japplphysiol.00778.2012`, Table 3: A-aDO2 `6.3 ± 3.7 Torr` at rest and `23.3 ± 5.3 Torr` during exercise.
- Harf et al., DOI `10.1152/jappl.1978.44.1.115`, indexed abstract: `40–150%` increase of apical blood flow at 50 W.

These values are included in `results.json` under `reference_anchors`; they were not used to calibrate the model.

### Frozen criteria and sensitivity

Of the five frozen criteria, only the placebo control passed. The heavy load did not give at least `10 mmHg` increase, did not fall in `10–40 mmHg`, did not achieve `5 mmHg` regional effect and did not give the perfusion ratio `0.784 ± 0.20`. No parameters were changed to make a criterion pass.

±50% sensitivity at heavy load (`results.json:sensitivity`):

- `hpv_gain` 0.5/1.0/1.5x: A-aDO2 `3.100682 / 2.441848 / 2.033177` mmHg.
- `rv_base_gradient` 0.5/1.0/1.5x: `3.254240 / 2.441848 / 1.852845` mmHg.
- `aw_resistance_gradient` 0.5/1.0/1.5x: `2.953539 / 2.441848 / 1.990285` mmHg.

### Verification and next steps

`python3 -m unittest -v test_model.py` gave 5/5 PASS; `python3 model.py` wrote `results.json`. All four case checks and the placebo difference `2.27e-13` passed. The first run revealed that the independent gas equation did not agree with homogeneous mass balance; the final version therefore uses the ventilation-weighted mass balance and reports the gas equation separately. The gap is `20.226557` mmHg at rest and `29.476793` mmHg under heavy load.

