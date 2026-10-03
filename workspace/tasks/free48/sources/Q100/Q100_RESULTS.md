BT-HX-Q100

## Slutsats

A first runnable reduced model is built in `model.py`: two lateral masses, homogeneous or series-connected tissue, unidirectional contact, capillary/mucus film and a compressed Bernoulli/resistive air stream. All Parameters are with unit and source/assumption in `results.json` and dimensional checks are `PASS`.

## Prediction against reference

| Variant | f0 (Hz) | Fel versus 120 Hz (Hz) | Kontaktandel (%) | Mean wall pressure (Pa) | Mean flow (m^3/s) |
|---|---:|---:|---:|---:|---:|
| baseline | 302.148 | 182.148 | 30.792 | 798.919 | 2.60577e-5 |
| layered | 138.910 | 18.910 | 14.221 | 788.429 | 7.60034e-5 |

The reference is about 120 Hz from N. Ishizaka and J. L. Flanagan, *Two-mass model of the vocal folds*, *Journal of the Acoustical Society of America* (1972). **UNVERIFIED, FROM MEMORY**; geometry and pressure matching are not verified. The frozen wide criterion is 60–240 Hz plus at least 5 um peak-to-peak movement. The layered model yields 138.910 Hz and 99.283 um: `PASS`. The improvement flag is `PASS` because the absolute f0 error decreases from 182.148 to 18.910 Hz. This is a broad mechanistic endorsement, not empirical validation; reference match and contact/mucus data is `UNKNOWN`.

Other layered outputs are glottal area 1.63570e-5 to 1.69527e-5 m^2, maximum contact indentation 14.339 um, wall pressure 778.798–800 Pa, maximum flow 9.94819e-5 m^3/s and declared mucus wave velocity proxy 0.948683 m/s. The proxy formula is `sqrt(T_sheet/(rho_m*h_m))`, not a measured wave velocity. The energy terms are air work 0.024448 mW, tissue loss 0.356229 mJ, contact work 0.015764 mJ, meniscal work -0.003950 mJ and final kinetic energy 0.100913 uJ.

## Sensitivity, ±50 %

| Parameter | -50 %: f0 / kontakt / amplitud | +50 %: f0 / kontakt / amplitud |
|---|---|---|
| `k_surface` | 127.875 Hz / 15.346 % / 232.974 um | 131.090 Hz / 0.000 % / 2.635 um |
| `k_contact_surface` | 137.331 Hz / 17.871 % / 106.332 um | 139.851 Hz / 11.022 % / 95.420 um |
| `p_sub` | 273.563 Hz / 55.711 % / 62.208 um | 137.246 Hz / 12.772 % / 292.468 um |

`p_sub` controls f0 and contact phase strongest. `k_surface` mainly controls the boundary between open oscillation and near stationary contact. `k_contact_surface` mainly affects contact proportion and less f0. No parametric sensitivity is calibrated against measurement data.

## Verification and builds on

`PREREG.md` was created before first run and `PREREG.sha256` verifies. `python3 test_model.py` passes the analytical undamped linear limit, placebo similarity (`placebo` == `baseline`) and dimensional checks. No previous model files or results were in this directory; only `BRIEF.md`, empty `ALLOW_WEB` and `inputs/` were present at inspection. `inputs/QUESTION.md` is scope specification, not measurement data. No internet source, external dataset or measured geometry was used.

## Next resolution step

1. Measure glottal geometry and layers tjocklek/fiberorientering in the same trial sequence.
2. Synchronize high-speed video with subglottal pressure, flow and acoustics over a pressure sweep.
3. Measure contact area, film thickness, contact angle and mucus rheology; replace the proxy parameters with identified material parameters.
4. Then compare the layered and two-mass model on the same recordings and calibrate only with train/validation-separerade data.

This is a mechanistic first model, not a claim about the phonation of saturated subjects. `results.json` contains all numbers, derivation and source status.
