# bodytwin_core

Run from the workspace root with Python 3.13, NumPy, SciPy and pytest:

```bash
OPENBLAS_NUM_THREADS=1 python -m bodytwin_core.demo
OPENBLAS_NUM_THREADS=1 python -m pytest -q bodytwin_core/tests
```



## Empirical femur geometry priors

`bodytwin_core.geometry` exports `prior_from_anthropometry(height_m, weight_kg, sex)`, `contralateral(observed, side, parameter)`, `prior_from_tibia(tibial_length_mm, plateau_width_mm)`, and `scalar_population_prior()`. The first three return a `ScalarPrior` with `mean`, interval endpoints, units, source, and an `in_range` indicator. They initialize missing geometry; they do not overwrite measured anatomy. Inputs outside the training covariate range are marked as extrapolations.

| Prior | Source and held-out error | Returned interval | Validity range |
|---|---|---|---|
| Femur length from height, weight, sex | VSD, 26 complete persons, person-held-out MAE 9.78 mm; median baseline 24.23 mm | 90% interval ±26.91 mm, calibrated on person-held-out absolute residuals | Height 1.52–1.87 m; weight 37.4–90 kg; sex F/M; predicts person mean length |
| Contralateral head radius | VSD, 30 paired persons, person-held-out MAE 0.533 mm; direct-copy MAE 0.524 mm | 95% paired-difference interval ±3.24 mm; held-out coverage 58/60 (96.7%) | Observed radius 19.06–27.14 mm; includes flagged fits |
| Contralateral CCD | VSD paired differences; Imperial-to-VSD location transfer MAE 2.34° | VSD 95% paired-difference interval ±9.96°; source sample is small | Observed CCD 122.39–141.27° |
| Head radius from tibial length and plateau width | Imperial, 35 persons/70 femurs, person-held-out two-input OLS RMSE 1.042 mm | 95% interval ±2.22 mm from held-out residuals | Tibial length 297.05–384.66 mm; plateau width 52.32–78.21 mm; no external VSD tibia check |
| Scalar population Gaussian (radius, CCD, signed anteversion) | Imperial fit; 51/54 unflagged VSD femurs inside the independently recalibrated 95% ellipsoid | Mahalanobis squared cutoff 11.432 from Imperial subject folds | Imperial observed ranges: radius 17.79–26.78 mm, CCD 120.29–139.44°, signed anteversion −38.94–38.69°; population distribution, not an individual prediction |

The DA-030 *incremental* result used 28 additional femur variables: Imperial person-held-out RMSE was 0.929 mm without tibia and 0.741 mm with tibia. `prior_from_tibia` is a simpler two-input model, with its own held-out error above; it does not inherit the 0.741 mm claim. The Gaussian off-diagonal covariance did not clearly beat a diagonal model on VSD. The Imperial bilateral difference gives useful radius and CCD locations on VSD, but its intervals did not transfer; the API uses VSD-calibrated widths. There is no anteversion transfer prior (DA-028), no contralateral absolute-anteversion prior on VSD (DA-029), and no priors from DA-031, DA-032, DA-036 or DA-037. Calibration and transfer details are in `results/CX-GEOMPRIOR/RESULTS.md`.
- `solver`: CX-SOLVER2 `cx_solver.solve_step` uses the N43 inner solve from `run_hybrid.py` (A295). It reached 141/141 KKT ≤ 1e-10 and 2,583.270 N right-hip peak on the identity the collaborator lift. `solve_k1(..., mode="certified")` applies XF4's active-set emulator with a numerical KKT/dual-gap certificate and exact fallback on its K1 scope. `solve_lift(..., mode="certified")` currently takes exact fallback on all the collaborator contact frames because XF4 was not trained on that system. The default is `exact`; see `results/CX-EMUCERT/RESULTS.md` for coverage and timing. `CX-SOLVER-ROBUST` reports scaling robustness in its result package.
- `solver.solve_bounds(A, b, joint_maps, offsets, lower, upper, A_ub, b_ub)` returns per-frame HiGHS min/max for each **affine force component** under equality, inequality and box constraints. It preserves infinite endpoints and raises on infeasibility. A three-dimensional force norm needs a separate nonlinear or conic optimization; the component intervals must not be read as exact norm bounds. See `results/CX-SETVALUED/`.
- `whatif`: CX-WHATIF2 `implicit_kkt.frame` differentiates a fixed active set after ±1° rebuilds (A369). The archived identity CCD peak is 8.351990 N/° versus 8.352019 N/° full finite difference. The local ±1° matrices are incomplete, so the demo labels archived derivatives; at active-set changes it selects the nearest stable frame. CX-D1PARITY (A907) verified Field U380 derivative parity at 3.0e-10–2.6e-8 relative error; `whatif.field_u380_shape_massprop` calls Field code at its source path, with the SHA-256 in `SOURCES.json`.
- `emulator`: BT-XF4 predicts within each active set (A303); archived maximum relative error is 0.146% on its K1 test set. Its model is scoped to that set and is not used as a the collaborator-lift emulator.
- `population`: N2b batched solver (A363) and BT-N50 band rule. N2b's 1,000-person, 9.5 mm attachment ensemble gives right-hip peak −8.49/+15.88%; N50's calibrated factor is 0.95857. The demo transfers N2b's relative band as a reference; an individual N50 band requires posterior-draw SD.
- `nullmodels`: B24 activity-specific bodyweight peak, N1 `k·|GRF|`, and N1g activity-group `k·|GRF|`. B24 hip/gen1 Walking LOSO RMSE is 70.374640 %BW. the collaborator lift lacks a fitted N1g activity coefficient, and these knee/gait or spine/lift scopes are not observed the collaborator hip-force truth.
- `determined`: `knee_load(trial, system=...)` returns per-frame inverse-dynamics inputs, the parameter-free law and its fixed lever-arm sensitivity envelope, exact affine total/medial/lateral contact LP intervals when an L1 system is supplied, plus measured/N1g/selected-model comparison and the local gait2392 sensor-closure protocol. `joint_load(system)` calls `solve_bounds`; infeasible frames are NaN by default and can be made strict with `on_infeasible='raise'`. Previously computed bounds can instead be passed to `knee_load(..., bounds=...)` after independent frame-alignment checks. `python -m bodytwin_core.demo --determined` prints all Grand Challenge trial rows and writes `results/CX-DETERMINED/determined.png`. See `results/CX-DETERMINED/RESULTS.md` for scope and held thresholds.

The demo times instantiation, all 141 solve frames, archived what-if lookup, band and null-model comparison. It does not claim clinical accuracy. Knee-force muscle models lost to N1g at 0.378 BW RMSE (A312/A324/A358). The DXA strength-law claim was corrected: incremental gain over mass, height, sex and age is 2.34%, below its 5% criterion (A359). Geometry-to-full-body force coupling beyond the frozen the collaborator systems remains unverified.

## Conditional design loop

`bodytwin_core.design.design_eval({'lengthening_mm': L}, subject='z001', task='collaborator_lift')` rebuilds Field femur anchors, N40b mechanics and projection, U380 bone mass properties, and the 141-frame N43 recruitment before reporting hip reaction. `L` is bounded to ±10 mm. `C01RFE` uses a TLEM proxy anatomy, not a same-person mapping. New intermediate files go to `external_media`; the pipeline uses frozen D1 and CX-WHATIF2 inputs, with source hashes in `results/CX-DESIGNLOOP/SOURCES.sha256`. See that result package for the exact preregistration, feasibility checks, timing, parity and limitations. The knee contact layer is outside this hip design call.

## L1 contact band batch

`bodytwin_core.contact_band_batch.contact_band_batch(trials, eps=None, tol=.01)` accepts L1 gait trial keys and returns frame-aligned lower/upper bands in BW. `eps=None` selects LP; a nonnegative epsilon selects the constrained minimum-effort SOCP band. For trial keys, the mask is the frozen CX-INVERSEOC mask. NumPy, SciPy and Clarabel are required; `highspy` enables the faster persistent LP basis path and SciPy HiGHS is the fallback. Pose sensitivity is diagnostic only: the function returns an endpoint from an identical problem, a numerically KKT-checked active set, or an exact solver. The complete benchmark, negative ≥20× speed result and frozen-band comparisons are in `results/CX-POSECACHE/RESULTS.md`.

Wrap evaluation requires an independently supplied surface evaluator and otherwise raises `NotImplementedError`. The landmark baseline reads `hip_joint_mm` from the caller-provided `geometry_inputs.json`; the reference H5 input is selected by `BODYTWIN_REFERENCE_H5`. These inputs are not bundled.
