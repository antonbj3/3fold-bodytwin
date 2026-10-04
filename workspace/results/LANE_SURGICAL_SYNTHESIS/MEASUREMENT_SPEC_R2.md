# Quantitative measurement specification, R2

PENDING_INDEPENDENT_REVIEW. This is a computed **planning analysis of unchanged synthetic closures**. Native strength uncertainty, full hemostasis time and edge necrosis are UNKNOWN/null. The table states prediction errors for model proxies; it is no biological posterior or achieved laboratory precision. No acquisition has been ordered.

256 frozen Sobol scenarios with seed6122 are propagated through the entire chain. [Design and precision budgets](MEASUREMENT_DESIGN_R2.json) state all intervals and costs. Biological and perfusion widths share a synthetic donor driver; other drivers are independent **only as a planning assumption**. Measurement error has separate random error and shared calibration error. Optimal linear prediction uses `V_after=V_y−C_yz(C_zz+R)^−1 C_zy`; interpretation as posterior variance requires a Gaussian approximation to these moments. We also show an estimator trial trained on128 scenarios and evaluated on128 new ones.

R1's six items remain. M1/M2/M3 take already represented physical port observations. M5 assumes in a conditional information scenario that chemistry/direction/connected crossings and mechanics can identify the synthetic Q observer, bridge and a constant reference scale. That assay mapping is missing in the world. M4 lacks linkage from chemistry/collagen moles to U/I/M and M6 from work to biological/perfusion damage; therefore their zero values mean only zero **implemented** downstream information value. The value when the unknown laws are identified is UNKNOWN.

Cost means chosen relative work units for preparation, instruments and analysis, with0,5–2× sensitivity. Actual hours/prices are null. The score is the equally weighted mean reduction of variance in strength28/90, min(hemostasis time,120) and damage proxy width, divided by cost. With a constant censored outcome its sampling variance is counted as0 without giving the unknown late time variance0.

| Planning rank | R1 item | Cost | SD strength28 before→after pp | SD strength90 before→after pp | SD min(T,120) before→after min | SD damage proxy before→after µm | Mean variance reduction/cost |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | M2: Oxygen/storage/boundary flux with fixed reference volume | 3 | 2.842→2.707 | 6.730→6.665 | 0.000→0.000 | 196.680→116.851 | 0.06325 |
| 2 | M1: Registered gap, independent injury/perfusion widths, radius/density/pressure/flow | 8 | 2.842→2.790 | 6.730→6.668 | 0.000→0.000 | 196.680→89.356 | 0.02649 |
| 3 | M5: Chemistry-age-connected bridge and independent traction/ultimate-load port | 16 | 2.842→1.279 | 6.730→4.265 | 0.000→0.000 | 196.680→188.730 | 0.02305 |
| 4 | M3: Early radius-matched plug and persistent-seal dynamics | 5 | 2.842→2.835 | 6.730→6.674 | 0.000→0.000 | 196.680→190.210 | 0.00430 |
| 5 | M4: Absolute precursor/HP chemistry and reference collagen (missing observer map) | 12 | 2.842→2.842 | 6.730→6.730 | 0.000→0.000 | 196.680→196.680 | 0.00000 |
| 6 | M6: Full force/work/new-area/interface and terminal-work ledger | 20 | 2.842→2.842 | 6.730→6.730 | 0.000→0.000 | 196.680→196.680 | 0.00000 |

256/256 scenarios do not reach hemostasis before120min and are censored; the closed fraction is0.0000. Sampling-SD for the bounded time is0.000min. This does not identify later hemostasis: the requested unbounded time's SD and information value remain null. Therefore measure r/Q/plug on a longer measured clock before a new time prior. Observations at120min alone cannot distinguish different late courses.

The damage proxy is the Q036 translation's width where `1−exp(−hazard)>0,20`; the hypoxia proxy uses10Torr. These declared thresholds are no necrosis assay. Native edge necrosis needs separate viability with a typed cell, time and local registration. For all six items its posterior-SD and SD reduction are null.

## Precision per subobservation

σ states a chosen future random-error target; “shared” is a calibration error shared across a package. These are requirements for the calculation, not verified assay accuracy. All subitems from R1 that lack an observation operator are in the design's `unrepresented` and are commented on above; they have UNKNOWN native information value.

| Package | Observation | σ / shared | ΔSD strength28 pp | ΔSD strength90 pp | ΔSD damage proxy µm |
|---|---|---:|---:|---:|---:|
| M1 | gap | 1e-05 / 0 | 0.0205 | 0.0363 | 0.0095 |
| M1 | biological_width | 2.5e-05 / 1e-05 | 0.0130 | 0.0138 | 52.7041 |
| M1 | perfusion_width | 2.5e-05 / 1e-05 | 0.0242 | 0.0181 | 102.9823 |
| M1 | radius_factor | 0.05 / 0.02 | 0.0003 | 0.0004 | 0.0904 |
| M1 | density_factor | 0.1 / 0.02 | 0.0001 | 0.0005 | 0.0054 |
| M1 | driving_pressure | 200 / 80 | 0.0042 | 0.0003 | 1.6240 |
| M1 | blood_viscosity | 0.000175 / 0 | 0.0011 | 0.0001 | 0.0518 |
| M1 | vascular_path | 5e-05 / 0 | 0.0012 | 0.0056 | 0.0300 |
| M2 | pressure20 | 2 / 1 | 0.0916 | 0.0096 | 20.1915 |
| M2 | pressure100 | 2 / 1 | 0.1141 | 0.0132 | 25.8553 |
| M2 | pressure300 | 2 / 1 | 0.1126 | 0.0017 | 63.1566 |
| M2 | face_flux | 0.02 / 0.005 | 0.0106 | 0.0002 | 4.8548 |
| M2 | oxygen_storage | 6 / 0 | 0.0005 | 0.0014 | 0.0344 |
| M2 | face_pressure | 2 / 1 | 0.0035 | 0.0000 | 1.7889 |
| M3 | coverage8_small | 0.03 / 0.01 | 0.0007 | 0.0036 | 1.3481 |
| M3 | coverage8_large | 0.03 / 0.01 | 0.0034 | 0.0005 | 1.6135 |
| M3 | seal8 | 0.05 / 0.02 | 0.0007 | 0.0105 | 5.3483 |
| M3 | initial_flow | 1e-10 / 0 | 0.0001 | 0.0005 | 0.3520 |
| M3 | aggregation_factor | 0.1 / 0 | 0.0001 | 0.0000 | 0.1975 |
| M3 | seal_rate | 2 / 0 | 0.0016 | 0.0014 | 0.0171 |
| M4 | chemistry_D10 | 0.005 / 0.002 | 0.0000 | 0.0000 | 0.0000 |
| M4 | chemistry_H10 | 0.002 / 0.001 | 0.0000 | 0.0000 | 0.0000 |
| M4 | chemistry_rate | 0.007 / 0 | 0.0000 | 0.0000 | 0.0000 |
| M4 | chemistry_yield | 0.03 / 0 | 0.0000 | 0.0000 | 0.0000 |
| M5 | trace_M21 | 0.015 / 0.005 | 0.6860 | 0.6533 | 7.2128 |
| M5 | bridge_fraction | 0.05 / 0 | 0.2409 | 0.7298 | 0.0101 |
| M5 | strength_reference_factor | 0.05 / 0.02 | 0.1996 | 0.6169 | 0.0289 |
| M5 | strength21 | 1 / 0.5 | 1.2168 | 1.6383 | 4.4708 |
| M6 | gamma_cut | 30 / 0 | 0.0000 | 0.0000 | 0.0000 |

SI units and the observations' exact meaning: gap/widths/path in m, pressure in Pa (driving) or Torr (pressure20/100/300/face), μ in Pa·s, r/density-factor and coverage/seal/bridge/reference/trace_M dimensionless, initial_flow m³/s, face_flux mlO2/(m² min), oxygen_storage mlO2/(m³ Torr), rates1/min or1/dygn according to the port table, D/H mol/molreferenskollagen, strength21 pp. M5's mechanical reference is an explicitly synthetic port acquisition at21dygn, not a late strength-target fit.

## Sensitivity, dependencies and useful prioritization

- n256_noise1.0_floorFalse: M2 > M1 > M5 > M3 > M4 > M6.
- n128_noise1.0_floorFalse: M2 > M1 > M5 > M3 > M4 > M6.
- n256_noise0.5_floorFalse: M2 > M1 > M5 > M3 > M4 > M6.
- n256_noise2.0_floorFalse: M2 > M1 > M5 > M3 > M4 > M6.
- n256_noise1.0_floorTrue: M2 > M1 > M3 > M5 > M4 > M6.

M2: new128-scenario trial gives variance/MSE reduction28/90/damage proxy 9.0/2.2/65.9%.

M1: new128-scenario trial gives variance/MSE reduction28/90/damage proxy 1.2/1.0/81.1%.

M5: new128-scenario trial gives variance/MSE reduction28/90/damage proxy 81.1/60.3/3.2%.

M3: new128-scenario trial gives variance/MSE reduction28/90/damage proxy -1.2/-2.5/2.8%.

M4: new128-scenario trial gives variance/MSE reduction28/90/damage proxy 0.0/0.0/0.0%.

M6: new128-scenario trial gives variance/MSE reduction28/90/damage proxy 0.0/0.0/0.0%.

Precision0,5×/2× and128/256 compare the same frozen prior. Separate cost intervals can change rank: for each pair dominance is counted only when the candidate's lower value/cost exceeds the control's upper one. The additional floors20pp and50µm are open sensitivity assumptions for an unresolved observation law, no measured error bounds. They can neither be identified nor reduced by the modeled port measurements.

The concrete next acquisition shall combine an amount inventory on a fixed reference with the early measurement ports that carry the planning table. Native strength still requires M4's absolute chemistry/reference **and** M5's connected mechanical observations with the same cohort. A low model score for M4 is an indicator of a missing consumer, not a reason to forgo chemistry. M6 is needed for new tools when the transfer law is filled. The parent goal's empirical prioritization is therefore still conditional.

All raw scenarios, subobservations, covariance/noise and costs are in [MEASUREMENT_VALUE.json](r2/measurement_v1/MEASUREMENT_VALUE.json) and sample000…255.json. Large field data are on `external_mount`. The strongest control receives the same acquisition, moments and optimal linear estimator: **TIE**.

## Executed nonlinear information trial after linear loss

The preserved128→128 trial worsened hypoxia prediction MSE by68,16% forM2 and32,09% forM1. The moment table's sampling-based contraction therefore cannot alone guarantee predictive benefit for this threshold observable. We changed the observation calculation to a positive Gaussian kernel with training-based, measurement-noised leave-one-out over five frozen bandwidths. Physics, biology and measurement budgets are unchanged.64 new Sobol scenarios (seed6124), each with32 independent future measurement-noise draws, evaluate prediction risk. Training uses only the earlier256 scenarios. The initial CV seed version is preserved; training6126 and future measurement noise6125 are separated in V2 before evaluation.

Here RMS is **prediction error on new scenarios**, which differs from the earlier Gaussian moment-SD. No native posterior intervals follow.

| Package | New64: RMS strength28 pp | RMS strength90 pp | RMS damage proxy µm | RMS hypoxia proxy µm | Hypoxia gate against prior | Strength/damage against linear control |
|---|---:|---:|---:|---:|---|---|
| M1 | 2.709 | 5.439 | 121.984 | 50.752 | PASS | FAIL |
| M2 | 2.669 | 5.446 | 132.339 | 44.018 | PASS | PASS |
| M3 | 2.683 | 5.252 | 193.687 | 56.689 | FAIL | PASS |
| M4 | 2.681 | 5.237 | 196.253 | 56.540 | PASS | PASS |
| M5 | 1.326 | 3.991 | 185.901 | 55.573 | PASS | PASS |
| M6 | 2.681 | 5.237 | 196.253 | 56.540 | PASS | PASS |

Frozen pivot gate: **FAIL**. The new observed value/cost order is M2 > M5 > M1 > M3 > M4 > M6. All original linear misses, the initial training, V2 coefficients, new64 fields and measurement noise/residuals remain in [VALIDATION.json](r2/nonlinear_information_v1/VALIDATION.json) and its raw/model files. The strongest conventional kernel control has the same input and operation: TIE. This is a trial of an information operator on synthetic data, no new biological capability.
