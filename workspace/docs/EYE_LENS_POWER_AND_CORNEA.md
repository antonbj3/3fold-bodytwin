# Intraocular lens power and corneal shape

Which lens power best reaches the target refraction, and does measuring the posterior cornea improve astigmatism prediction?
In this retrospective cohort, power choices from preoperative inputs have a mean miss of **0.5345 D**, below that of implanted powers against a hindsight reference, and measured posterior cornea improves the recorded cylinder scores.

The power reference depends on the optical calculation itself. It does not measure the outcome of an alternative implant, and the comparison does not establish prospective clinical performance. [1, 2]

The astigmatism advantage depends on the axis convention. Using another device's printed posterior axis gives a worse result than the population estimate. Thickness and central curvature also fail as general summaries of corneal optics. [4, 8]

## Lens power and postoperative depth

Postoperative lens depth is predicted from anterior chamber depth, axial length, lens thickness and central corneal thickness. Each eye is excluded from its own linear fit. The retrospective power calculation instead uses measured postoperative depth and removes bias using the other patients' outcomes. [1, 2]

The hindsight power combines measured postoperative spherical-equivalent refraction with the calculation's local sensitivity to implant power. MAE means mean absolute error, and D denotes dioptres. [1, 2]

| Test | Chain or measured-input result | Comparator or limitation | Source |
|---|---|---|---|
| Retrospective power, 20 eyes | Power MAE 0.6883 D; closer in 11 eyes | Implanted power MAE 1.0699 D; closer in 2; postoperative depth supplied | [1] |
| Preoperative-input power, 89 eyes | Power MAE 0.5345 D; closer in 36 eyes | Implanted power MAE 1.1304 D; closer in 22 | [2] |
| Measured postoperative depth, same power cohort | Power MAE 0.569 D | Predicted-depth result above is smaller, without establishing why | [2] |
| Predicted postoperative depth | Reported MAE 0.114 mm; maximum error 0.5647 mm | Measured depth's sample standard deviation 0.2842 mm | [2] |
| Recommendation changes when depth is predicted | 44 of 89 grid choices change | This compares predicted-depth with measured-depth recommendations, not with implanted powers | [2] |
| Power within the script's 0.25 D threshold | 29 of 89 eyes | Threshold is applied to power miss, not to measured refraction error | [2] |
| Refraction target and power grid | Target 0.0 D | Recommendations rounded to a 0.5 D grid | [1, 2] |

Power inversion returns no finite recommendation when the target lies outside its search bracket. Implant labels are treated as exact, while axial-length device bias and actual implant material and geometry remain unresolved. [1, 2, 8]

## Posterior cornea and astigmatism

Cylinder scores compare predictions with postoperative cylinder, using either magnitude alone or magnitude with axis. The population control is a fixed posterior-cylinder estimate, with the orientations shown below. It is not a fitted clinical nomogram. [3]

| Test | Chain or measured-input result | Comparator or limitation | Source |
|---|---|---|---|
| Cylinder magnitude, 69 toric eyes | Measured posterior MAE 0.2808 D | Population estimate MAE 0.6209 D; alternative estimate axis 0.3636 D | [3] |
| Cylinder including axis, same toric eyes | Measured posterior error 0.4985 D | Population estimate 0.8923 D; alternative estimate axis 0.6051 D | [3] |
| Posterior cylinder from another device | Aligned-axis MAE 0.2932 D | Printed-axis MAE 0.7525 D: worse than the population estimate | [4] |
| Total corneal cylinder magnitude | Assembled mean 2.0291 D; measured mean 1.8314 D | Mean absolute difference 0.2918 D | [5] |
| Measured total cornea substituted | Printed-axis MAE 0.3112 D; rotated-axis MAE 0.3109 D | Assembled-cornea MAE 0.2808 D; mean spherical power still supplied by the chain | [5] |
| Total corneal axis agreement | Output median gap 86.45° | Script description says 3.55°; the output rows support the larger gap | [5] |
| Population-control orientation | Posterior cylinder 0.30 D at 90° | Alternative orientation 0° | [3] |

The device comparisons use the same eyes and postoperative outcomes. Similar residuals after rotating the total-cornea axis do not establish axis agreement, because the assumed implant axis rotates with it. The output and its accompanying description disagree about the axis gap. [3, 4, 5]

These scores assume implant alignment. They omit incision-induced astigmatism, subsequent implant rotation and postoperative corneal shape change. The scalar power calculation does not answer an astigmatic planning question. [1, 3]

## Can predictions replace preoperative measurements?

Each biometric prediction beats its cohort-mean control but misses by more than the devices disagree. The lens-depth result therefore does not justify replacing these preoperative measurements. [2, 6]

Instrument disagreement is the mean absolute difference between devices on the same 89 eyes, not error against a known true value. The predictions use leave-one-out fits from the remaining biometric quantities, axial length and anterior keratometry. [6]

| Preoperative quantity | Instrument disagreement | Prediction MAE | Leave-one-out cohort-mean MAE | Source |
|---|---|---|---|---|
| Central corneal thickness | 9.91484 µm | 28.51378 µm | 30.70271 µm | [6] |
| Anterior chamber depth | 0.15484 mm | 0.17782 mm | 0.32525 mm | [6] |
| Lens thickness | 0.04084 mm | 0.28498 mm | 0.3877 mm | [6] |
| Corneal white-to-white diameter | 0.16164 mm | 0.24368 mm | 0.33375 mm | [6] |

## What corneal summaries leave out

Thickness, index and central curvature do not determine the optical effect of a surface. The recorded counterexamples concern synthetic geometry, with clinical direction unknown. [7, 8]

| Test | Chain or measured-input result | Comparator or limitation | Source |
|---|---|---|---|
| Thickness, index and surface-allocation summary | General sufficiency fails in synthetic surface counterexamples | Clinical direction remains unknown | [8] |
| Central curvature as a wavefront summary | Sufficiency fails when higher surface terms change | This is a geometric counterexample, not a patient outcome | [8] |
| Curvature-fit and aperture summaries | General sufficiency is unresolved; sampled counterexamples are reported | A pupil-specific sampled check does not establish general sufficiency | [7] |

Surface position, deformation laws and pupil-resolved geometry remain unspecified. The primary surface arrays are absent from this snapshot, so these failures support no quantitative patient prediction. [7, 8]

The power ranking remains conditional on its hindsight sensitivity, and the posterior-measurement advantage remains conditional on the population control and axis convention. Neither comparison supplies the separate-cohort result, registered surfaces or actual implant alignment needed to settle those interpretations. [1–5, 8]

## References

1. Stored retrospective implant-power decision and its scoring script.
2. Stored preoperative-input implant-power decision and its scoring script.
3. Stored toric predictions, magnitude and axis-sensitive scores, and their calculation script.
4. Stored posterior-cylinder device substitution and its calculation script.
5. Stored total-cornea comparison and its calculation script; their axis accounts disagree.
6. Stored biometric prediction-versus-device-disagreement test and its calculation script.
7. Recorded corneal-fit and pupil-support tests; primary arrays are unavailable here.
8. Recorded optical surface counterexamples and domain constraint edges; primary surface arrays are unavailable here.

## Where the numbers come from

Keys below are relative to each JSON root; script keys name constants or the function/docstring.

- [1] `results/ASSEMBLY_IOL_DECISION/DECISION_V1.json :: summary.{patients,mean_power_miss_twin_D,mean_power_miss_implanted_D,twin_closer_to_target,surgeon_closer_to_target}` → 20, 0.6883, 1.0699, 11, 2, respectively; `summary.target_refraction_D` → 0.0.
- [1, 2] `tasks/assembly/iol_decision.py :: POWER_STEP`; `tasks/assembly/prospective_power.py :: GRID_D` → 0.5; the latter's `TARGET_D` → 0.0 and `TOL_D`, used by `main` for power miss, → 0.25.
- [2] `results/ASSEMBLY_PROSPECTIVE_POWER/PROSPECTIVE_V1.json :: eyes` → 89; `scored_against_hindsight_correct_power.{mean_miss_with_predicted_position_D,mean_miss_of_implanted_power_D,prospective_closer_than_implanted,implanted_closer_than_prospective,mean_miss_with_measured_position_D,within_tolerance_prospective}` → 0.5345, 1.1304, 36, 22, 0.569, 29, respectively.
- [2] `results/ASSEMBLY_PROSPECTIVE_POWER/PROSPECTIVE_V1.json :: lens_position_prediction.{mean_abs_error_mm,max_abs_error_mm,spread_of_the_quantity_itself_mm}` → 0.114, 0.5647, 0.2842, respectively; `power_recommendation_shift_D.eyes_where_the_grid_choice_changes` → 44.
- [3] `results/ASSEMBLY_TORIC_DECISION/TORIC_V1.json :: summary.toric_patients` → 69; `summary.{residual_cyl_mae_with_measured_posterior_D,residual_cyl_mae_with_population_estimate_D,residual_cyl_mae_with_population_estimate_alternative_axis_D}` → 0.2808, 0.6209, 0.3636, respectively.
- [3] `results/ASSEMBLY_TORIC_DECISION/TORIC_V1.json :: summary.{residual_vector_mae_with_measured_posterior_D,residual_vector_mae_with_population_estimate_D,residual_vector_mae_with_population_estimate_alternative_axis_D}` → 0.4985, 0.8923, 0.6051, respectively.
- [3] `tasks/assembly/toric_decision.py :: {POP_POSTERIOR_CYL_D,POP_POSTERIOR_AXIS_DEG,POP_POSTERIOR_AXIS_ALTERNATIVE_DEG}` → 0.30, 90, 0, respectively; `to_vector` and `main` define the scoring and assumed alignment.
- [4] `results/ASSEMBLY_TORIC_CROSS_DEVICE/CROSS_DEVICE_V1.json :: residual_cylinder_mae_D.{casia_posterior,casia_posterior_as_printed}` → 0.2932, 0.7525, respectively; `tasks/assembly/toric_cross_device.py :: main` defines the held inputs.
- [5] `results/ASSEMBLY_TOTAL_CORNEA_CHECK/TOTAL_CORNEA_V1.json :: cylinder_magnitudes_D.{assembled_mean,measured_total_mean,mean_abs_difference}` → 2.0291, 1.8314, 0.2918, respectively.
- [5] `results/ASSEMBLY_TOTAL_CORNEA_CHECK/TOTAL_CORNEA_V1.json :: residual_cylinder_mae_D.{measured_total_as_printed,measured_total_rotated,assembled_cornea}` → 0.3112, 0.3109, 0.2808, respectively; `tasks/assembly/total_cornea_check.py :: main.variants` retains assembled spherical power and aligns the implant to each substituted axis.
- [5] `results/ASSEMBLY_TOTAL_CORNEA_CHECK/TOTAL_CORNEA_V1.json :: axis_relationship_tested_not_assumed.median_axis_gap_deg` and `rows[*].axis_gap_deg` → 86.45; `tasks/assembly/total_cornea_check.py :: module docstring` → 3.55. These conflict.
- [6] `results/ASSEMBLY_PREDICT_VS_MEASURE/PREDICT_VS_MEASURE_V1.json :: results[0].eyes` → 89; each following triple uses `{measurement_disagrees_with_itself_mean,prediction_error_from_other_quantities_mean,cohort_mean_control_error}`.
- [6] `results/ASSEMBLY_PREDICT_VS_MEASURE/PREDICT_VS_MEASURE_V1.json :: results[0]` → 9.91484, 28.51378, 30.70271; `results[1]` → 0.15484, 0.17782, 0.32525; `results[2]` → 0.04084, 0.28498, 0.3877; `results[3]` → 0.16164, 0.24368, 0.33375, respectively.
- [6] `tasks/assembly/predict_vs_measure.py :: PAIRS, EXTRA, main` defines quantities, units, predictors and mean control.
- [7] `results/LANE_CORNEA_SHAPE/night_rounds/r2.json :: subgates, obstacle, files` records sufficiency failures, unresolved scope and unavailable primary artifacts; no numerical performance claim is taken from this summary.
- [8] `results/LANE_EYE_OPTICAL_TWIN/night_rounds/r13.json :: outcome.surface_sufficiency, obstacle`; `results/LANE_EYE_OPTICAL_TWIN/night_rounds/r14.json :: outcome.eta_sufficiency, outcome.curvature_sufficiency, clinical_gate, files` records synthetic failures and unknown clinical direction; no numerical performance claim is taken from these summaries.
- [8] `data/CONSTRAINT_NET_TISSUE.json :: bodytwin.tissue_constraint_net.edges[6,19,20,21]` supplies domain limitations; the grid-change interpretation in edge `[21]` conflicts with the comparator implemented in reference [2].

Status: pending independent review. No claim of clinical validation.
