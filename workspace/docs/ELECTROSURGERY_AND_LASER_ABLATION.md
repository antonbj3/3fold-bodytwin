# Electrosurgical incision and laser ablation

Can a damage-depth tolerance determine electrosurgical dwell time, or a removal target determine laser fluence?
The calculations return conditional settings, but the dwell interval misses published damage means and the laser rule loses to a fixed fluence at **400 µg/cm²**.

Neither result supplies a patient-specific setting. The incision calculation lacks setting-dependent injury calibration, and the laser comparison measures agreement with a published law. [1, 2]

## Electrosurgical dwell and damage depth

Microscopy of fresh-frozen human meniscus supplies the damage-depth means. They distinguish apparatus and exposure time, but not individual output settings. [1]

| Incision quantity | Calculated answer | Check against the source or calculation |
|---|---|---|
| Mean damage depth, bipolar apparatus A, at **1 / 3 s** | — | **564 / 648 µm** [1] |
| Mean damage depth, bipolar apparatus B, at **1 / 3 s** | — | **1444 / 1697 µm** [1] |
| Mean damage depth, monopolar apparatus C, at **1 / 3 s** | — | **895 / 1057 µm** [1] |
| Apparent diffusivity, depth squared divided by time | **3.1360e-7–9.6333e-7 m²/s** | All source means give **1.39968e-7–2.085136e-6 m²/s**; the stored interval misses **2** means |
| Damage tolerance **0.5 mm**, setting **30 or 60** | Dwell interval **0.2595–0.7972 s**; returns lower endpoint | Same answer at both settings; no setting-dependent coefficient |
| Damage tolerance **1.0 mm**, setting **30 or 60** | Dwell interval **1.0381–3.1888 s**; returns lower endpoint | Bounds inherit the incomplete apparent-diffusivity interval |
| Damage tolerance **1.5 mm**, setting **30 or 60** | Dwell interval **2.3356–7.1747 s**; returns lower endpoint | No measured injury guarantee |
| Setting **90** | No answer | Outside the accepted setting range |

The dwell bounds are tolerance squared divided by apparent diffusivity, with the lower endpoint returned. They enclose parameter choices rather than statistical confidence, and the diffusivity interval fails to cover all the source means. [1]

A damage boundary does not measure thermal conductivity or identify an evaporative mechanism. The energy costs below also depend on assumed material properties, rather than measured cutting work.

| Assumed cut geometry | Vaporisation cost | Fracture cost |
|---|---|---|
| **1.0 mm²** kerf | **2260 J/m** | **4 J/m** |

The accepted setting range is an implementation limit, not an established injury boundary. Its supporting correlations mix coagulation depth and cut volume across operating modes, so they do not demonstrate deterioration of a single observable at high settings.

| Setting range | Stored correlations | Consequence |
|---|---|---|
| **10–60** | **0.95 / 0.98 / 0.92** | Accepted by the calculation |
| **10–120** | **0.73** | No answer outside the accepted range |

A setting-dependent correction can reproduce the same depth curve as a change in mechanism. A slope break therefore cannot identify that mechanism. The conditional setting–duration comparison below tests separability, provided its observation-error bound holds independently.

| Setting–duration comparison | Conditional criterion |
|---|---|
| Absolute log depth cross-ratio | Reject separability above **0.08** |
| Bounded log-depth error per observation | **0.02** |

## Laser fluence and removed mass

The reference is a held-out mass-loss regression in porcine reticular dermis. Its threshold is the regression intercept, rather than the first local boiling event. The Gaussian and top-hat beam conventions give different calibration values, with no measured beam map to decide between them. [2]

| Laser calibration | Published law | Gaussian beam proxy | Top-hat diagnostic |
|---|---|---|---|
| Regression threshold, J/cm² | **1.15** [2] | **1.4548** | **1.9081** |
| Regression slope, µg/J | **267** | **200.33** | **325.46** |

Fluence equals threshold plus target mass per area divided by slope. The Gaussian calculation gives the choices below. “Removed mass” is the published law's prediction at that fluence, not a measured tissue outcome.

| Target, µg/cm² | Chosen fluence, J/cm² | Removed mass under published law, µg/cm² | Mass error | Fixed-fluence mass error |
|---|---|---|---|---|
| **25** | **1.5796** | **114.71** | **+358.85%** | **+1550.00%** |
| **100** | **1.9540** | **214.67** | **+114.67%** | **+312.50%** |
| **400** | **3.4516** | **614.51** | **+53.63%** | **+3.13%** |
| **800** | **5.4483** | **1147.64** | **+43.45%** | **−48.44%** |

The fixed control uses **2.6949 J/cm²**. It is closer at the target identified in the opening, although the calculated rule is closer at the other listed targets. Using the published law directly matches that law by construction, so this score cannot establish an advantage over that control.

Transfer to another pulse, tissue preparation or device is unsupported without the raw forcing vector and measured beam profile.

## Power and temperature are different readouts

Nominal generator settings do not equal the sampled calculated power. The surface-temperature comparison also fails, and its source current and cohort temperatures are unpaired.

| Source or readout question | Computed result | Limitation or failed comparison |
|---|---|---|
| Nominal electrical setting **25 / 50 W** | Sampled model power **15.8980–21.9038 / 33.2991–45.5367 W** | Calculated instantaneous power; not measured time-mean delivery |
| Surface-temperature transfer, baseline | **8 of 24** graphical reading boxes overlap; maximum centre error **32.2066 K** | Failed comparison; source current and cohort temperatures are unpaired |
| Same **10 J** and **4 s**, reversed heating order | At **1 mm** depth: **31.8681 / 27.6955 °C** | Analytical counterexample: energy and duration do not determine depth temperature |

Reversing the heating order changes depth temperature despite equal energy and duration. This is an analytical heat solution with fixed material and source area, no spatial grid and no requirement that surface temperatures match. It is not a wound-temperature measurement.

Absolute temperature at injury depth is unavailable. The depth-readout archive reports failed bench calibration without the underlying calibration outputs. Surface temperature, depth temperature and histological damage remain separate quantities.

## Focal delivery and membrane transport

The receptor averaging length comes from measured membrane hop scales and assumes diffusion in the membrane plane. It is a transport length, not an optical spot size. [3]

| Source or readout question | Computed result | Limitation or failed comparison |
|---|---|---|
| Patterned ligand delivery, **10 s** pulse | Receptor averaging length **2.6480–3.1305 µm** | Conditional transport length from measured hop scales [3], not a laser spot size |
| Patch spacing **1 / 3 / 5 µm** | Reject / cannot decide / necessary geometry condition passes | Passing does not establish an activation advantage over a uniform bath |

Patch spacing alone does not establish an activation advantage over a uniform bath at equal total dose. Ligand gradients and receptor response lack an observation model, and focal delivery does not certify optical focus.

## References

1. [Radiofrequency thermal effects on the human meniscus: an in vitro study of systems with monopolar and bipolar electrodes](https://pubmed.ncbi.nlm.nih.gov/12642261/). Apparatus-specific microscopy means.
2. [Comparison of pulsed carbon-dioxide laser ablation at different wavelengths](https://pubmed.ncbi.nlm.nih.gov/9694144/). Mass-loss regression; the slope is transcribed in the decision output.
3. [Rapid hop diffusion of a G-protein-coupled receptor in the plasma membrane as revealed by single-molecule techniques](https://pubmed.ncbi.nlm.nih.gov/15681644/). Membrane confinement and residence times.

## Where the numbers come from

Paths are relative to the repository root. Displayed decimals are rounded; intervals retain both endpoints.
- Kerf **1.0**, energy costs **2260 / 4**: `results/ASSEMBLY_INCISION_SETTING/DECISION_V1.json :: energy_balance_recomputed_by_coordinator.{kerf_mm2,vaporisation_J_per_m,fracture_J_per_m}`.
- Exposure **1 / 3**, means **564 / 648 / 1444 / 1697 / 895 / 1057**: `results/LANE_SETTING_TO_DIFFUSIVITY/r1/sources/PMID12642261_EPMC_V2.json :: resultList.result[0].abstractText`; extracted as `results/LANE_SETTING_TO_DIFFUSIVITY/r1/OUTCOME_V1.json :: parent_challenge.held_rows[*].{time_s,mean_depth_um}`. Apparatus A/B/C follow that row order.
- Original diffusivity **3.1360e-7–9.6333e-7**: `results/ASSEMBLY_INCISION_SETTING/DECISION_V1.json :: facit.d_eff_m2_s`; all-means **1.39968e-7–2.085136e-6**, missed count **2**: `results/LANE_SETTING_TO_DIFFUSIVITY/r1/OUTCOME_V1.json :: parent_challenge.{source_descriptive_D_extrema_m2_s,parent_missed_means}`.
- Tolerances **0.5 / 1.0 / 1.5**, settings **30 / 60 / 90**, dwell bounds **0.2595–0.7972 / 1.0381–3.1888 / 2.3356–7.1747**: `results/ASSEMBLY_INCISION_SETTING/DECISION_V1.json :: rows[*].{tolerance_mm,setting,max_dwell_s_conservative,max_dwell_s_optimistic}`.
- Gate **10–60**, reported range **10–120**, correlations **0.95 / 0.98 / 0.92 / 0.73**: `results/ASSEMBLY_INCISION_SETTING/DECISION_V1.json :: refusal.{calibrated_setting_range,reported_setting_range,R_restricted_depth_radius_volume,R_full_range_volume}`; scope correction: `data/CONSTRAINT_NET_TISSUE.json :: bodytwin.tissue_constraint_net.edges[161].correlation_scope_corrected` and `tasks/assembly/incision_setting_decision.py :: module documentation`.
- Thresholds **1.15 / 1.4548 / 1.9081**, slopes **267 / 200.33 / 325.46**: `results/ASSEMBLY_LASER_FLUENCE_DECISION/LASER_DECISION_V1.json :: external_anchor.{threshold_J_cm2,slope_ug_per_J}` and `decision_rows[0,6].{our_threshold_J_cm2,our_slope_ug_per_J}`.
- Targets **25 / 100 / 400 / 800**, chosen fluences **1.5796 / 1.9540 / 3.4516 / 5.4483**, masses **114.71 / 214.67 / 614.51 / 1147.64**, errors **358.85 / 114.67 / 53.63 / 43.45%**: `results/ASSEMBLY_LASER_FLUENCE_DECISION/LASER_DECISION_V1.json :: decision_rows[0,2,4,5].{target_ug_cm2,fluence_we_would_choose_J_cm2,mass_actually_removed_under_published_law_ug_cm2,mass_error_percent}`.
- Fixed fluence **2.6949**, errors **1550.00 / 312.50 / 3.13 / −48.44%**: `results/ASSEMBLY_LASER_FLUENCE_DECISION/LASER_DECISION_V1.json :: control_fixed_fluence.fluence_J_cm2` and `control_fixed_fluence.rows[0,2,4,5].mass_error_percent`.
- Settings **25 / 50**, model power **15.8980–21.9038 / 33.2991–45.5367**: `results/LANE_SETTING_TO_DIFFUSIVITY/r6/OUTCOME_FINAL_V1.json :: source_model_power_W.{25,50}`. Sampling scope: `results/LANE_SETTING_TO_DIFFUSIVITY/r5/WAVEFORM_INPUTS_V1.json :: {use,averaging_class}`; these extrema are not temporal averages.
- Overlap **8 of 24**, error **32.2066**: `results/LANE_SETTING_TO_DIFFUSIVITY/r8/OUTCOME_V1.json :: held_transfer.baseline.{overlap_count,count,max_center_residual_K}`; the boxes are graphical reading allowances, not confidence intervals.
- Energy **10**, duration **4**, depth **1**, temperatures **31.8681 / 27.6955**: `results/LANE_SETTING_TO_DIFFUSIVITY/r1/OUTCOME_V1.json :: energy_duration_sufficiency.pairs[0,1].{energy_J,duration_s,temperatures_C["1"]}`; analytical solution: `results/LANE_SETTING_TO_DIFFUSIVITY/r1/thermal_query.py :: temperature`.
- Pulse **10**, averaging bounds **2.6480–3.1305**, spacing **1 / 3 / 5**: `results/ASSEMBLY_FOCAL_DELIVERY/decision.json :: gate_table[2].{pulse_duration_s,averaging_length_um}` and `worked_verdicts[0,2,1].inputs.patch_spacing_um`.
- Depth-readout withholding: `results/LANE_DEPTH_TEMPERATURE/night_rounds/r4.json :: {gate,gate_scope,files}`; the listed bench-calibration outputs are absent from this archive. Their numerical summaries are not used here.
- Rejection threshold **0.08**, error bound **0.02**: `results/ASSEMBLY_INCISION_SETTING/DECISION_V1.json :: witness_thresholds.{four_corner_rejects_above,log_depth_error_per_observation_eta}`; replacement test: `data/CONSTRAINT_NET_TISSUE.json :: bodytwin.tissue_constraint_net.edges[161].replacement_witness`.

Status: pending independent review. No claim of clinical validation.
