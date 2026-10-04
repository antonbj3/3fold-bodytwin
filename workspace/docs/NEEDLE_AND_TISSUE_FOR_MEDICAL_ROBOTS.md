# Needle and tissue mechanics for medical robots

Can a robot recover needle-tip load from shaft measurements and predict how tissue will resist it?
Shaft force alone cannot identify tip load, and adding moment gives a conditional interval only when contact positions and sensor-error bounds are known.

The measured needle pull-off mean lies outside both friction-only capacity intervals. The skin-tearing calculation also fails its mechanism comparison, so neither result supplies a joint prediction of penetration, tearing and grasp. [5, 6]

## First and repeated needle passages

The recorded post-penetration measurements show different friction on first and second passage, where cutting is absent. Treating every passage as a repeated passage over-allows the count on some budgets, but changes none of the selected budget decisions. [1]

| Quantity or test | Chain or recorded result | Comparator and consequence | Source |
|---|---|---|---|
| Post-penetration friction | First passage **0.349 N**, scatter **0.018 N**; second **0.241 N**, scatter **0.004 N** | Equal passage friction misses the difference | [1] |
| First/repeated friction ratio | **1.44813**, propagated scatter **0.07846** | Describes this measurement, not a depth-independent material constant | [1] |
| Selected friction budgets | **0** changed decisions among **9** tested budgets | No demonstrated count advantage on this selection; includes a refusal | [1] |
| Recorded budget scan | **739 / 1670** returned counts differ from equal friction | The control over-allows a passage on those budgets; scan increment **0.001 N** | [1] |
| Budget refusal | **0.3 N** returns no count | Below the recorded lower dispersion endpoint **0.331 N** | [1] |

The budget sums friction forces across passages. It is a cumulative score in force units, not peak instantaneous force or mechanical work, and assumes the measured second-passage friction for every later passage. Neither that assumption nor a depth-independent ratio follows from the recorded measurements. [1]

Primary friction figures are absent from the supplied directories. The values are traceable to the stored output, but cannot be re-extracted from the figures here. Reported scatter is neither a deterministic bound nor a confidence interval. [1]

## What force and moment identify

Tip load follows `y = L + S` and `M = lL + rS`, giving `L = (M − ry)/(l − r)`. Here `S` is other contact resistance, and the arms are projected moment arms. A friction-coefficient comparison has the same force-only ambiguity. [2]

| Quantity or test | Chain or recorded result | Comparator and consequence | Source |
|---|---|---|---|
| Force-only ambiguity | Identical shaft force admits tip loads **0.2–0.8 N** | Force alone does not identify tip load | [2] |
| Force and moment, nearby nuisance contact | At **1.0 N**, **0.02 N·m**, tip arm **0.15 m**, nuisance arm **0.01 m**: tip load **[0.06357, 0.07929] N** | Conditional interval under the declared geometry and sensor errors | [2] |
| Same readings, changed nuisance contact | Nuisance arm **0.075 m**: tip load **[−0.75667, −0.71000] N** | Even the sign changes; contact location is an input | [2] |
| Sensor-error assumptions | Moment **0.001 N·m**, force **0.01 N** | Declared bounds, not measured device performance | [2] |
| Arm-separation refusal | Fraction below **0.1** | Equal arms make the inversion singular | [2] |
| Boundary example | Stored fraction **0.09999999999999991**, refused | Accompanying note gives **0.09999999999999999** instead | [2] |

Changing the nuisance-contact position can reverse the recovered load's sign at identical shaft readings. The intervals propagate declared sensor bounds with fixed arms, and do not cover uncertainty in contact position or isolate pure adhesion or cutting force. [2]

There is no transferable pressure calibration for adhesion dominance in the supplied contact results. Separating tip load from other resistance requires an independent constraint on that resistance. [2]

## Grasp capacity and needle pull-off

The flat-handle interval uses an envelope of published material means. Its empirical status is unknown, and it does not describe rotational wrench, peak pressure or rate-dependent friction. [6]

| Quantity or test | Chain or recorded result | Comparator and consequence | Source |
|---|---|---|---|
| Flat handle/glove grasp | Total tangential capacity **[11.6, 13.2] N** at **10 N** per finger | Mathematical interval using an envelope of published material means; empirical status unknown | [6] |
| Needle-holder pull-off | Friction-only **[12, 16] N**; alternative two-surface convention **[24, 32] N** | Recorded measured mean **20.11 N** lies outside both; neither closure covers it | [6] |

The held-out needle pull-off mean defeats both force conventions. Without matched jaw normal loads and slip measurements, their disagreement does not identify whether curved-needle interlocking explains the failure. [6]

## Bending and tearing

The bending ratio agrees with an analytic reference under imposed bending displacement. That boundary condition helps enforce the agreement, so the comparison does not establish the ratio under moment or transverse loading. [4]

| Quantity or test | Chain or recorded result | Comparator and consequence | Source |
|---|---|---|---|
| Bending/stretch response | Analytic **0.25**; solver **0.2501331541804101** at span/thickness **20.0**, **0.2501679265348504** at **2.0** | Agreement under imposed bending kinematics; not a test of unrestricted loading | [4] |
| Skin opening/shear toughness ratio | Fixed inventory **1.0**; predicted pooled ratio **1.0238068702894416** | Recorded observed ratio **1.474757281553398**; mechanism test fails | [5] |
| Conditional skin shear toughness | **14.32683469588746 kJ/m²**, relative error **0.3045225875782787** | Conditional fixed-inventory model fails; conventional comparison is a tie | [5] |

A mechanical solver requires a volume mesh, elastic constants, boundary displacement fields, a bending definition, a field response functional and units. The tissue's lateral contraction ratio and the specimen's end conditions are missing. [4]

The skin result fails despite numerical agreement checks, and the conventional comparison is a tie. Connectivity, fibre recruitment, terminal opening and anchor survival remain unknown. Original assay rows and mechanism code are unavailable here, so the recorded toughness values support no joint tearing-and-penetration prediction. [5]

## Instrument uncertainty: an optical comparison only

Paired optical measurements illustrate how device disagreement and covariance affect a margin. They do not measure needle-sensor error or calibrate a robot instrument. [3]

| Quantity or test | Chain or recorded result | Comparator and consequence | Source |
|---|---|---|---|
| Paired-device corneal thickness | **89** eyes; mean difference **9.89721 µm**, SD **6.25897 µm** | A device offset and spread in an optical dataset; no needle transfer established | [3] |
| Optical instrument-component margin | **69** eyes; covariance-aware **1.6893**, quadrature **1.8227**, summed-spread **1.3221** | Dimensionless margins against **0.5 D**; output flags a heavy tail | [3] |

## References

1. Recorded post-penetration force source, passage-count output and decision script.
2. Shaft force/moment decision output and script; recorded friction-coefficient comparator.
3. Paired-device optical instrument-component output and computation script.
4. Finite-element bending/stretch output, computation script and recorded boundary-condition limitation.
5. Skin toughness comparison output; original assay rows and mechanism code unavailable here.
6. Robot contact interval records, exact grasp computation, published-input extraction and needle pull-off check.

## Where the numbers come from

Each entry maps the displayed numbers to their source file and key. Rounded endpoints
are distinguished from raw source values and exact rational endpoints.

- [1] **0.349, 0.018; 0.241, 0.004 N**: `results/ASSEMBLY_PASSAGE_FRICTION/DECISION_V1.json :: first_passage_N; second_passage_N`.
- [1] **1.44813; 0.07846**: `results/ASSEMBLY_PASSAGE_FRICTION/DECISION_V1.json :: passage_friction_ratio; passage_friction_ratio_sd` (rounded).
- [1] **0 / 9; 739 / 1670; 0.001 N**: `results/ASSEMBLY_PASSAGE_FRICTION/DECISION_V1.json :: budgets_where_the_decision_changes; budgets_tested; budget_axis_scan.budgets_where_the_decision_differs_from_the_control; budget_axis_scan.budgets_with_a_decision; budget_axis_scan.step_N`.
- [1] **0.3; 0.331 N**: `results/ASSEMBLY_PASSAGE_FRICTION/DECISION_V1.json :: rows[0].budget_N; rows[0].reason`.
- [2] **0.2–0.8 N**: `results/ASSEMBLY_TIP_LOAD/DECISION_V1.json :: why_force_alone_fails`.
- [2] **1.0 N; 0.02 N·m; 0.15, 0.01 m; [0.06357, 0.07929] N**: `results/ASSEMBLY_TIP_LOAD/DECISION_V1.json :: rows[0].shaft_force_N; rows[0].shaft_moment_Nm; rows[0].tip_arm_m; rows[0].nuisance_arm_m; rows[0].tip_load_enclosure_N` (interval rounded).
- [2] **0.075 m; [−0.75667, −0.71000] N**: `results/ASSEMBLY_TIP_LOAD/DECISION_V1.json :: rows[1].nuisance_arm_m; rows[1].tip_load_enclosure_N` (interval rounded).
- [2] **0.001 N·m; 0.01 N**: `tasks/assembly/tip_load_decision.py :: decide(moment_noise_Nm, force_noise_N)` defaults `1e-3`, `1e-2`.
- [2] **0.1; 0.09999999999999991; 0.09999999999999999**: `results/ASSEMBLY_TIP_LOAD/DECISION_V1.json :: refusal.min_arm_separation_fraction; rows[2].arm_separation_fraction; boundary_note`.
- [3] **89; 9.89721, 6.25897 µm**: `results/ASSEMBLY_INSTRUMENT_COMPONENT/INSTRUMENT_V1.json :: length_quantities_both_devices[0].eyes_with_both_devices; length_quantities_both_devices[0].mean_difference; length_quantities_both_devices[0].sd_of_difference`.
- [3] **69; 1.6893; 1.8227; 1.3221; 0.5 D**: `results/ASSEMBLY_INSTRUMENT_COMPONENT/INSTRUMENT_V1.json :: covariance_aware_margin.n_units; covariance_aware_margin.margin; covariance_aware_margin.rss_margin; covariance_aware_margin.sum_margin; spec_D`.
- [4] **0.25; 0.2501331541804101; 20.0; 0.2501679265348504; 2.0**: `results/ASSEMBLY_INTEROP_BENDING/INTEROP_V1.json :: ours_analytic; rows[0].external_solver_ratio; rows[0].span_over_thickness; rows[3].external_solver_ratio; rows[3].span_over_thickness`.
- [5] **1.0; 1.0238068702894416; 1.474757281553398; 14.32683469588746; 0.3045225875782787**: `results/LANE_SKIN_TOUGHNESS_GAP/night_rounds/r2.json :: outcome.metrics.fixed_inventory_I_III_ratio; outcome.metrics.predicted_pool_I_III_ratio; outcome.metrics.observed_I_III_ratio; outcome.metrics.conditional_modeIII_mean_kJ_m2; outcome.metrics.conditional_modeIII_relative_error`.
- [6] **[11.6, 13.2]; 10 N**: `results/LANE_ROBOT_PORT_CONNECT/HUM_ROBOTPORT_20261003/PORT_HUM_20261003.json :: ports[1].computed_value.{lower,upper}; ports[1].certificate.payload.bounds.normal_force_per_finger[0].{lower,upper}` (capacity stored as `58/5`, `66/5`).
- [6] **[12, 16]; [24, 32]; 20.11 N**: `results/LANE_ROBOT_PORT_CONNECT/HUM_ROBOTPORT_20261003/needle_external_check.json :: paper_friction_only_capacity_N.{lower,upper}; two_surfaces40N_per_jaw_capacity_N.{lower,upper}; external_mean_pull_off_force_N` (mean stored as `2011/100`; also `results/LANE_ROBOT_PORT_CONNECT/HUM_ROBOTPORT_20261003/published_inputs.json :: needle_driver.held_mean_pull_off_N`).
- Constraint discrepancies and status: `data/CONSTRAINT_NET_TISSUE.json :: bodytwin.tissue_constraint_net.edges` (passage relation: `constraint`, `third_correction_2026_10_04`; skin relation: `status`, `evidence`). Bending limitation: `results/ASSEMBLY_INTEROP_BENDING/INTEROP_V1.json :: limitation_found_in_my_own_design`.

Status: pending independent review. No claim of clinical validation.
