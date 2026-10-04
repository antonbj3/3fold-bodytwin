# Spine and joint load: pressure, motion and force sharing

What load can a disc or joint carry without exceeding a pressure limit?
These calculations give conditional pressure and displacement ranges, but the available force, pressure and geometry measurements do not establish a permissible load for a person.

The disc-pressure prediction exceeds the walking measurement, with no matched disc force to identify the cause. The finite-element meniscus areas and peak pressures also cannot carry the stated whole-contact load. Those mismatches prevent interpreting either calculation as a load limit. [1, 6]

## Disc pressure and axial load

The pressure law is `P = kF/A`. Its walking force is a lumped lumbar contact resultant from another subject, rather than an independently measured force through the observed disc. The implanted transducer's walking range describes means around maxima and minima, rather than instantaneous extrema. [1]

| Question or input | Chain or calculation | Comparator and interpretation |
|---|---|---|
| Disc pressure inputs | Area **1800 mm²**; nominal stress **0.9045 MPa**; multiplier **1.3–1.5** | Stored reference load **1628.0 N**; explicit replay **1628.1 N**. The displayed load difference is rounding. |
| Pressure at the explicit replay load | **1.17585–1.35675 MPa** | Measured walking range **0.53–0.65 MPa**: above-range arithmetic; attribution to multiplier error is unresolved. |
| Load giving any overlap with that pressure range | **636–900 N** | Inverse compatibility set under the retained area and multiplier; no measured walking force is supplied. |
| Axial load at a **0.5 MPa** pressure target | **600–692 N** | **1252–1536 N** using the pressure/reference-stress quotient; that quotient assumes a matched force which is absent. |
| Axial load at a **0.65 MPa** pressure target | **780–900 N** | **1628–1997 N** under the same unestablished assumption. These are conditional inversions, not tolerable loads. |

Pressure alone cannot distinguish an error in force from an error in the multiplier. A walking label does not resolve the missing disc-force measurement, segment mapping or matched pressure statistic. The recorded comparison of multiplier derivations is withdrawn, with no admitted replacement.

## Disc height and facet load

The displacement calculation combines ex-vivo disc stiffness with a separately sourced facet-load fraction. Their common spinal level, posture and specimen are unestablished. Stiffness dispersion supplies a sensitivity range, not a statistical confidence interval. [2, 3]

| Question or input | Chain or calculation | Comparator and interpretation |
|---|---|---|
| Disc/facet inputs | Disc stiffness **1734 ± 446 N/mm**; facet fraction **0.03–0.25** | The stiffness dispersion is used as a sensitivity range, not a statistical confidence interval. |
| Disc height loss at **200 N** total load | **0.0688–0.1506 mm** | Varying facet fraction alone gives width **0.0254 mm**; stiffness alone **0.0546 mm**. |
| Disc height loss at **1000 N** total load | **0.344–0.7531 mm** | Sensitivity interval from the same load split; no matched displacement measurement scores it. |

Stiffness contributes more interval width than facet fraction in this calculation. The output's recommendation to prioritise facet fraction contradicts its own numerical sensitivities. Neither interval has a matched displacement measurement.

## Ligament limits are sampled, not measured motion limits

The tension lookup compares all represented ligaments with the anterior longitudinal ligament (ALL) alone. The supraspinous ligament (SSL) limits flexion in the case where the intervals differ, while ALL limits extension.

| Question or input | Chain or calculation | Comparator and interpretation |
|---|---|---|
| Ligament tension ceiling **3 N** | Lower sampled endpoint **−10°**; upper endpoint unresolved | No complete sampled interval around neutral is returned. |
| Ligament tension ceiling **5 N** | All represented ligaments: **[−30°, +10°]**, width **40°** | ALL alone: **[−90°, +10°]**, width **100°**; overlooks **60°** of flexion restriction. SSL binds flexion, ALL extension. |
| Ligament tension ceiling **7 N** | **[−90°, +10°]** | Same sampled interval using ALL alone; no additional restriction from the other ligaments at this ceiling. |
| Angular sampling | **20°** spacing; no neutral sample | Values between samples are not certified |

The force-per-strain coefficients are chosen values, not measured spring stiffnesses. Unmatched ligament slack and attachment registration prevent interpreting the sampled endpoints as a measured range of motion.

The work-sign scenario has no measured anatomical threshold. A disc-boundary response also cannot identify forces in individual ligaments.

## Muscle-force redistribution

Uniform scaling does not reproduce the static-optimisation force redistribution. The comparator is another calculation, not a human force measurement.

| Question or input | Chain or calculation | Comparator and interpretation |
|---|---|---|
| Muscle-force scaling | **187** shared muscles; mean error **≈4.0 N**, maximum **≈45.3 N** | Explains fraction **≈0.702** of the redistribution. The output rejects uniform scaling as sufficient. |
| Repeated muscle-force rows | **10** rows per arm; **1** effective pose per arm | Repeating the score across rows supplies no evidence across different poses. |
| Empirical load/work prediction | **0 of 92** held cases jointly pass a chronological linear response | Held disc-boundary response, not individual ligament forces; uniform empirical prediction fails. |

Repeated rows represent the same effective pose, so they do not establish performance across poses. The failed load/work prediction concerns disc-boundary response and supplies no individual-ligament force law.

## Meniscus area, pressure and force must describe the same contact

For compressive contact over a common area, `load = integral(pressure dA) ≤ peak pressure × area`. The products below are capacity bounds, not measurements of local force or load share.

| Question or input | Chain or calculation | Comparator and interpretation |
|---|---|---|
| Meniscus intact versus removed, cadaver | At **1000 N**, area **1150 → 520 mm²**, peak **3 → 6 MPa** | Peak/mean factors **3.45** and **3.12** describe these observations; they do not establish an invariant shape. |
| Intact contact, finite-element source | At **1000 N**, **110 mm² × 1.2 MPa** gives capacity **132 N** | Capacity/load **0.132**: the reported area and peak cannot describe the complete contact carrying the stated load. |
| Partial radial tear, finite-element source | **80 mm² × 2.1 MPa**, capacity **168 N** | Capacity/load **0.168** under the same stated load; same contact-domain failure. |
| Complete radial tear, finite-element source | **35 mm² × 3.3 MPa**, capacity **115.5 N** | Capacity/load **0.1155**; same contact-domain failure. |
| Repaired contact, finite-element source | **95 mm² × 1.4 MPa**, capacity **133 N** | Capacity/load **0.133**; same contact-domain failure. |
| Constant concentration after partial meniscectomy | Area ratio **0.9** predicts peak ratio **≈1.111** | Reported peak ratio **1.65** [5]. Mean-pressure-only control also predicts **≈1.111**; the added concentration assumption earns no improvement. |
| Pressure/area/load label sweep | **8 of 12** pressure-and-force edges have capacity below stated load | Counts are edges, including paired area/pressure records, rather than independent experiments. |

The finite-element source states the axial load itself. Its area and peak-pressure values therefore need a resolved contact domain before they can describe how that load is carried. Reported dispersions do not establish paired hard bounds, and the pressure definition and local reaction force remain unresolved. [6]

Cadaver peak/mean factors describe their own observations, without establishing a constant contact shape. The partial-resection comparison fails to transfer that assumption and gains nothing over the mean-pressure-only control. Its approximate source summaries lack paired dispersion. [4, 5]

## References

1. *New in vivo measurements of pressures in the intervertebral disc in daily life*. [Primary paper](https://doi.org/10.1097/00007632-199904150-00005). Retained full text supplies pressure, area and statistic.
2. *Comparison of animal discs used in disc research to human lumbar disc: axial compression mechanics and glycosaminoglycan content*. [Publication record](https://pubmed.ncbi.nlm.nih.gov/18344845/). Stiffness is from the retained extraction, not a fresh reading of the primary table.
3. *Mechanism of facet load transmission as a hypothesis for low-back pain*. [Publication record](https://pubmed.ncbi.nlm.nih.gov/6238423/). Fraction is from the retained extraction; transfer to the disc-stiffness specimens is unestablished.
4. *The contact area and pressure distribution pattern of the knee. A study of normal and osteoarthrotic knee joints*. [Publication record](https://pubmed.ncbi.nlm.nih.gov/6894212/). Retained abstract and extracted data.
5. *Meniscal tears: the effect of meniscectomy and of repair on intraarticular contact areas and stress in the human knee. A preliminary report*. [Publication record](https://pubmed.ncbi.nlm.nih.gov/3755296/). Retained abstract; changes are approximate.
6. *Radial tears of the lateral meniscus reduce contact area by 70% and triple condylar stress: A physics-based finite element analysis*. [Publication record](https://pubmed.ncbi.nlm.nih.gov/41539441/). Retained abstract states the applied load and reported contact quantities.

## Where the numbers come from

- Disc inputs **1800, 0.9045, 1.3–1.5, 1628.0**: `results/ASSEMBLY_DISC_LOAD_DECISION/DISC_LOAD_V1.json :: measured_inputs.{disc_area_mm2,nominal_stress_MPa,our_multiplier_range,reference_load_N}`.
- Replay **1628.1**, pressure **1.17585–1.35675**, overlap **636–900**: `results/LANE_DISC_LOAD_LABEL/r1/OUTPUT_MODEL_REPLAY_V1.json :: {load_n,predicted_pressure_mpa.interval,compatibility_load_n.interval}`; rounding operation: `tasks/assembly/disc_load_decision.py :: main/reference_load_n and measured_inputs.reference_load_N`.
- Walking pressure **0.53–0.65**: `results/LANE_DISC_LOAD_LABEL/r1/SOURCE_LOAD_V1.json :: walking_pressure_MPa`; primary check: `results/LANE_DISC_LOAD_LABEL/r1/WILKE_1999_V1.txt :: Results, walking paragraph`.
- Target **0.5**, loads **600–692, 1252–1536**: `results/ASSEMBLY_DISC_LOAD_DECISION/DISC_LOAD_V1.json :: decision_rows[0].{pressure_limit_MPa,load_our_model_N,load_facit_implied_N}`.
- Target **0.65**, loads **780–900, 1628–1997**: `results/ASSEMBLY_DISC_LOAD_DECISION/DISC_LOAD_V1.json :: decision_rows[1].{pressure_limit_MPa,load_our_model_N,load_facit_implied_N}`.
- Ceiling **3**, endpoint **−10**: `results/ASSEMBLY_LIGAMENT_ROM/LIGAMENT_ROM_V1.json :: decision_rows[3].{tension_limit_N,all_five_ligaments.low_deg}`; missing upper endpoint: `all_five_ligaments.high_deg`.
- Ceiling **5**, angles **−30, +10, −90**, widths **40, 100**, difference **60**: `results/ASSEMBLY_LIGAMENT_ROM/LIGAMENT_ROM_V1.json :: decision_rows[4].{tension_limit_N,all_five_ligaments,using_ALL_alone,flexion_side_lost_by_using_one_ligament_deg}`.
- Ceiling **7**, angles **−90, +10**: `results/ASSEMBLY_LIGAMENT_ROM/LIGAMENT_ROM_V1.json :: decision_rows[5].{tension_limit_N,all_five_ligaments,using_ALL_alone}`; spacing **20**: `scope`; absent neutral sample: `tasks/assembly/ligament_rom_decision.py :: TABLE`.
- Stiffness **1734 ± 446**, facet fraction **0.03–0.25**: `results/ASSEMBLY_LOAD_PATH/LOAD_PATH_V1.json :: inputs_both_from_harvested_edges.{disc_axial_stiffness_N_per_mm,disc_axial_stiffness_sd,facet_compressive_fraction}`.
- Load **200**, interval **0.0688–0.1506**, widths **0.0254, 0.0546**: `results/ASSEMBLY_LOAD_PATH/LOAD_PATH_V1.json :: decision_rows[0].{load_N,displacement_low_mm,displacement_high_mm,width_from_facet_fraction_alone_mm,width_from_stiffness_alone_mm}`; conflicting recommendation: `which_uncertainty_governs`.
- Load **1000**, interval **0.344–0.7531**: `results/ASSEMBLY_LOAD_PATH/LOAD_PATH_V1.json :: decision_rows[2].{load_N,displacement_low_mm,displacement_high_mm}`.
- Muscles **187**, mean **≈4.0** (stored **4.015890405799041**), maximum **≈45.3** (stored **45.3419189654162**), fraction **≈0.702** (stored **0.7024520336728026**): `results/ASSEMBLY_LOAD_REDISTRIBUTION/decision.json :: {shared_muscles,surrogate_mae_N,surrogate_max_error_N,fraction_of_redistribution_explained}`.
- Rows **10**, poses **1**: `results/ASSEMBLY_LOAD_REDISTRIBUTION/decision.json :: {rows,effective_distinct_poses_per_arm}`; force-spread qualification: `largest_muscle_force_spread_across_rows_N`.
- Joint passes **0**, cases **92**: `results/LANE_SPINE_LIGAMENT_STIFFNESS/r5/CHRONOLOGICAL_RESPONSE_RESULTS_V1.json :: summary.chronological_linear.{joint_pass_count,cases}`; independently count `records[*].held[*].predictions.chronological_linear.{load_pass,work_pass}`.
- Cadaver load **1000**, areas **1150, 520**, peaks **3, 6**, factors **3.45, 3.12**: `source_documents/MECHANISM_MENISCUS_LOAD_DISTRIBUTION_evidence.json :: computed_evidence.raw_data.fukubayashi_kurosawa_1980.{load_N,intact,meniscectomy_total}` and `computed_evidence.computed_ratios.{fk_peak_over_mean_with,fk_peak_over_mean_without}`.
- Finite-element load **1000**, area/peak/capacity/ratio **110/1.2/132/0.132; 80/2.1/168/0.168; 35/3.3/115.5/0.1155; 95/1.4/133/0.133**: `results/LANE_MENISCUS_LOAD_LABEL/SOURCE_AUDIT_R1_V2.json :: FE_rows[0:4].{axial_load_N,reported_area_mm2,reported_peak_MPa,upper_force_at_central_values_N,local_load_fraction_upper_at_central_values}`.
- Partial-resection ratio **0.9**, predicted and control **≈1.111** (both stored **1.1111111111111112**), observed **1.65**: `results/LANE_MENISCUS_LOAD_LABEL/SOURCE_AUDIT_R1_V2.json :: baratz_transfer.{held_area_ratio,fixed_concentration_prediction,mean_pressure_ratio_control,held_peak_ratio}`.
- Sweep counts **8, 12**: `results/ASSEMBLY_LOAD_LABEL_SWEEP/SWEEP_V1.json :: {edges_with_a_product_below_its_stated_load,edges_with_pressure_and_force}`.
- Reference metadata and the percentage in reference [6]'s title: `results/LANE_DISC_LOAD_LABEL/r1/WILKE_1999_V1.txt :: title`; `results/LANE_DISC_LOAD_LABEL/r1/EUROPEPMC_10222525_V1.json :: resultList.result[id=10222525].{title,pmid,doi}`; `source_documents/MECHANISM_INTERVERTEBRAL_DISC.md :: Literature anchors, stiffness entry`; `source_documents/MECHANISM_SPINE_MOTION_SEGMENT_FACET.md :: Literature anchors, facet-load entry`; `results/LANE_MENISCUS_LOAD_LABEL/EUROPEPMC_PRIMARY_R1.json :: resultList.result[id=6894212|3755296|41539441].{title,pmid,doi,abstractText}`.
- Scope and withdrawals: `data/CONSTRAINT_NET_TISSUE.json :: bodytwin.tissue_constraint_net.edges[id=HARVEST-E0001|HARVEST-E0002|HARVEST-E0003|HARVEST-E0004|HARVEST-E0057|HARVEST-E0058|HARVEST-E0059|HARVEST-E0060|HARVEST-E0061|T-E23-meniscus_patch_share-vs-stated_axial_load|T-E31-nucleus_pressure_multiplier-vs-route_of_derivation|T-E37-ligament_slack-vs-work_sign]`; walking force mismatch: `results/LANE_DISC_LOAD_LABEL/r1/SOURCE_AUDIT_V1.json :: {activity_label,force_domain,same_subject,same_statistic}`; ligament admission: `results/LANE_SPINE_LIGAMENT_STIFFNESS/night_rounds/r5.json :: outcome.human_ligament_admission`.

Status: pending independent review. No claim of clinical validation.
