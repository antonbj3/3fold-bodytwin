# Final port table, including RESPONSE R4

| Stage | Port | Value / unit | Uncertainty | Status | Source / scope |
|---|---|---|---|---|---|
| tool→mechanics | hypodermic_16G.effective_J | 1153.4569983136594 J/m² | [1103.4569983136594, 1203.4569983136594] | MEASURED | LANE_SURGICAL_INCISION/PORTS_R4_FINAL.json; first-pass minus same-hole repeat effective work; intrinsic cleavage and terminal-process components unseparated; fullskin, tool-specific, figure reading not CI |
| tool→mechanics | hypodermic_18G.effective_J | 1801.0118043844857 J/m² | [1751.0118043844857, 1851.0118043844857] | MEASURED | LANE_SURGICAL_INCISION/PORTS_R4_FINAL.json; first-pass minus same-hole repeat effective work; intrinsic cleavage and terminal-process components unseparated; fullskin, tool-specific, figure reading not CI |
| tool→mechanics | hypodermic_21G.effective_J | 944.3507588532883 J/m² | [894.3507588532883, 994.3507588532883] | MEASURED | LANE_SURGICAL_INCISION/PORTS_R4_FINAL.json; first-pass minus same-hole repeat effective work; intrinsic cleavage and terminal-process components unseparated; fullskin, tool-specific, figure reading not CI |
| tool→mechanics | hypodermic_25G.effective_J | 2543.001686340641 J/m² | [2493.001686340641, 2593.001686340641] | MEASURED | LANE_SURGICAL_INCISION/PORTS_R4_FINAL.json; first-pass minus same-hole repeat effective work; intrinsic cleavage and terminal-process components unseparated; fullskin, tool-specific, figure reading not CI |
| incision→response | Gamma_cut_empirical_J_m2 | None J/m² | None | UNKNOWN | LANE_SURGICAL_INCISION/PORTS_R4_FINAL.json; Cannot fill from another assay |
| incision→response | Gamma0_intrinsic_J_m2 | None J/m² | None | UNKNOWN | LANE_SURGICAL_INCISION/PORTS_R4_FINAL.json; Cannot fill from another assay |
| incision→response | wound_gap_empirical_m | None m | None | UNKNOWN | LANE_SURGICAL_INCISION/PORTS_R4_FINAL.json; Cannot fill from another assay |
| incision→response | bleeding_flow_m3_s | None m³/s | None | UNKNOWN | LANE_SURGICAL_INCISION/PORTS_R4_FINAL.json; Cannot fill from another assay |
| incision→response | vessel_map | None documented map/law | None | UNKNOWN | LANE_SURGICAL_INCISION/PORTS_R4_FINAL.json; Cannot fill from another assay |
| incision→response | healing_strength_law | None documented map/law | None | UNKNOWN | LANE_SURGICAL_INCISION/PORTS_R4_FINAL.json; Cannot fill from another assay |
| injury zone | mechanical_strain_zone_halfwidth_m | None m | None | UNKNOWN | LANE_SURGICAL_INCISION/PORTS_R4_FINAL.json; Distinct widths; no mechanical→biological transfer |
| injury zone | cell_viability_zone_halfwidth_m | None m | None | UNKNOWN | LANE_SURGICAL_INCISION/PORTS_R4_FINAL.json; Distinct widths; no mechanical→biological transfer |
| injury zone | perfusion_loss_zone_halfwidth_m | None m | None | UNKNOWN | LANE_SURGICAL_INCISION/PORTS_R4_FINAL.json; Distinct widths; no mechanical→biological transfer |
| incision work | Gamma_cut_scenario | [150, 380] J/m² | scenario envelope, not CI | SYNTHETIC | LANE_SURGICAL_BINDINGS/PORTS_R4_FINAL.json; Explicit opt-in; not native cut |
| micro→macro | fibril_inventory_upper | 6965.365853658535 J/m² | borrowed interface bound and imposed V/A | SYNTHETIC | LANE_SURGICAL_BINDINGS/INTERFACE_BUDGET_RESULTS_R4.json; Conditional ceiling, not universal native upper bound |
| bridging | mode_work_reference | [30380, 20600] J/m² | SD 4900/2150; one juvenile porcine fullskin cohort | MEASURED | LANE_SURGICAL_BINDINGS/PORTS_R4_FINAL.json; Tear assay I/III, not cut or healing stress |
| bridging | native_terminal_traction | None Pa | None | UNKNOWN | LANE_SKIN_TOUGHNESS_GAP/PORTS_R2_FINAL_V2.json; Mode recruitment/anchor survival and cut connectivity absent |
| bridging | mm_bridge_work | [24500, 19137.5] J/m² | L5mm,eps.5 strength-as-traction; not CI | SYNTHETIC | LANE_SKIN_TOUGHNESS_GAP/MECHANISM_TABLE_R2_FINAL.json; Corrected crack/opening axes; anisotropy FAIL |
| bleeding | vessel_radius_histogram | None m and 1/m² per layer | None | UNKNOWN | BT-FW48-AUTO-55a5902a27c3f6/vessel_map_flow_results.json; Source histogram explicitly assumed; count control OVERpredicts25.3757x |
| hemostasis | shear_adhesion_law | None 1/min versus 1/s | None | UNKNOWN | BT-FW48-AUTO-27a6dac22a7a23/platelet_hemostasisis_shear.py; Qualitative anchors support synthetic rate closure; no native shear/closure assay |
| oxygen edge | dry_mass_density | None kg/m³ | None | UNKNOWN | LANE_SURGICAL_RESPONSE/RESULTS.md; Same-specimen quantities not identified by resting human skin proxies |
| oxygen edge | vascular_flow_Hb_map | None m³/s and mlO2/m³ | None | UNKNOWN | LANE_SURGICAL_RESPONSE/RESULTS.md; Same-specimen quantities not identified by resting human skin proxies |
| oxygen edge | oxygen_boundary | None Torr and mlO2/(m² min Torr) | None | UNKNOWN | LANE_SURGICAL_RESPONSE/RESULTS.md; Same-specimen quantities not identified by resting human skin proxies |
| oxygen edge | cell_damage_law | None 1/min | None | UNKNOWN | LANE_SURGICAL_RESPONSE/RESULTS.md; Same-specimen quantities not identified by resting human skin proxies |
| healing→strength | physical_collagen_birth_turnover | None reference mass/day | None | UNKNOWN | LANE_SURGICAL_RESPONSE/RESPONSE_PORTS_R3.json; Tracer≠mass, HP≠total maturity; species/cohort joint posterior null |
| healing→strength | catalytic_LOX | None 1/day | None | UNKNOWN | LANE_SURGICAL_RESPONSE/RESPONSE_PORTS_R3.json; Tracer≠mass, HP≠total maturity; species/cohort joint posterior null |
| healing→strength | joint_Q_U_Q_I_Q_M | None reference mass tensor | None | UNKNOWN | LANE_SURGICAL_RESPONSE/RESPONSE_PORTS_R3.json; Tracer≠mass, HP≠total maturity; species/cohort joint posterior null |
| healing→strength | chemistry_to_strength_law | None Pa or %intact | None | UNKNOWN | LANE_SURGICAL_RESPONSE/RESPONSE_PORTS_R3.json; Tracer≠mass, HP≠total maturity; species/cohort joint posterior null |
| healing→strength | chemistry_61_90 | None mol/molcollagen | None | UNKNOWN | LANE_SURGICAL_RESPONSE/RESPONSE_PORTS_R3.json; Tracer≠mass, HP≠total maturity; species/cohort joint posterior null |
| healing→strength | R3_HP_input_strength_RMSE | 8.306367094671637 percentage points | conditional observation proxy, nominal gate FAIL | SYNTHETIC | LANE_SURGICAL_RESPONSE/r3/crosslink_port/summary.json; HP42 input, not predicted from early chemistry; Levenson heldout diagnostic |
| final chemistry R4 | native precursor D10 | None mol/molreference collagen | None | UNKNOWN | LANE_SURGICAL_RESPONSE/RESPONSE_PORTS_R4.json; Same reference and absolute glyco pools missing |
| final chemistry R4 | cell-free first-order control k,eta | [0.06949679103768976, 0.2833333333333334] 1/day,1 | Apparent endpoint-derived analogue; first/second order nonidentified | MEASURED | LANE_SURGICAL_RESPONSE/RESPONSE_PORTS_R4.json; External cell-culture typeI assay, native transfer synthetic |
| final chemistry R4 | HP inventory42/10 | 2.220261437908497 1 | [1.747159090909091, 2.8467987804878048] | MEASURED | LANE_SURGICAL_RESPONSE/r4/SCORES_R4_v2.json; Ratio of group means; fixed assay reference; graphical box, not CI |
