# Running

Status values: VERIFIED-FRESH = the module's own gates passed in the recorded run; OWN-GATE-FAIL = one or more own gates failed, result retained and unpromoted; CUDA-ONLY = GPU execution not yet verified in the recorded scope; SYNTHETIC-ONLY = passing evidence limited to synthetic inputs. Earlier attempts are notes, not additional current module rows.

    python -m venv .venv && .venv/bin/pip install -r requirements.txt          # numpy, scipy, sympy, pytest
    BODYTWIN_OUT=./outputs .venv/bin/python src/bodytwin/cells/nervous/na_k_atpase.py   # any cell: prints its gates, writes $BODYTWIN_OUT/<cell>/<cell>_results.json
    cp -r examples/synthetic_inputs/* ./outputs/                                    # seeds the six producers that are not in this repository
    BODYTWIN_OUT=./outputs .venv/bin/python src/bodytwin/cells/cardiovascular/cardiac_output.py   # then any consumer, producers first (tests/cell_dependencies.json)
    BODYTWIN_ROOT=./examples/mini_graph BODYTWIN_DESIGNS=./examples/mini_graph/designs/cell_designs.json .venv/bin/python src/bodytwin/framework/bodytwin_build_anchor_graph.py
    BODYTWIN_ROOT=./examples/mini_graph .venv/bin/python src/bodytwin/framework/bodytwin_hardening_check.py   # exit 1 on any schema or fold-path violation
    .venv/bin/python -m pytest -q tests    # every cell once, producers before consumers, in one temporary BODYTWIN_OUT; the fold layer; the compiler example

A cell that reads another cell's result looks for it under `$BODYTWIN_OUT/<producer>/`; `tests/cell_dependencies.json` is the edge list. Cells marked SYNTHETIC-ONLY below are rooted in one of the synthetic producers under `examples/synthetic_inputs/`; their numbers are not the source project's numbers.

## Status

| Module | Status | Evidence |
|---|---|---|
| src/bodytwin/geometry/optics/tissue_photon_segmented_v1.py | VERIFIED-FRESH | 5/5full transport gates on L4(128photons) and H100(1024); all5arrays exact twice vs original, API28.4s to1.44-1.78s and30.6s to0.888-0.893s;10boundary controls exact. |
| src/bodytwin/geometry/optics/tissue_photon_prepared_v1.py | VERIFIED-FRESH | 5/5 on L4/H100; full five arrays exact versus segmented baseline and repeats. Reuse1434-1445->410ms L4,942-958->230ms H100; cold H1001038ms retained. docs/photon_prepared.md. |
| src/bodytwin/framework/coupling_evidence_v1.py | SYNTHETIC-ONLY | 5/5 controls; explicit PASS/FAIL/ABSTAIN,10malformed input/declaration refusals; absent evidence never accepted. |
| `src/bodytwin/framework/bodytwin_fold.py` | VERIFIED-FRESH |
| `src/bodytwin/framework/bodytwin_build_anchor_graph.py` | VERIFIED-FRESH |
| `src/bodytwin/framework/bodytwin_hardening_check.py` | VERIFIED-FRESH |
| `src/bodytwin/framework/bodytwin_graph_io.py` | VERIFIED-FRESH |
| `src/bodytwin/framework/bodytwin_fold_skip_rules.py` | VERIFIED-FRESH |
| `src/bodytwin/framework/bodytwin_scorecard.py` | VERIFIED-FRESH |
| `src/bodytwin/framework/bodytwin_provenance_scan.py` | VERIFIED-FRESH |
| `src/bodytwin/framework/bodytwin_writeback_provenance_scan.py` | VERIFIED-FRESH |
| `src/bodytwin/framework/mt_id_alias.py` | VERIFIED-FRESH |
| `src/bodytwin/framework/mt_extract.py` | VERIFIED-FRESH |
| `examples/cell_graph/build_cell_graph.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/aging/cellular_senescence.py` | VERIFIED-FRESH | Current wording-bound replay 4/4 gates; two complete receipts exact, existing numerical gate values unchanged. reports/wording_refresh/aging_cells.json; archived scope in docs/wording_verification_refresh.md. |
| `src/bodytwin/cells/aging/telomere_attrition.py` | VERIFIED-FRESH | Current wording-bound replay 4/4 gates; two complete receipts exact, existing numerical gate values unchanged. reports/wording_refresh/aging_cells.json; archived scope in docs/wording_verification_refresh.md. |
| `src/bodytwin/cells/cardiovascular/arterial_nonlinear_compliance_pressure_sweep_deferred_arithmetic.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/cardiovascular/arterial_pressure.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/cardiovascular/athlete_heart_remodeling.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/cardiovascular/baroreflex.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/cardiovascular/baroreflex_sigmoid_symmetry_vs_asymmetric_anchor_deferred_arithmetic.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/cardiovascular/bp_map_bikia_crosscheck.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/cardiovascular/capillary_starling.py` | VERIFIED-FRESH | Current wording-bound replay 4/4 gates; two complete receipts exact, existing numerical gate values unchanged. reports/wording_refresh/explicit_all_gates.json; archived scope in docs/wording_verification_refresh.md. |
| `src/bodytwin/cells/cardiovascular/cardiac_cicr_ecc.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/cardiovascular/cardiac_cicr_ode_model.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/cardiovascular/cardiac_detraining_reversibility_tau.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/cardiovascular/cardiac_output.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/cardiovascular/cardiac_output_exercise_heat_oversubscription_deferred_arithmetic.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/cardiovascular/cardiac_output_geometric.py` | VERIFIED-FRESH | Current wording-bound replay 4/4 gates; two complete receipts exact, existing numerical gate values unchanged. reports/wording_refresh/pipeline_cells.json; archived scope in docs/wording_verification_refresh.md. |
| `src/bodytwin/cells/cardiovascular/cardiac_pacemaker_funny_current.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/cardiovascular/cardiac_output_resting_spread_gate.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/cardiovascular/coronary_absolute_flow.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/cardiovascular/coronary_blood_flow.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/cardiovascular/exercise_bloodflow_redistribution.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/cardiovascular/murray_law_vascular_branching.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/cardiovascular/raas.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/cardiovascular/venous_return.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/cardiovascular/ventricular_ap_core_tt04style.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/cardiovascular/wall_tension_myogenic_generator.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/cardiovascular/windkessel_ode_rebuild.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/endocrine/adrenal_four_syndrome_state_machine.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/endocrine/betacell_bursting_chay_keizer.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/endocrine/betacell_reserve_auc_gap_and_direct_or_crosscheck.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/endocrine/calcium_pth_vitd.py` | VERIFIED-FRESH | Current wording-bound replay 4/4 gates; two complete receipts exact, existing numerical gate values unchanged. reports/wording_refresh/coupled_activation.json; archived scope in docs/wording_verification_refresh.md. |
| `src/bodytwin/cells/endocrine/calcium_pth_vitd_peterson_shrestha_crosscheck.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/endocrine/circadian_rhythm.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/endocrine/endo_receptor_reserve_elasticity_rebuild.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/endocrine/gh_igf1_pulsatile_chareq_dead_cell_rebuild_2026_07_28.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/endocrine/glucagon_counterregulation.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/endocrine/glucose_insulin_minimal_model.py` | VERIFIED-FRESH | Current wording-bound replay 4/4 gates; two complete receipts exact, existing numerical gate values unchanged. reports/wording_refresh/pipeline_cells.json; archived scope in docs/wording_verification_refresh.md. |
| `src/bodytwin/cells/endocrine/glucose_meal_dallaman2007.py` | OWN-GATE-FAIL | Fresh typed replay12/19 required gates,7failures despite exit0; complete3069-byte records repeat. |
| `src/bodytwin/cells/endocrine/glucose_meal_egp_gate_reanchored.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/endocrine/hpa_cortisol_axis.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/endocrine/hpg_male_axis.py` | VERIFIED-FRESH | Current wording-bound replay 4/4 gates; two complete receipts exact, existing numerical gate values unchanged. reports/wording_refresh/strict_endocrine.json; archived scope in docs/wording_verification_refresh.md. |
| `src/bodytwin/cells/endocrine/hpg_testosterone_dose_response_bhasin2001.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/endocrine/insulin_portal_peripheral_egp_trajectory_deferred_arithmetic.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/endocrine/insulin_ultradian_sturis_hopf.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/endocrine/lactation_letdown.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/endocrine/ovarian_cycle_axis.py` | VERIFIED-FRESH | Current wording-bound replay 4/4 gates; two complete receipts exact, existing numerical gate values unchanged. reports/wording_refresh/strict_endocrine.json; archived scope in docs/wording_verification_refresh.md. |
| `src/bodytwin/cells/endocrine/parturition_myometrium.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/endocrine/sperm_motility_axoneme.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/endocrine/spermatogenesis_sertoli.py` | OWN-GATE-FAIL | Fresh21/21 numerical gates and finite flag pass; complete19928byte records differ at one timestamp leaf, so full-record repeat fails. |
| `src/bodytwin/cells/endocrine/thyroid_metabolic_axis.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/endocrine/vitamin_d_activation.py` | OWN-GATE-FAIL | Absent input28/29; empty and partial synthetic sibling blocks incorrectly emit30/30 with0/3 and1/3comparisons, full repeats exact. |
| `src/bodytwin/cells/energy/cori_hepatic_atp_affordability.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/energy/cori_lactate_hgp_crosscheck.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/energy/enzyme_eyring_ceiling.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/energy/exercise_heat_shared_budget.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/energy/glycolytic_oscillations_selkov.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/energy/krogh_bmr_voidfloor_gate.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/energy/mitochondrial_oxphos.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/energy/muscle_energetics.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/energy/organ_co2_fick_rq_partition.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/energy/organ_o2_consumption_partition.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/energy/organ_o2_gate_power_montecarlo.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/energy/protein_turnover_atp_energy_ledger.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/energy/protein_turnover_recycled_fraction_deferred_arithmetic.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/energy/warburg_metabolism.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/energy/whole_body_lactate_ledger.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/eye/aqueous_humor_iop.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/eye/corneal_transparency.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/eye/lens_accommodation.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/eye/tear_film_model.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/gastrointestinal/bile_acid_postprandial_zscore_crosscheck.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/gastrointestinal/bile_enterohepatic.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/gastrointestinal/deglutition_swallowing.py` | VERIFIED-FRESH | Current wording-bound replay 5/5 gates; two complete receipts exact, existing numerical gate values unchanged. reports/wording_refresh/embedded_overall.json; archived scope in docs/wording_verification_refresh.md. |
| `src/bodytwin/cells/gastrointestinal/emesis_reflex.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/gastrointestinal/gastric_acid_secretion.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/gastrointestinal/gi_motility_slow_waves.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/gastrointestinal/gut_microbiome.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/gastrointestinal/hepatic_clearance.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/gastrointestinal/hepatic_splanchnic_series_topology_check.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/gastrointestinal/hepatocyte_count_clearance_factorization.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/gastrointestinal/pancreatic_exocrine.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/gastrointestinal/urea_cycle_flux_jin_recalibration.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/growth/actin_treadmilling_pollard_kinetics.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/growth/angiogenesis_vegf.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/growth/cytoskeleton_critical_concentration.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/growth/cellcycle_restriction_switch_biomd265_rebuild.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/growth/cytoskeleton_dynamics.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/growth/dna_repair_kinetics.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/growth/dna_replication_fork.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/growth/dna_replication_kinetic_proofreading_drake.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/growth/dna_supercoiling_topology_rebuild.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/growth/fetal_circulation.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/growth/genetic_code_error_minimization.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/growth/gut_crypt_littles_law_renewal.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/growth/hedgehog_gli_8state_ode.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/growth/liver_regeneration.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/growth/mrna_translation_kinetics.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/growth/neonatal_transition.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/growth/notch_delta_lateral_inhibition_rebuild.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/growth/notch_lateral_inhibition_dead_cell_rebuild_2026_07_28.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/growth/placental_transfer.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/growth/somite_hes7_lewis_zeiser_dde_rebuild.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/growth/spindle_sac_biomd186_rebuild.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/growth/tumor_adaptive_therapy_lotka_volterra.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/haematology/b12_serum_functional_dissociation.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/haematology/bilirubin_metabolism.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/haematology/blood_oxygen_transport.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/haematology/epo_hematocrit_regime_sweep_deferred_arithmetic.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/haematology/erythropoiesis.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/haematology/hema_f8_intron22_frequency_crosscheck.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/haematology/hematopoiesis.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/haematology/hemoglobin_hill_regime_crosscheck.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/haematology/hemoglobin_mwc_allostery_rebuild.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/haematology/iron_hepcidin.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/haematology/iron_hepcidin_conservation_blind_check.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/haematology/iron_hepcidin_dynamic_4state.py` | VERIFIED-FRESH | Current wording-bound replay 4/4 gates; two complete receipts exact, existing numerical gate values unchanged. reports/wording_refresh/explicit_all_gates.json; archived scope in docs/wording_verification_refresh.md. |
| `src/bodytwin/cells/haematology/iron_oxygen_carriage_conservation.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/haematology/platelet_hemostasis.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/haematology/rbc_osmotic_volume_crosscheck.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/immune/acute_phase_inflammation.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/immune/bcell_affinity_maturation.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/immune/complement_ap_bistability_eculizumab.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/immune/complement_cascade.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/immune/ferroptosis_cert.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/immune/macrophage_polarization.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/immune/necroptosis_cert.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/immune/tcell_activation_exhaustion.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/musculoskeletal/achilles_elastic_energy_return_running.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/musculoskeletal/bmu_formation_period_definition_check.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/musculoskeletal/bmu_turnover_kinetics.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/musculoskeletal/bone_calcium_flux_reconciliation.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/musculoskeletal/bone_composite_modulus_bounds.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/musculoskeletal/bone_fracture_toughness_lefm.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/musculoskeletal/bone_rankl_opg_lemaire2004.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/musculoskeletal/bone_remodeling.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/musculoskeletal/bone_remodeling_transient_arithmetic.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/musculoskeletal/bone_strain_bounded_classifier_proposal.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/musculoskeletal/bone_wolff_law_mechanostat.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/musculoskeletal/cartilage_poroelastic_relaxation.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/musculoskeletal/collagen_triple_helix_thermal_stability.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/musculoskeletal/force_velocity_power_rfd.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/musculoskeletal/meniscus_hoop_stress_model.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/musculoskeletal/motor_unit_recruitment.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/musculoskeletal/muscle_memory_substrate_function_dissociation.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/musculoskeletal/muscle_pcsa_crossbridge_specific_tension_rebuild.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/musculoskeletal/reproduce_motor_unit_force_factorization.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/musculoskeletal/sarcomere_force_rice2008.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/musculoskeletal/sprint_spring_mass.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/musculoskeletal/sprint_topspeed_force_limit_rebuild.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/musculoskeletal/tendon_collagen_hierarchical_mechanics.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/musculoskeletal/tendon_incidence_cagr_crosscheck.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/musculoskeletal/tmj_lever_model.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/musculoskeletal/trabecular_stress_orthogonality.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/nervous/alphasyn_prion_aggregation.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/nervous/amygdala_fear_circuit.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/nervous/amyloid_beta_aggregation.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/nervous/basal_ganglia_gating.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/nervous/blood_brain_barrier.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/nervous/cerebral_autoregulation.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/nervous/corticospinal_motor_command.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/nervous/csf_davson_icp_flux_deferred_arithmetic.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/nervous/dopamine_kinetics.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/nervous/fingertip_tactile.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/nervous/glymphatic_peclet_bound.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/nervous/hairbundle_corner_frequency_gating_spring.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/nervous/hh_spike_overlap_reshape_sim.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/nervous/hodgkin_huxley_squid_axon_ap_rebuild.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/nervous/intracranial_pressure.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/nervous/myelin_gratio_rebuild.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/nervous/na_k_atpase.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/nervous/neurovascular_coupling.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/nervous/nmj_synaptic_delay_adjudication.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/nervous/olfactory_transduction.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/nervous/reward_prediction_error.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/nervous/rhodopsin_thermal_dark_noise_arrhenius_gate.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/nervous/rod_darknoise_poisson_pooling_threshold.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/nervous/serotonin_system.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/nervous/sleep_two_process_recovery_rebound.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/nervous/spinal_cpg.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/nervous/synaptic_tsodyks_markram_ppr.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/nervous/synaptic_vesicle_release.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/nervous/taste_transduction.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/nervous/tau_pathology.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/nervous/vestibular_vor.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/nervous/wilson_cowan_ei_seizure_bifurcation.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/organ_systems/countercurrent_multiplier.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/organ_systems/hair_follicle.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/organ_systems/melanin_photoprotection.py` | VERIFIED-FRESH | Current wording-bound replay 4/4 gates; two complete receipts exact, existing numerical gate values unchanged. reports/wording_refresh/additional_cells.json; archived scope in docs/wording_verification_refresh.md. |
| `src/bodytwin/cells/organ_systems/nephron_transport_cell.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/organ_systems/skin_pulp_mechanics.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/organ_systems/thermoregulation.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/renal/acid_base_co2.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/renal/countercurrent_crossspecies_adjudication.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/renal/donnan_electroneutrality_osmolalgap.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/renal/mineral_calcium_gross_net_po4_ledger_deferred_arithmetic.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/renal/nephron_tgf_oscillation_rebuild.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/renal/nitrogen_urea_cycle_unit_rda_nh4_deferred_arithmetic.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/renal/renal_adh_osmolality_robertson_threshold.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/renal/renal_ammoniagenesis_nae_grid_deferred_arithmetic.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/renal/renal_ammoniagenesis_nitrogen_coupling.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/renal/renal_c3g_joint_lesion_probability_deferred_arithmetic.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/renal/renal_countercurrent_allometric_exponent_crosscheck.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/renal/renal_egfr_hr_compounding_crosscheck.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/renal/renal_filtration.py` | VERIFIED-FRESH | Current wording-bound replay 4/4 gates; two complete receipts exact, existing numerical gate values unchanged. reports/wording_refresh/coupled_activation.json; archived scope in docs/wording_verification_refresh.md. |
| `src/bodytwin/cells/renal/renal_medullary_washout_gradient_tradeoff.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/renal/renal_nephron_segment_fraction_grid_deferred_arithmetic.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/renal/sodium_reabsorption_fe_na_chain_deferred_arithmetic.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/renal/water_osmotic_ledger_leaveoneout_deferred_arithmetic.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/respiratory/altitude_acclimatization.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/respiratory/carbon_co2_conservation_closure.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/respiratory/co2_chemoreflex_rebuild.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/respiratory/cough_reflex.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/respiratory/laplace_alveolar_surfactant_stability.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/respiratory/mucociliary_clearance.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/respiratory/pulmonary_gas_exchange.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/respiratory/pulmonary_surfactant_alveolar_stability.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/respiratory/reproduce_alveolar_surface_factorization.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/respiratory/respiratory.py` | SYNTHETIC-ONLY |
| `src/bodytwin/cells/signalling/apoptosis_intrinsic_extrinsic.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/signalling/autophagy_mitophagy.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/signalling/calcium_ip3_lirinzel_oscillation.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/signalling/circadian_q10_compensation_goodwin.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/signalling/genetic_toggle_switch_gardner_collins.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/signalling/insulin_pi3k_akt_signaling.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/signalling/mapk_erk_cascade.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/signalling/mtorc1_feedback_rebuild.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/signalling/nfkb_hoffmann2002_biomd140_rebuild.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/signalling/nfkb_signaling_dynamics.py` | VERIFIED-FRESH | Current wording-bound replay 4/4 gates; two complete receipts exact, existing numerical gate values unchanged. reports/wording_refresh/additional_cells.json; archived scope in docs/wording_verification_refresh.md. |
| `src/bodytwin/cells/signalling/smad_nucleocytoplasmic_shuttling_biomd173.py` | VERIFIED-FRESH |
| `src/bodytwin/cells/signalling/wnt_betacatenin_rebuild.py` | VERIFIED-FRESH |
| `src/bodytwin/chains/metabolic_chain_v1.py` | VERIFIED-FRESH |
| `examples/anatomy/inspect_mesh.py` | OWN-GATE-FAIL | Two raw atlas loads identical; original1497 vertices/2982 faces fail closed gate with10 unpaired edges, no degenerate faces. |
| `src/bodytwin/geometry/mesh_ingest_v1.py` | VERIFIED-FRESH | Atlas exact welding preserves all triangle coordinates; tetrahedron boundary volume1/3, six faces, eight invalid-input rejections; explicit SI array conversion exact. |
| `src/bodytwin/geometry/__init__.py` | VERIFIED-FRESH | Imported by the complete CPU compose run. |
| `examples/anatomy/compose_mesh_to_field.py` | VERIFIED-FRESH | CPU six gates PASS; voxel volume103640 vs mesh103741.15500169418mm3, relative error0.0009750710958686, two full array/results identical. |
| `src/bodytwin/geometry/optics/tissue_photon_mc_v3.py` | OWN-GATE-FAIL | Relocated cube full integer hashes unchanged; wider boundary failures remain 1467 leaks and slab 6 leaks/860 caps. |
| `src/bodytwin/geometry/optics/tissue_photon_mc_v1.py` | OWN-GATE-FAIL | Byte-identical relocation of the frozen failed comparator; source SHA in optics_package_move.json. Earlier 578026 failed photons remain uncorrected. |
| `src/bodytwin/geometry/optics/tissue_photon_mc_v2.py` | OWN-GATE-FAIL | Byte-identical relocation; earlier one failed photon remains uncorrected. |
| `src/bodytwin/geometry/optics/__init__.py` | VERIFIED-FRESH | Package imports in all six completed L4 runs; all six moved file hashes unchanged. |
| `src/bodytwin/geometry/optics/tissue_photon_mc_v4.py` | OWN-GATE-FAIL | Original cube6/6 PASS; refractive cube4/6 with1 leak; slab6/8 with6 leaks and0 caps; all full repeats exact. |
| `src/bodytwin/geometry/optics/tissue_photon_mc_v5.py` | OWN-GATE-FAIL | Two full L4 runs per case;999993/999994/986347 failed packets. Failed variant frozen. |
| `src/bodytwin/geometry/optics/tissue_photon_mc_v6.py` | VERIFIED-FRESH | Cube6/6, refractive cube6/6, slab8/8 PASS; zero residual/leaked/capped packets in all six complete runs; full arrays exact twice. |
| `src/bodytwin/geometry/layered_box_v1.py` | SYNTHETIC-ONLY | Four/four CPU gates PASS, all array hashes exact twice;12 vertices/22 faces, per-region volumes36000/180000. |
| `src/bodytwin/geometry/optics/tissue_photon_paths_v1.py` | VERIFIED-FRESH | All7 observation/detector gates PASS, canonical full transport unchanged, complete arrays exact twice. |
| `src/bodytwin/geometry/optics/nirs_forward_v1.py` | SYNTHETIC-ONLY | Five engineering gates PASS twice; nirs_forward_v1.json; clinical slope certification remains false. |
| `src/bodytwin/chains/nirs_chain_v1.py` | SYNTHETIC-ONLY | All8 conditional chain gates PASS twice; joint span0.055127841665889354<=0.06987483101741154; nirs_chain_v1.json. |
| `src/bodytwin/geometry/layered_cylinder_v1.py` | SYNTHETIC-ONLY | Curved fixture geometry gates PASS; region volume error0.001605606964382, closed oriented regions and exact repeats. |
| `src/bodytwin/chains/nirs_curved_chain_v1.py` | SYNTHETIC-ONLY | Eight/eight conditional gates PASS twice; joint sample span0.05616480545855<=0.07217381811722 isolated sum; exact hashes/nulls. |
| src/bodytwin/geometry/label_interfaces_v1.py | OWN-GATE-FAIL | 4/5 gates; edge_contact has1 nonmanifold edge with4 incidents instead of2; all volumes/orientations/full arrays exact. |
| src/bodytwin/geometry/tetra_interfaces_v1.py | SYNTHETIC-ONLY | 6/6 gates on four valid/two invalid controls; reports/tetra_interface_probe.json; cell inversion gives identical output. |
| src/bodytwin/geometry/optics/ocular_photon_v1.py | SYNTHETIC-ONLY | 6/6 gates at two wavelengths x two1e6 runs; exact packet energy, zero leaks/caps/residuals; reports/ocular_transport_1000000.json. |
| src/bodytwin/geometry/optics/ocular_receiver_v1.py | SYNTHETIC-ONLY | 3/3 normalization gates; zero integral residual, full maps repeat, four invalid controls rejected; reports/ocular_receiver.json. |
| src/bodytwin/chains/nirs_spec_v1.py | SYNTHETIC-ONLY | 6/6 spec gates, both frozen chains byte-identical; typed seams, deterministic source and six independent generated child observations exact; reports/nirs_spec_probe.json. |
| src/bodytwin/geometry/tetra_ray_walk_v1.py | SYNTHETIC-ONLY | 5/5 straight-ray gates,512 rays with exact cell sequences and maximum4.440892098500626e-16 parameter error; reports/tetra_ray_walk.json. |
| src/bodytwin/geometry/tetra_ray_batch_v1.py | SYNTHETIC-ONLY | 6/6 gates;512 rays,2735 full records bit-identical to frozen scalar and repeats;4.221-4.277x paired CPU gain,32/32 controls. docs/tetra_ray_batch.md. |
| reports/tetra_partition_external_summary.json | OWN-GATE-FAIL | Pinned public benchmark:24 cell components,170 outer edges with4 incidents and1 with6,812 regional edge defects; full repeated counts exact, volume error1.4877210207958217e-16. |
| src/bodytwin/geometry/tetra_ray_cuda_v1.py | SYNTHETIC-ONLY | 4/4 CUDA straight-ray gates; all512 cell sequences and interval values exactly equal frozen CPU results, full repeats exact; reports/tetra_ray_cuda.json. |
| src/bodytwin/geometry/optics/tetra_attenuation_v1.py | SYNTHETIC-ONLY | 5/5 scoring gates; full int64 cell-dose/terminal arrays exactly equal sequential CPU reference and repeat; reports/tetra_attenuation.json. |
| src/bodytwin/geometry/optics/pulse_path_inverse_v1.py | SYNTHETIC-ONLY | 4/4 homogeneous/ownership/refusal selftest gates twice; supplied optical-model inversion only, no clinical calibration. |
| src/bodytwin/geometry/optics/profiled_information_v1.py | VERIFIED-FRESH | 7/7 numerical selftest controls; reusable SVD projection handles redundant nuisance columns, empty nuisance and underdetermined interests. |
| src/bodytwin/geometry/optics/grouped_detector_selection_v1.py | VERIFIED-FRESH | 5/5 selftest gates,11 invalid/unresolved experiments refused; canonical tie order and no input mutation. |
| examples/verification/verify_bodytwin.py | OWN-GATE-FAIL | Thirteen explicit recipes bind27/237fresh declarations; embedded-overall full reports match Modal/local and repeat exactly. Full coverage fails210unmapped; one timestamp-dependent source demoted separately. Earlier isolated host1/4 exact report matches retained; no universal portability claim. |
| Makefile | OWN-GATE-FAIL | verify-detectors passes; full verify exits2 with0numerical children because explicit coverage is incomplete. |
| src/bodytwin/framework/numerical_gate_contract_v1.py | VERIFIED-FRESH | 6/6 explicit gate-contract controls,7record/6declaration refusals; archived16record parity4/4, including12/19failure, exact repeats. |
| src/bodytwin/geometry/optics/profiled_information_bound_v1.py | VERIFIED-FRESH | 5/5 scalar perturbation controls,64corners/6refusals; archived integration4/4 gives evaluated gain lower1.4728425209401963, full repeats exact; not interval arithmetic. |



## Coupled chains

`src/bodytwin/chains/metabolic_chain_v1.py` is the first chain whose certificate is composed across
cells instead of asserted per cell. Four shipped cells are wired along one physical quantity each,
using only functions those cells already export:

    BODYTWIN_OUT=./outputs .venv/bin/python src/bodytwin/chains/metabolic_chain_v1.py   # 30 s, CPU, exit 1 on any FAIL
    .venv/bin/python -m pytest -q tests/test_metabolic_chain.py

`glucose_insulin_minimal_model.ogtt_leg` (meal load -> plasma glucose -> suprabasal glucose disposal
by mass balance on its own state) -> `mitochondrial_oxphos.variant_grid` P/O with
`warburg_metabolism.Y_OXPHOS_MODERN` (disposal -> ATP flux and the O2 its oxidation costs) ->
`blood_oxygen_transport.sao2_hill`/`cao2_ml_dl` (blood gases -> a-vO2 content difference) ->
`cardiac_output.fick_q_l_min`/`hr_from_q_sv` (required cardiac output and heart rate). Input is the
declared synthetic reference body (mass from `examples/synthetic_inputs/metabolic_cost`, every other
value the owning cell's own published constant); no individual's measurements enter anywhere.

Nominal, 75 g glucose load: peak glucose 184.5 mg/dL, 2 h glucose 7.50 mmol/L, suprabasal disposal
367.7 mg/min = 2.041 mmol/min; ATP flux 61.59 mmol/min at 30.18 ATP/glucose; O2 required 257.8
mL/min on top of the 250 mL/min resting anchor; SaO2 96.9 %, SvO2 75.1 %, a-vO2 difference 4.55
mL/dL; required cardiac output 11.154 L/min (baseline 5.491, +5.663) at 143.2 bpm. The chain
oxidises the whole disposed load inside the window, so that cardiac output is an upper bound (no
shipped cell splits oxidation from glycogen storage, and no such split is invented here).

Uncertainty propagation, seed 20260912, 256 draws, each stage's inputs drawn inside the band its own
cell states for them: final cardiac output 10.537 L/min mean, range [8.034, 14.136] L/min,
half-width 3.051, SD 1.276 L/min. Per-stage half-widths: S1 glucose/insulin 0.3403, S2 oxphos/ATP
0.0000, S3 blood O2 2.3235, S4 cardiac output 1.0969 L/min. **The blood-oxygen stage dominates** --
haemoglobin concentration over the transport cell's own 14-18 g/dL reference band moves the required
cardiac output more than every other declared tolerance in the chain combined. S2 contributes
exactly zero because the P/O variant changes the ATP return on a glucose, not the O2 that oxidising
it costs; the P/O spread shows up in the reported ATP flux only.

Cert composition: full-chain half-width 3.0509 L/min <= linear sum of the stage half-widths 3.7607
L/min, so composition is sub-additive and no stage is flagged nonlinear. Null cases: a zero glucose
load leaves glucose at 90 mg/dL, disposal, ATP flux and O2 demand at exactly 0 and the cardiac
output at its baseline; a doubled (150 g) load raises glucose peak, disposal, ATP flux, O2 demand
and cardiac output, and leaves the a-vO2 difference unchanged (it does not depend on the load).

## Long-running cells

`fingertip_tactile`, `cellcycle_restriction_switch_biomd265_rebuild` and `gi_motility_slow_waves`
run 86 s, 160 s and 303 s here (CPU, `nice -n 15`). None of them has a size or iteration flag, so
`tests/test_cells.py` gives exactly those three a 400 s cap (`TIMEOUTS`) instead of reducing their
work; the default cap for every other cell stays 120 s. All three passed here at full size, which is
what their VERIFIED-FRESH rows above mean.

## Scrub gate

A case-insensitive grep over the repository for absolute home paths, internal project names, personal names, device threads and non-redistributable asset names, plus a secondary sweep for individual-trial labels and scratch/agent paths, leaves these justified classes:

| hit | justification |
|---|---|
| `operat*` | the ordinary words operating point, operates, operation, cooperativity, and the `"operator"` field of a gate comparison (`"<="`, `">"`) |
| `proton` | the particle (H+) |
| `coil` | coiled-coil domain, DNA supercoiling, elastic recoil |
| author surnames (Antonarakis, Antony, Davidson) | literature citations |
| the source module prefix | the source file map |
| `walking1` in `cardiac_output_resting_spread_gate.py` | JSON key `reference_walking1_fick_chain`: keys are frozen, so only the person-denoting prefix of the source key was renamed; `walking1` is the source project's own trial label and carries no personal measurement -- the anchor it names is read from the `cardiac_output` cell result |

No home path, individual-trial subject label, scratch path or personal name remains under
`src/`, `examples/` or `tests/`.

## Public anatomy mesh seam: preregistered mechanism measurement

Inspect the unmodified atlas triangle array before adapter design: vertex/face
counts, finite coordinates, degenerate triangles, edge incidence, consistent
winding, bounds, surface area and signed volume. Explicit fixture coordinates
are declared in mm; no units are inferred from an OBJ header or a body size.
Two independent loads must produce identical diagnostics and array hashes.
Fixed gates: every undirected edge has exactly two incident triangles, consistent
winding, no zero-area face, positive signed volume and finite coordinates.
A failed original mesh is reported before any optional separate repaired fixture.
No automatic repair, rescaling, decimation, or orientation reversal in this probe.


Mechanism table before ingest design:

| representation | vertices | faces | unpaired edges | winding consistent | signed volume mm3 |
|---|---|---|---|---|---|
| original atlas indices | 1497 | 2982 | 10 | true | 103741.15500169418 |
| exact coordinate deduplication, diagnostic | 1493 | 2982 | 0 | true | 103741.15500169418 |

The four duplicate coordinates explain this topology mismatch: `unique(axis=0)`
and index remapping retain every face corner exactly. No hole filling or surface
motion is needed. Original closed gate FAIL is retained in
`examples/anatomy/mesh_mechanism.json`; the original OBJ is untouched.

Ingest design: strict validation by default, optional explicitly requested exact
coordinate welding with reported count; triangle arrays in mm for field, explicit
0.001 scale for motion arrays. Tetrahedra must have positive signed volume, no
repeated cell or nonmanifold face; export oriented boundary triangles. No implicit
unit guessing, shape repair, or surface resampling.

Compose gates, fixed before run: original mesh rejected without explicit welding;
welded face-corner arrays exactly unchanged; tet boundary closed/oriented and
volume agrees with tet sum to1e-12 relative; malformed indices/nonfinite input,
unknown units, degenerate/inverted/duplicate tetrahedra rejected; two complete
seam/field arrays bit-identical; frozen field gate OK; sparse/dense occupancy
mismatch zero; atlas voxel volume error<=5% (existing compose band) AND within a
conservative raster bound. The raster bound counts the union of voxel cubes
intersecting triangle AABBs, expanded by half a voxel plus the frozen ray jitter;
its volume bounds potentially misclassified boundary cubes. Fixed pitch2mm,
margin0; do not tune the pitch or either acceptance bound after the run.

Initial compose submission exited1 before any field sample: the new index guard
called `iinfo.max` as a function. Corrected the property access only; zero accepted
field results in that submission, no gate or geometry change. The hidden-device
CPU discovery warning is expected and is not new GPU fault evidence.


Measured compose table, fixed2mm pitch, margin0:

| quantity | result | unchanged gate |
|---|---|---|
| relative volume error | 0.0009750710958686022 | <=0.05 PASS |
| absolute volume error mm3 | 101.15500169417646 | <=75616 conservative raster bound PASS |
| sparse/dense occupancy mismatch | 0 | =0 PASS |
| deep independent parity disagreement | 0 | frozen field gate OK |
| all-probe parity disagreement | 0.041666666666666664 | retained surface-sensitive diagnostic |
| field shape | 34x26x58 | two complete fields identical |
| explicit coordinate weld | 4 duplicate vertices | all triangle-corner bytes unchanged |
| invalid input cases rejected | 8/8 | all required |

Evidence: `examples/anatomy/mesh_to_field_result.json` and
`examples/anatomy/mesh_to_field_arrays.npz`. Field SHA256:
`b1f261063e5043877162625ba63aed34bf4e5b0cc59cbec84fa1310bbbd3d7d9`.
The conservative boundary-AABB bound is intentionally loose (75616mm3); the
separate fixed5% band prevents that bound from being the sole quality criterion.
The seam validates connectivity/orientation, not arbitrary surface/tetrahedron
self-intersection or anatomical accuracy. Motion arrays are supplied in metres;
no triangle-contact engine integration has been certified here. The sampled
field is the frozen raster distance transform, not exact distance to triangles.

Run from this repository:

    pip install -r requirements-geometry.txt
    python examples/anatomy/inspect_mesh.py examples/anatomy/left_kidney_atlas.obj --units mm --out examples/anatomy/mesh_mechanism.json
    CUDA_VISIBLE_DEVICES='' OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python examples/anatomy/compose_mesh_to_field.py

The first command-line probe intentionally exits1 for the retained raw topology
negative. The compose exits0. `THREEFOLD_ROOT` optionally selects the parent of
the sibling staging repositories. No sibling file is modified. Next: preregister
photon transport reference cases and published calibration limits before design.

## Photon reference measurement (preregistered)

Before designing the new transport kernel, measure the external MCX `cube60`
benchmark using its own `mcxcreate` configuration unchanged, twice:1e6 photons,
60mm cube, pencil source, absorption0.005/mm, scattering1/mm, anisotropy0.01,
index1.37, reflection disabled, fixed seed1648335518 and5ns window.
Source: [official MCX configuration](https://github.com/fangq/mcx/blob/master/pmcx/pmcx/utils.py).
Fixed diagnostic gates: same configuration, finite nonnegative fluence, complete
fluence arrays bit-identical twice. Record any external floating-tally repeat
failure rather than weakening the exact-repeat criterion. Report integrated
absorption, field hash and maximum repeat difference. No speed claim in this
reference probe. The RT-MMC paper supplies stochastic agreement/leakage results,
not a universal5% acceptance guarantee for arbitrary tissue/source geometries.
A published calibration band for a760/850nm intensity ratio remains unverified;
it must not be invented or transferred from another wavelength pair.


Measured reference table before photon solver design:

| leg | integrated absorbed fraction | maximum flux | exact full-field repeat |
|---|---|---|---|
| 1 | 0.177973707130105 | 174313536 | FAIL |
| 2 | 0.177974013526867 | 174323760 | FAIL |

Packet energy absorption is0.1779742482680604 in both reference runs; integrated
voxel tally differs by3.063967620e-7 in absorbed fraction. Maximum flux difference
10224 corresponds to5.86655e-5 of the second field peak, but is not bit identity.
The external implementation is retained unchanged. Evidence:
`reports/mcx_reference_l4.json`, `reports/mcx_reference_fields_l4.npz`.
This motivates an integer energy/fluence tally for the separate Warp candidate;
it does not justify claiming the external reference deterministic.

Published reference gaps remain explicit: cube60 is the official MCX benchmark,
not the heterogeneous RT-MMC paper's cubesph case. A fresh reference comparison
may be measured against a preregistered engineering band, but that band cannot
be labelled a paper-specified tolerance. A raw760/850 intensity ratio is also
not the pulsatile ratio-of-ratios used in pulse-oximetry calibration.

### Warp photon candidate gates (before implementation/run)

Separate mesh-BVH transport, continuous absorption along scattering paths,
Henyey-Greenstein directions, face-adjacent region labels, Fresnel reflection/
refraction when enabled. Fixed per-photon seed and integer packet energy2^30;
integer voxel absorption sums. Split path deposits at voxel boundaries rather
than assigning a whole segment to its end voxel. Record all-time energy and
separately the first5ns field to match the external reference window. Continue
packets after5ns for all-time energy accounting; cap50000 physical events and
report residual packet weight separately. Never book an unexplained mesh miss or
out-of-grid point as a successful escape.

First full gate run: two1e6-photon cube60 configurations, mesh12 triangles,
properties/source matching the external table, reflection disabled as in MCX.
Require full integer tally/output arrays bit-identical; absorbed+escaped equals
launched exactly; zero event-cap residual and zero numerically leaked weight;
finite nonnegative field. Compare the5ns field to both frozen MCX arrays:
absorbed fraction within2% relative and5mm-binned fluence L1 relative error<=5%.
Those two comparison bands are preregistered engineering gates, not tolerance
values claimed to come from a paper. Record native event/query counts, boundary
hits and numerical leakage. This first cube test does not satisfy the separate
slab-diffusion, refractive/layered, or physiological calibration gates.

Initial Warp submission: code generation failed before launching photons because
a loop-index name was reused for a vector. Zero photon samples accepted; renamed
that new local variable only, with no formula or tolerance change.

### Oxygen-to-optics mechanism preregistration

Before a NIRS cell design, compose the existing `sao2_hill` and `HB_G_DL` exports
with the [published extinction table](https://omlc.org/spectra/hemoglobin/summary.html):
760nm oxygenated/deoxygenated586/1548.52,850nm1058/691.32 in cm^-1/M,
molar mass64500g/mol, natural attenuation factor ln(10). Explicitly convert
concentration g/dL to mol/L and inverse cm to inverse mm. Sample the frozen Hill
cell at40/60/80/100mmHg, twice. Fixed gates: identical tables; positive absorption
and saturation strictly between0 and1;760 absorption derivative negative and850
positive. The equal-path log(I760/I850) derivative is a mechanism diagnostic,
not a pulsatile calibration or a tissue scattering result.

Historical pre-relocation observation for `src/bodytwin/cells/optics/tissue_photon_mc_v1.py` (OWN-GATE-FAIL): L4 two1e6 runs exact, but578026 leaked packets and43.1298% absorption error; energy/leak/reference gates FAIL, unpromoted. Current module status is recorded for `src/bodytwin/geometry/optics/tissue_photon_mc_v1.py` below.
Historical pre-relocation observation for `src/bodytwin/cells/optics/__init__.py` (VERIFIED-FRESH): Imported in the full L4 cube probe. Current module status is recorded for `src/bodytwin/geometry/optics/__init__.py` below.

Measured first photon candidate table, both runs identical:

| quantity | result | gate |
|---|---|---|
| leaked photons | 578026/1000000 | zero FAIL |
| event-cap photons | 0 | zero PASS |
| absorbed integer energy | 108677877113643 | retained |
| escaped integer energy | 421091989344847 | retained |
| residual integer energy | 543971957541510 | zero FAIL |
| launched integer energy | 1073741824000000 | absorbed+escaped FAIL |
|5ns absorbed fraction | 0.10121415246804058 | reference0.17797370713010524/0.17797401352686726 |
| worst relative absorption error | 0.4312981403166415 | <=0.02 FAIL |
| worst binned fluence error | 0.4312981403166416 | <=0.05 FAIL |

All energy is accounted for when failed-packet residual is included; failed
packets were not relabelled escaped. Evidence: `reports/photon_cube_l4.json`
and `reports/photon_cube_fields_l4.npz`. Next: identify failure sites in a
separate observational copy before a numerical change; preserve this baseline.


Measured Hill/extinction seam, no tissue or clinical calibration claim:

| oxygen pressure | saturation | mu760/mm | mu850/mm |
|---|---|---|---|
| 40.0 | 0.750543721682 | 0.442367645034 | 0.517561889209 |
| 60.0 | 0.899914058416 | 0.365379943219 | 0.54689099458 |
| 80.0 | 0.951339882015 | 0.338874305892 | 0.556988537587 |
| 100.0 | 0.972761280059 | 0.327833397495 | 0.561194663476 |

Evidence: `reports/hemoglobin_mechanism.json`. Existing oxygen cell unchanged.

### Photon failure-site observer preregistration

Use a separate copy with distinct integer reason codes: region-label mismatch,
invalid ray distance, scoring-grid exit, incomplete segment loop. Record whether
a boundary hit was already known and the remaining segment length when stopped.
Do not change paths, energy or control predicates. Gates: two complete diagnostic
arrays identical; canonicalized terminal status and all original tally/counter
hashes match the frozen failed run exactly. No relaxed physical gate.

Historical pre-relocation observation for `src/bodytwin/cells/optics/tissue_photon_failure_probe.py` (VERIFIED-FRESH): L4 canonical field/terminal/counter hashes equal frozen failed variant exactly; two diagnostic arrays identical. Current module status is recorded for `probes/optics/tissue_photon_failure_probe.py` below.

Measured failure mechanism table (both1e6-photon legs identical):

| reason | photons | already had boundary hit | remaining length median mm | maximum mm |
|---|---|---|---|---|
| region label | 0 | 0 | 0 | 0 |
| ray distance | 0 | 0 | 0 | 0 |
| grid exit | 307890 | 307890 | 2.38e-7 | 6.3195e-5 |
| segment limit | 270136 | 21971 | 0.624459648 | 11.615131648 |

Evidence: `reports/photon_failure_sites_l4.json`. All grid exits occur while a
correct boundary is already known; segment-limit failures occur despite a1024
voxel-segment allowance. This localizes the defect to scoring traversal rather
than an unexplained BVH miss or region-label transition.

Next separate v2 design: parameterize voxel crossings from the fixed segment
origin in float64, advance integer voxel indices at each crossed grid plane,
and refine a reported mesh hit against that face's plane in float64. Avoid
repeated float32 position increments and look-ahead epsilon for selecting the
scoring voxel. Keep all original physical/reference gates,1e6 photons,50000
event cap and1024 segment cap. No acceptance tolerance changes; original v1
and the observational copy remain frozen. No performance claim.

Historical pre-relocation observation for `src/bodytwin/cells/optics/tissue_photon_mc_v2.py` (OWN-GATE-FAIL): L4 two1e6 runs exact; leaks578026->1, absorption error0.66677% and binned fluence error1.11931% PASS; zero-leak/exact escaped-energy gates still FAIL. Current module status is recorded for `src/bodytwin/geometry/optics/tissue_photon_mc_v2.py` below.

Measured v2 table, both full runs identical:

| quantity | value | gate |
|---|---|---|
| events | 83158852 | reported |
| boundary hits | 999999 | reported |
| leaked photons | 1/1000000 | zero FAIL |
| failed-packet residual energy | 480925106 | zero FAIL |
|5ns absorbed fraction | 0.17678734392561762 | relative error<=2% PASS |
| worst absorption relative error | 0.006667656573752762 | <=0.02 PASS |
| worst binned fluence L1 error | 0.011193055639850501 | <=0.05 PASS |

Evidence: `reports/photon_cube_v2_l4.json` and its field NPZ. Original v1 negative
retained. One failed photon is still failure under the unchanged zero-leak gate;
it is not waived using a published leakage rate. Next separate observer records
that packet's failure site and last position/direction/ray distance. Require
canonical field/terminal/counter hashes identical to v2 and complete diagnostic
arrays identical twice. No numerical changes in that observer.

Historical pre-relocation observation for `src/bodytwin/cells/optics/tissue_photon_last_failure.py` (VERIFIED-FRESH): Two full diagnostic arrays identical; canonical v2 field/terminal/counter hashes unchanged. Current module status is recorded for `probes/optics/tissue_photon_last_failure.py` below.

Last-packet mechanism before further design:

| photon | position mm | ray distance mm | boundary hit | residual integer |
|---|---|---|---|---|
|431164 |12.919497489,24.640228271,13.076808929 |0 |0 |480925106 |

This is an interior zero scattering distance, not an undetected tissue boundary.
Warp1.13's local `native/rand.h` returns a24-bit uniform on[0,1), including0.
The existing `-log(1-u)` therefore returns zero at that endpoint. Evidence:
`reports/photon_last_failure_l4.json`.

Separate v3 design: map that same24-bit draw to its bin midpoint by adding2^-25
in float64, then evaluate the optical-depth logarithm in float64 before casting
to float32. Endpoints become strictly interior without retries or additional
random draws. Keep the zero-distance rejection, all physical gates and both caps
unchanged; do not waive the failed packet. Retain v1/v2 and both observers.
The minimum resulting optical depth is2.980232283178454e-8, rather than zero.

### Next boundary/slab reference gates (before runs)

Measure unchanged MCX `cube60b` (refractive external boundary) twice and a
separate homogeneous slab approximation:100mm cube, pencil source[50,50,0],
mua0.01/mm, mus10/mm, g0.9, n1, reflection disabled,1e6 photons,5ns window.
Use the semi-infinite diffusion Green function and source/mirror convention of
[MCX cwdiffusion](https://github.com/fangq/mcx/blob/master/utils/cwdiffusion.m).
Fixed scoring ROI: radius10mm central cylinder, two-mm depth bins5..21mm;
compare voxel-centre averaged analytic fluence to the identical MCX voxel ROI.
Require maximum relative bin error<=5%, a declared engineering MC/diffusion band
for this geometry, not a universal band beyond one transport mean free path.
Also require finite nonnegative fields and report exact-repeat failure without
relaxation. No new transport implementation is certified by the external run.

Historical pre-relocation observation for `src/bodytwin/cells/optics/tissue_photon_mc_v3.py` (VERIFIED-FRESH): L4 cube six unchanged gates PASS, two1e6 full arrays identical, zero leaked/capped/residual packets and exact absorbed+escaped integer energy. Current module status is recorded for `src/bodytwin/geometry/optics/tissue_photon_mc_v3.py` below.

Measured v3 cube table (both legs identical):

| quantity | result | fixed gate |
|---|---|---|
| photons | 1000000 | >=1000000 PASS |
| leaked / capped photons | 0 / 0 | zero PASS |
| absorbed integer energy | 189853201239019 | exact accounting |
| escaped integer energy | 883888622760981 | sum=1073741824000000 PASS |
| residual integer energy | 0 | zero PASS |
| events / boundary hits | 83158057 / 1000000 | reported |
|5ns absorbed fraction | 0.17678744561749604 | reference error0.66671% PASS |
| coarse fluence L1 error | 0.01119298472586431 | <=0.05 PASS |

Evidence: `reports/photon_cube_v3_l4.json` and its field NPZ. Integer absorption
SHA256 `202bc4c21c8718277683dc6149356c6f6e51f70839f128fb27b879dba4d433e1`.
All preceding negatives remain unchanged. This validates only the declared
homogeneous nonreflecting cube at the fixed seed, not arbitrary tissues,
refractive interfaces, diffusion, clinical calibration, or an RT speedup.
The module's direct script entry runs this same full gate suite.

First boundary-reference submission: two external cube runs printed27.31924%
absorption, then the slab setup rejected an empty detector list before launching
slab photons. No complete reference report or slab sample was accepted. Omit the
optional detector key for the detector-free slab; optical/source/tally parameters
and all gates are unchanged. This is an API setup correction, not a physics fix.


Boundary reference table:

| case | absorbed fraction leg1 | leg2 | maximum full-fluence repeat difference |
|---|---|---|---|
| cube60b | 0.273191075811929 | 0.273190536874355 | 0.000106930732727051 |
| slab | 0.25427508157572 | 0.254270217663437 | 0.000484257936477661 |

Evidence: `reports/mcx_boundary_references_l4.json` and the corresponding NPZ.
The fixed slab ROI has maximum0.502268% diffusion error; the chosen approximation
is measured for this source/material/ROI only. Fresh external cube60b energy
normalization uses1e6 incident-in-tissue packets (normalizer200), matching the
new kernel's inside-tissue launch convention. No source-energy scale is fitted.

Regression finding: the existing CPU cell-discovery test also collects every
experimental transport backend placed under `cells/optics`; a filtered run of
the v1 library fails1 test because it is not a standalone physiology cell.
Do not weaken that test or relabel the failed transport as a passing cell.
Place experimental solvers under `geometry/optics` instead, preserving each
backend file byte-for-byte; companion probe imports follow that package move.
The actual future NIRS physiology cell remains a separate integration step.
Record moved source hashes and verify the certified cube hashes after relocation.

### Boundary certification and package-regression gates (before run)

The imported backend package is now `bodytwin.geometry.optics`; all six backend
files, including every frozen failed comparator, keep their exact source bytes
(`reports/optics_package_move.json`). The CPU physiology-cell tests remain
unchanged. Companion probes only change their package imports.

Run the relocated frozen v3 at1e6 photons twice for each of: original cube60,
refractive cube60b, and the measured100mm slab. Original cube full integer hashes
must match the pre-move certified hashes. All cases require exact repeat arrays,
absorbed+escaped=launched, zero residual/leaked/capped photons, finite nonnegative
fluence. Refractive cube and slab each require absorbed fraction within2% and
5mm-binned fluence L1 within5% of both fresh external references. Slab's declared
central-cylinder profile must be within5% of the analytic profile and both MCX
profiles in every fixed bin. No gate/source/ROI/seed/cap tuning after the run.

### Extended boundary result (two complete L4 runs per case)

| case | failed packets | capped packets | residual integer energy | max absorption error | max binned fluence error | gates passed |
|---|---|---|---|---|---|---|
| cube60 | 0 | 0 | 0 | 0 | 0 | 5/5 |
| cube60b | 1467 | 0 | 149320786816 | 0.00368212586529 | 0.00867129892313 | 4/6 |
| slab | 6 | 860 | 4543095627 | 0.00341040140843 | 0.0084754226389 | 6/8 |

All full repeat arrays are bit-identical. Original cube hashes equal the pre-move
certified hashes exactly. The refractive cube and slab FAIL exact absorbed+escaped
energy and zero failed/residual/capped gates. Slab maximum analytic-profile error
is 0.011594185751949928; maximum MCX profile error is 0.015030295983538233,
both below the preregistered 0.05. These passing profile comparisons do not waive
failed packet accounting. No general boundary certification or promotion.
Evidence: `reports/photon_boundary_l4.json`. No throughput claim.


CPU regression after relocation: 9 passed, 228 deselected in the existing cell-graph,
metabolic-chain and blood-oxygen subset, with CUDA hidden before imports. No test
code changed. The earlier library-as-cell discovery failure remains recorded.

Next measurement: a separate observer must preserve all canonical v3 hashes and
repeat arrays while counting failure sites and residual/event distributions for
the refractive cube and slab. Both observational gates must pass before a new
transport design. The original 50000 event cap and zero-failure gates stay fixed.




Arithmetic measurement gates, before execution: two CPU tables exactly equal and
finite loss thresholds. Measure the frozen direction-offset expression in float32
at the fixture boundaries and per-segment rounding thresholds at slab properties.
This diagnostic alone does not attribute actual photon failures.

Arithmetic table (CPU, two exact repeats):

| mechanism | measured value |
|---|---|
| direction offset at coordinate60 or100, normal component -0.1,-0.01,-0.001 | represented displacement0 in all6 cases |
| slab weight1, segment0.1, mua0.01 | continuous loss0.000999500166625; rounded loss0 |
| slab weight32, segment0.1 | rounded loss0; first nonzero loss needs segment1.574835696814 |
| slab weight1 | first nonzero loss needs segment69.314718055995 |


### Frozen-v3 failure attribution

| case/site | count | measurement |
|---|---|---|
| refractive cube / ray distance | 1467 | All known boundary hits at distance0; sampled positions lie exactly on coordinate60 faces with inward reflected direction. |
| slab / event cap | 860 | All50000 events; residual weights66..81 (median73), sum63144; zero-loss segments40814..41692, trailing maximum16721. |
| slab / ray distance | 1 | Boundary distance0 at position x=100, outgoing direction; residual243542. |
| slab / grid exit | 5 | No boundary hit reported; residual4542788941; sampled surface distances range near1e-7 to0.01455. |
| region labels / negative segments / segment cap | 0 | No observed failures at these sites in either fixture. |


Separate v4 design after these tables: preserve continuous float64 packet weight
across segments, quantize remaining weight cumulatively, and tally its integer
difference. Thus sub-integer absorption carries forward instead of being discarded
each segment. Partial-window scoring uses the same continuous prefix weight.
At handled interfaces move along the oriented surface normal into the selected
region using the existing1e-5 displacement magnitude, rather than projecting that
magnitude onto a possibly grazing ray. This does not certify arbitrarily scaled
meshes and does not fix or forgive unreported boundaries.

Gates before v4 execution: exactly the same full cube/refraction/slab suite,
1e6 photons twice each, zero leaks/residual/caps, exact integer conservation,
finite nonnegative fluence, external absorption2%, binned fluence5%, slab profiles5%.
The original cube hash comparison is specific to the package-only move; v4 is a
new algorithm and is compared by these unchanged physical gates and exact repeats.
All v1/v2/v3 files and failed results remain frozen. No cap or tolerance increase.




Ray-replay measurement gates before execution: two CPU replay tables exactly equal,
all six sampled slab failure rays have a double-precision triangle hit. Compare
bounded/unbounded Warp queries to a strict double triangle test. Diagnostic positions
were truncated at1e-9, so this replay is not a canonical full-photon-path claim.

Ray replay table (CPU; rounded diagnostic positions, not original full paths):

| photon | bounded/unbounded Warp hit | double distance |
|---|---|---|
| 343906 | False/True | -0 |
| 128514 | False/False | 1.959779207054538e-07 |
| 410895 | True/True | 6.950915921212805e-07 |
| 434896 | True/True | 0.02389985595141492 |
| 806850 | True/True | 4.208978223224778e-07 |
| 864837 | True/True | 0.005010847277854945 |

Both replay gates PASS. One replay ray misses even with unbounded Warp query,
while the strict double triangle test finds distance1.9597792070545384e-7.
Two rays are reported at signed zero instead of their positive double distances.
CPU replay differs from some original GPU observations; no identity inference is
made across rounded inputs or execution devices.


### Cumulative-weight and normal-offset v4 result

| case | leaks | event caps | residual integer | max reference absorption error | gates |
|---|---|---|---|---|---|
| cube60 | 0 | 0 | 0 | 0.0066670848581938 | 6/6 |
| cube60b | 1 | 0 | 567818321 | 0.0035944784524388 | 4/6 |
| slab | 6 | 0 | 4543032486 | 0.0034103943636651 | 6/8 |

Two complete arrays/results per fixture are bit-identical. The original cube passes all six unchanged physical gates. Refractive failures
fall1467->1 (residual567818321), so its accounting still fails. Slab event caps
fall860->0, but all six
geometric failures remain; exact energy and zero-failure gates still FAIL.
The two measured fixes do not resolve the separate slab ray-query defect.
No general transport promotion or speed claim; photon_boundary_v4_l4.json.


Documentation correction: the first v4 prose summary incorrectly said both cubes
passed. The measured table and JSON always record one refractive failure; the
summary is corrected above. No failed gate was changed.

### Double-precision triangle candidate variant, gates before execution

The v3 site table and scalar replay show boundary-distance zero/misses in the
float32 ray path; v4 independently removes caps but retains6 slab and1 refractive
failures. Build v5 beside v4: Warp BVH AABB queries identify candidate triangles;
strict double-precision ray/triangle tests choose the nearest hit. Keep position
and remaining scattering depth in double precision. Refine a selected distance
against the stored face plane as before. A conservative float32 broad-phase box
padding of2^-21*(maximum endpoint magnitude+1) only expands candidate retrieval;
it never changes the exact triangle acceptance or physical scoring interval.
The same refractive, absorption, source, cap, seed and physical tolerances apply.

Run the complete three-fixture physical suite twice at1e6 photons, same six/eight
gates as v4. Zero failed/residual/capped packets and exact integer energy remain
required; no cube-specific analytic boundary substitution. No speed claim until
separate idle timing legs exist. This is still a software BVH backend, not RT.




### V5 negative result

| case | leaks | residual integer | boundary hits | gates |
|---|---|---|---|---|
| cube60 | 999993 | 883888491253633 | 19 | 4/6 |
| cube60b | 999994 | 883888546645509 | 19 | 2/6 |
| slab | 986347 | 801647146974389 | 0 | 6/8 |

All repeat arrays exact and event caps0; physical accounting fails severely.
Refractive reference absorption/fluence error0.3528798612794555 also FAIL.
The near-total absence of boundary hits requires an accessor-contract measurement
before interpreting this as failure of double-precision geometry.


Accessor measurement gates before execution: two complete CPU vertex arrays exact;
direct face-corner access must equal all36 source triangle corners byte-for-byte.
Also measure the candidate's extra index-indirection expression without changing it.

Accessor result table (two identical complete CPU arrays):

| expression | source corners mismatched /36 |
|---|---|
| mesh_get_point(mesh,corner) | 0 |
| mesh_get_point(mesh,mesh_get_index(mesh,corner)) | 36 |


The installed Warp1.13 native accessor already applies mesh.indices internally.
V5 passed an already dereferenced vertex index, so its exact triangle test used
incorrect triangles. Separate v6 changes only the three corner-access expressions
and the companion entrypoint. V5 stays frozen. Before running v6, retain every
v5 physical gate, all sources/parameters, two1e6 photons per fixture and exact
integer energy/zero-failure requirements. No tolerance or geometry repair.




### V6 full boundary result

| case | absorbed integer | escaped integer | residual / leaks / caps | max reference absorption error | max binned fluence error | gates |
|---|---|---|---|---|---|---|
| cube60 | 189853201327492 | 883888622672508 | 0 / 0 / 0 | 0.0066670847211908 | 0.011192983020515 | 6/6 |
| cube60b | 292340208610692 | 781401615389308 | 0 / 0 / 0 | 0.0035943196245712 | 0.0086721706779167 | 6/6 |
| slab | 272094677025611 | 801647146974389 | 0 / 0 / 0 | 0.0034103760889518 | 0.0084754398510559 | 8/8 |

All20 fixed case gates PASS, two1e6 full arrays per case bit-identical. Each
integer absorbed+escaped sum equals1073741824000000 exactly. Maximum slab
analytic-profile error0.01159414463830433, maximum external
profile error0.01503025501282079, both below0.05.
Evidence: `reports/photon_boundary_v6_l4.json`; Warp1.13 on L4. Correctness only,
no certified throughput/idle timing or hardware RT claim.


Remaining scope: the fixtures are homogeneous with an exterior refractive boundary.
Internal region labels are supported by the interface but heterogeneous-layer
transport is not certified yet. Next measure an external layered reference before
a layer/detector adapter; then connect760/850 absorption to the existing SaO2 cell
and propagate declared uncertainty in a chain. The brief's published calibration
slope requirement remains open: these engineering fluence bands are not a pulse
calibration certificate. A shared hardware RT backend and idle timing are also
unmeasured. Earlier failed variants and exact-repeat-negative external fields
remain recorded; v6 success does not retroactively relabel them.

### Layered reference measurement, gates before execution

Continue from certified v6; no backend edits. External MCX cube60b is changed only
by a z=10 split in the60-cube, a second optical row, and one radius2 detector at
(39,29,0). Region1 is (mua .005, mus1, g .01, n1.37); region2 is
(.01,2,.8,1.4). Launch/source/seed/1e6 photons/5ns/exterior reflection stay fixed.
Record integrated voxel fluence, per-layer absorption, detector partial paths and
weights before designing the mesh-layer seam. The external weight helper is
checked against exp(-sum(mua*partial_path)) with fixed1e-12 maximum difference.
Other preregistered gates: finite nonnegative arrays, absorbed fraction in(0,1),
nonempty detector and exact external full-field repeat. Run twice, preserve any
float-atomic repeat failure. No tolerance tuning or clinical calibration inference.
Detector counts and path statistics are measurements; no fitted detector acceptance.



First layered reference submission completed one external1e6 simulation
(absorbed28.58477%,4217 detector photons in the external log), then failed Python
postprocessing because pmcx.run returns a raw detector matrix rather than the
parsed dictionary. Zero completed two-leg reference reports were accepted. Use
the official detphoton parser with2 media and the unchanged dpxv flag. This is
an API decoding correction; source, geometry, physics, counts and gates stay fixed.

Layer geometry measurement gates before execution: two exact CPU tables, both
input boxes closed/consistently oriented, analytic volumes36000 and180000 exactly.
Count duplicate coordinates and the four input interface triangles before any
region-labelled seam construction. This geometry measurement does not certify
photon transport or detector response.

Layer source geometry table (two exact CPU runs):

| quantity | region1 | region2 |
|---|---|---|
| closed / consistent orientation | True / True | True / True |
| volume |36000|180000|
| interface triangles / area |2 /3600|2 /3600|
| interface normal |+z|-z|

Combined input16 vertices has12 unique coordinates. All4 interface triangles are
geometrically distinct because opposite box faces use different diagonals; simply
removing duplicate faces cannot construct the seam.


Separate labelled-box adapter design: keep both sets of exterior faces and just
the first region's interface triangulation, labelled back1/front2. Exact-coordinate
welding merges the4 repeated interface vertices. For each region, assemble its
outward boundary using front/back labels and check it with the unchanged strict
surface validator. Preregistered CPU gates: full arrays repeat exactly; both region
boundaries closed and consistently oriented; exact analytic volumes36000/180000;
12 vertices/22 faces with exactly2 shared interface triangles and area3600.

### Layered external and geometry results

| measurement | first leg | second leg |
|---|---|---|
| absorbed fraction |0.2858475061524023|0.2858476059011319|
| first-layer absorption |0.17706063701499675|0.17706073678792114|
| second-layer absorption |0.1087868691374055|0.10878686911321071|
| detector photon count |4217|4217|
| detected fraction |0.00281081862729486|0.00281081862729486|
| detected paths visiting layer2 |669|669|
| weighted layer1/layer2 path |62.61902390723209 /2.188499103138981|same|
| independent/helper weight max difference |1.1102230246251565e-16|same|

External full-field repeat FAIL, maximum difference1.9669532775878906e-5.
All four other external gates PASS. Detector aggregate repeats are identical;
raw detector path row hashes differ, so no full raw-detector repeat claim.
The launch position is the unmodified cube60b source(29,29,0); the report's source
key holds its upstream-code URL. Evidence: mcx_layer_reference_l4.json and NPZ.


Preregistered layered transport gates before execution: frozen v6 plus the measured
labelled mesh, two1e6 complete outputs exactly equal; exact absorbed+escaped energy;
zero leaked/residual/capped packets; finite nonnegative fluence; total absorption
within2% and each layer's absorption within2% of both external runs;5mm-binned
fluence L1 within5% of each reference. Same source, optical rows,5ns,seed and caps.
These are the existing engineering field/energy bands applied to the new fixture,
not a clinical or detector-response certificate. No backend change in this step.



### Detector observer design and gates before execution

The external detector table gives4217 paths,669 visiting layer2, detected fraction
0.00281081862729486 and weighted mean paths62.61902390723209/2.188499103138981.
A separate v6 observational copy will record per-region path lengths and each
escaped packet's position, outgoing direction and clock. It must preserve the
frozen layered transport's three canonical hashes and repeat every full array
exactly twice. No new random draw, termination rule or transport parameter.

Fixed observer/detector gates: canonical hashes exact; all arrays repeat exact;
nonempty radius2 detector at(39,29,0), accepting only exits within5ns; reconstructed
escaped integer energy from partial paths differs by at most1 integer unit per
escaped photon; accepted exits satisfy abs(z)<=1e-8. Compare detected fraction to
each external reference within5 combined estimated Monte Carlo standard errors,
computed from per-launched-photon weights including zero weights for nondetection.
Compare both weighted mean layer paths within5 combined ratio-estimator standard
errors. These statistical gates are fixed before execution, are not a clinical
calibration band, and must not replace the layered field/energy gates.




### Frozen-v6 layered transport result

| quantity | value |
|---|---|
| absorbed_integer | 305851234493267 |
| escaped_integer | 767890589506733 |
| residual_integer | 0 |
| leaks | 0 |
| event_caps | 0 |
| events | 164624214 |
| boundary_hits | 3653758 |
| window_absorbed_fraction | 0.284845973468123 |
| layer_absorbed_fractions | [0.17651375285443383, 0.10833222061368916] |
| maximum total absorption error | 0.003504078440157611 |
| maximum per-layer absorption error | 0.004179259200318404 |
| maximum binned fluence error | 0.008249553977462782 |


Detector energy reconstruction uses the actual float32-stored optical rows cast
back to float64, matching the frozen kernel's coefficients rather than decimal
input spelling. The one-integer reconstruction gate is unchanged. External
reference weights still follow the external helper's declared property table.

### Detector observation result

| quantity | candidate | external reference | fixed comparison bound |
|---|---|---|---|
| detected packets |4119|4217|count reported, no fitted count gate|
| detected fraction |0.0027638959100153297|0.00281081862729486|absolute5-SE bound0.00032035569504741294|
| layer1 weighted mean path |61.27105205513939|62.61902390723209|absolute5-SE bound3.9853409593017917|
| layer2 weighted mean path |2.2617713692072168|2.188499103138981|absolute5-SE bound0.7546786524763944|
| integer energy reconstruction max error |0|not applicable|<=1|
| detector exit abs(z) maximum |4.2229452272161743e-16|1.1008133071754855e-7|candidate<=1e-8|

All7 detector/observation gates PASS, including unchanged canonical layered hashes
and two complete absorption/terminal/counter/path/exit arrays exactly equal.
Detected fraction difference4.6922717279530065e-5 is within the preregistered
statistical bound against each reference. Evidence: photon_detector_l4.json and
photon_detector_paths_l4.npz. No clinical calibration claim.


Before a wavelength model, measure path-reuse completeness: with fixed scattering
and refraction, attenuation does not draw random numbers or change trajectories
until integer weight reaches zero. Compute the largest possible baseline optical
depth before5ns from max(mua/n)*c*window and compare to log(2*PACKET), the cumulative
rounding-to-zero threshold. Gates: exact two CPU tables, the conservative depth
bound strictly below that threshold, all saved detected paths within that bound.
This proves no baseline absorption termination before the detector time window;
it does not claim completeness after5ns or under changed scattering/refraction.

Path-reuse mechanism table, two exact CPU runs:

| quantity | value |
|---|---|
| maximum baseline depth before5ns |10.706873443019695|
| cumulative rounding-to-zero depth |21.487562597358306|
| minimum continuous weight before5ns |24041.62480141925|
| maximum saved detected-path depth |9.39673976567699|
| maximum saved detected-path time |4.9847370984031425ns|


### Synthetic two-wavelength forward model, gates before execution

Use the certified fixed-window detector path library with unchanged scattering,
refraction, geometry and detector. The existing Hill oxygen cell supplies SaO2
at PO2=40,60,80,100; its Hb=15 supplies whole-blood absorption from the previously
measured published760/850 extinction coefficients. Declare synthetic layer blood
fractions0.02/0.04 and background absorption0.002/0.004; these are engineering
fixture inputs, not measured anatomy or clinical priors. Intensity is the sum of
exp(-sum(layer_mua*path)) divided by the original1e6 launched photons.

Place the experimental forward model with optical solvers, outside certified
physiology-cell discovery until clinical calibration is resolved. CPU gates:
complete weight arrays/reports repeat exactly; finite intensities in(0,1); ratio
I760/I850 strictly increases over the fixed oxygen grid; analytic d(log ratio)/dS
positive; central-difference relative derivative error<=1e-7 with step1e-5.
Clinical slope calibration is separately marked unverified:0 accepted matching
calibration anchors, not silently replaced by the derivative identity.

### Synthetic NIRS forward result

| PO2 | SaO2 | I760/I850 | d(log ratio)/dS |
|---|---|---|---|
| 40.0 | 0.7505437216819002 | 1.083025935881287 | 0.7614437782953757 |
| 60.0 | 0.899914058415897 | 1.215323496343131 | 0.7823739172783407 |
| 80.0 | 0.951339882014719 | 1.265481748658765 | 0.7905730101195656 |
| 100.0 | 0.9727612800585383 | 1.287144663618554 | 0.7941639581271539 |

All5 engineering gates PASS; full weight arrays and results exact twice.
Maximum derivative relative error1.6337799058914892e-10 <=1e-7. The ratio rises
1.0830259358812873->1.2871446636185537 over the fixed oxygen grid. This is a
synthetic raw-intensity forward model; matching clinical calibration anchors0,
clinical calibration not certified.


Direct transport validation gates before execution: for each of the four oxygen
values and both wavelengths, run the unchanged path-recording transport twice at
1e6 photons with only absorption set by the new model. Use actual float32-stored
coefficients in the path prediction. All full arrays must repeat exactly, all
energy/accounting/no-failure gates remain unchanged, finite nonnegative arrays.
Compare summed detected integer energy with rounded path prediction: absolute
difference <=4119 integer units (one per original detected path), corresponding
to3.836117684841156e-12 of launched energy. Require direct I760/I850 to increase
over the same four oxygen points. No stochastic band is used for this same-path
comparison; no clinical-calibration inference or timing claim.



### NIRS chain sensitivity measurement, before chain design

Use the frozen forward model and existing Hill cell at synthetic nominal
PO2=60, Hb=15, blood-fraction multiplier1. Declare independent engineering input
bands PO2 +/-2, Hb +/-0.75 and blood-fraction multiplier +/-0.1. These are
illustrative input uncertainties, not published clinical confidence intervals.
Measure each input alone and all27 combinations of lower/nominal/upper values,
with the certified4119-path library and fixed scattering/refraction/background.
Record nominal ratio, each isolated span, joint corner span and whether the joint
span exceeds their sum; retain that composition negative if it occurs. Also
measure finite-difference sensitivities and Hb*fraction degeneracy.
Preregistered observation gates: exact two full tables; finite positive outputs;
normalized Hb and fraction sensitivities agree within1e-7 relative. These
observation gates do not imply the separate composition inequality passes.

### NIRS sensitivity mechanism table

| quantity | value |
|---|---|
| PO2 isolated ratio span |0.015442893677285907|
| Hb isolated ratio span |0.018323087693881135|
| blood fraction multiplier isolated span |0.0366558708431719|
| joint27-point span |0.07036612253119023|
| isolated span sum |0.07042185221433894|
| composition excess | -5.572968314870863e-5 |
| normalized Hb/fraction sensitivity relative difference |0|

All3 observation gates PASS, two full27-point tables exact. The sampled-grid
composition inequality passes. Normalized Hb and blood-fraction sensitivities
are identical; this forward model depends on their product and does not identify
them separately.


### Conditional uncertainty chain, gates before execution

Build a separate chain of the unchanged Hill cell, extinction/tissue mapping and
fixed-window detector response. Draw256 independent uniform samples of the three
declared input bands with seed20260912; repeat the entire chain twice. Reuse the
same draws with one stage perturbed at a time, as in the existing chain pattern.
No additional physiological coupling or calibration is fitted.

Preregistered gates: full input/intermediate/output arrays repeat exactly; finite
positive physical outputs; zero bands return nominal exactly; zero blood fraction
gives ratio1 exactly; increasing oxygen raises the ratio; joint Monte Carlo span
<= sum of three isolated Monte Carlo spans; the unchanged measured27-point joint
span <= isolated span sum; all random outputs inside the measured grid extrema.
Report the conditional central95% sample interval. These numerical composition
checks cover only the declared synthetic inputs; geometry, extinction errors and
finite-photon path-library uncertainty are not included and clinical calibration
remains false. They are not a rigorous global or clinical confidence certificate.

### Direct NIRS transport and conditional chain results

| PO2 | direct I760/I850 | maximum integer prediction error across wavelengths/repeats |
|---|---|---|
| 40.0 | 1.083025984235581 | 0 |
| 60.0 | 1.215323488776917 | 0 |
| 80.0 | 1.265481720608977 | 0 |
| 100.0 | 1.287144663899917 | 0 |

All6 direct-transport gates PASS over16 complete1e6 runs. Every pair repeats all
arrays exactly, all energies balance exactly with zero leaks/residual/caps.
Every detected integer sum equals its path prediction exactly (error0, fixed
bound4119). The full transport therefore verifies the reused-path wavelength
response for these fixed synthetic inputs. nirs_direct_transport_l4.json.

Conditional chain:8/8 preregistered gates PASS, two complete256-draw chains exact.
Nominal ratio1.215323496343131; joint sample span0.055127841665889354 versus isolated
sum0.06987483101741154. Central95% conditional sample interval
[1.1916161693493625,1.238607335149803]. Measured27-point span0.07036612253119023
versus isolated sum0.07042185221433894; zero-band and zero-blood nulls exact.
No clinical confidence or global-bound inference from these finite samples.


Outstanding: matching clinical slope calibration anchors0; optical/Hb blood-fraction
nonidentifiability retained; geometry, extinction-coefficient and finite-photon
uncertainty excluded from the conditional chain. The experimental forward module
is not promoted into certified clinical physiology cells. Next shared-backend work
starts with the unmodified headless OptiX toolchain proof, then measured ray-query
seam comparisons. No RT throughput/correctness claim yet.

### Earlier module attempts

Earlier attempt: `probes/optics/tissue_photon_boundary_failure.py`; CUDA-ONLY; Observer preregistered: preserve all three canonical hashes of each frozen boundary fixture, full diagnostic repeats exact, two1e6 runs each; result pending..

Earlier attempt: `probes/optics/photon_boundary_failure_probe.py`; CUDA-ONLY; Same preregistered observation gates; no changed physical gates or transport parameters..

Earlier attempt: `src/bodytwin/geometry/optics/tissue_photon_mc_v4.py`; OWN-GATE-FAIL; Separate variant; fixed boundary/slab gates above, execution pending..

Earlier attempt: `probes/optics/photon_boundary_v4_probe.py`; OWN-GATE-FAIL; Same reference fields, tolerances, event cap and two1e6 runs per fixture..

Earlier attempt: `src/bodytwin/geometry/optics/tissue_photon_mc_v5.py`; OWN-GATE-FAIL; Separate candidate, fixed full physical suite preregistered above; no accepted photon samples yet..

Earlier attempt: `probes/optics/photon_boundary_v5_probe.py`; OWN-GATE-FAIL; Unchanged physical gates, two1e6 runs per case; pending..

Earlier attempt: `src/bodytwin/geometry/optics/tissue_photon_mc_v6.py`; CUDA-ONLY; Corrected accessor variant, preregistered full physical gates; no accepted photon samples yet..

Earlier attempt: `probes/optics/photon_boundary_v6_probe.py`; CUDA-ONLY; Unchanged two1e6 runs per fixture, sources/caps/tolerances; pending..

Earlier attempt: `probes/optics/mcx_layer_reference.py`; OWN-GATE-FAIL; Two external layered measurements preregistered above; pending..

Earlier attempt: `probes/optics/photon_layer_probe.py`; CUDA-ONLY; Seven preregistered layered gates above, two1e6 runs; pending..

Earlier attempt: `src/bodytwin/geometry/optics/tissue_photon_paths_v1.py`; CUDA-ONLY; Separate observer; seven preregistered detector/identity gates, pending layered baseline result..

Earlier attempt: `probes/optics/photon_detector_probe.py`; CUDA-ONLY; Two1e6 full observations with unchanged fixture; not run yet..

Earlier attempt: `probes/optics/nirs_direct_transport_probe.py`; CUDA-ONLY; Six preregistered validation gates over16 full1e6 runs; pending..

### Transport evidence reconciliation

The frozen v1 failure578026/1000000 and43.1298% absorption error remains an
OWN-GATE-FAIL comparator. The separate v6 backend already passes cube6/6,
refractive cube6/6 and slab8/8 with zero leaks/caps/residual and identical complete
hashes across both runs. Reinspection of photon_boundary_v6_l4.json confirms all20
gates and zero failed packets. Layer transport7/7, direct NIRS transport6/6 and
conditional NIRS chain8/8 recorded gates also pass. This audit reads saved evidence;
it is not a new GPU run. Full numeric tables above remain authoritative.

The NIRS fixture is a synthetic two-layer benchmark, not an accepted layered
anatomical head or forearm. Clinical slope calibration remains uncertified.
Those scope gaps are still open; passing transport does not close them. Existing
v1 and intermediate failed variants remain unchanged for comparison.

### NIRS observable contract gates before execution

The frozen forward model emits static I760/I850. A published primary study defines
pulse modulation through AC/DC and an inter-wavelength ratio:
https://pmc.ncbi.nlm.nih.gov/articles/PMC6211094/ . These are distinct observables;
a calibration for one cannot be transferred solely because both depend on oxygen.
Before designing a pulse model, measure a counterexample using frozen paths and
four existing oxygen inputs. Apply detector gains0.5/1/2 at760 while850 stays fixed.
For an explicit synthetic blood-fraction perturbation of+/-1%, compare static
ratios with modulation ratios (I_low-I_high)/mean(I_low,I_high) per wavelength.
This is an algebraic observer, not a physiological pulse model or clinical dataset.

Frozen gates: two complete tables/weight hashes exact; static ratio scales exactly
with detector gain; normalized modulation ratio is exactly gain-invariant; at
least one accepted matching calibration anchor exists. Last gate is expected to
FAIL with0 current anchors; do not substitute a different observable or tolerance.


Measured counterexample at oxygen input40: static ratio0.541512968/1.083025936/
2.166051872 at gains0.5/1/2, while normalized modulation ratio remains exactly
0.886374669. At input100 static ratio0.643572332/1.287144664/2.574289327,
normalized modulation remains0.658471923. All four oxygen cases obey exact gain
scaling/invariance in both runs; no clinical anchor accepted. Evidence
reports/nirs_observable_contract.json. This preserves the failed calibration gate
and prevents applying a pulse-ratio slope bound to the existing static ratio.

### Curved two-layer fixture mechanism gates before execution

Before a new region-labelled curved fixture, measure concentric closed cylinder
surfaces with outer radius30/length120 and inner radius27/length114 at32/64/128
angular sections. These are declared synthetic shape parameters, not measured
anatomy. Observe closedness/orientation, analytic volume error, shell volume,
radial polygon sagitta and minimum radial/end clearance. Three observer gates:
two full array/metric hashes repeat, both individual boundaries closed/oriented,
positive shell volume and strictly positive conservative nesting clearance.
Predeclare candidate geometry volume-error bound0.002 per region; retain failing
coarse tessellation numbers. No transport or anatomical fidelity claim.


| angular sections | max region relative volume error | conservative clearance | candidate volume gate |
|---|---|---|---|
|32|0.006413148855795|2.855541800166|FAIL|
|64|0.001605606964382|2.963863686155|PASS|
|128|0.000401546850321|2.990964560886|PASS|

Evidence reports/curved_layer_mechanism.json, complete geometry hashes repeat.
Before design choose64 as the first measured resolution meeting unchanged0.002.
New synthetic curved two-layer fixture uses outer front0/back1 and inner
front1/back2, each interface once. Gates before transport: both reconstructed
region boundaries closed/oriented, volume errors<=0.002, all arrays repeat.
Transport v6 stays unchanged, two1e6 runs, gates exact absorbed+escaped energy,
zero residual/leaks/caps, full arrays repeat and finite nonnegative absorption.
No MCX equivalence, detector/clinical or measured-anatomy gate asserted for this
curved synthetic transport control. Existing slab and MCX certificates stay with
their own fixtures. GPU only on cloud; no local GPU reservation required.


Evidence reports/photon_curved_layer_l4.json; frozen v6 source hash unchanged.
Cloud simulation exit0. Bulk local artifact download failed with storage errno28;
the small result JSON was recovered separately from the persistent cloud volume,
without rerunning simulation. This storage failure is not a transport gate failure.

### Curved detector and NIRS observer gates before execution

Use the unchanged path-recording backend with the same curved fixture/source/
properties/seed, two1e6 runs. Detector radius2 centred at[42,32,2] with5ns window;
geometry is a declared synthetic curved control. Five observer gates: canonical
transport hashes equal frozen curved transport, full path/exit arrays repeat,
nonempty detector, reconstructed escape energy error<=1 integer (same earlier
path gate), and finite positive NIRS ratios increasing across the four frozen
oxygen inputs. No external detector or calibration gate claimed here. Record
all results, including empty-detector failure if it occurs, without moving it.


Curved static ratios at the four frozen oxygen inputs:1.0914721238114424,
1.2399585366270993,1.2971776200085616,1.3220532990516518. Layer mean paths
27.25023186014166/25.87265433760058. Reports photon_curved_detector_l4.json and
photon_curved_detector_paths_l4.npz. No external detector or clinical claim.

### Curved direct NIRS gates before execution

Before direct transport, measured maximum optical depths on3590 saved paths at
oxygen40/60/80/100 are25.4262499521/26.6337874914/27.0495233988/27.2226979584;
minimum rounded packet weight is0 in every case. The earlier box no-rounding-zero
argument does not transfer; retain this negative rather than assuming identical
live paths under changed absorption. Direct transport comparison is required.

Separate probe: four oxygen values/two wavelengths/two full1e6 runs, frozen v6
path backend and curved geometry/detector. Six gates before running: full repeats,
exact total energy, zero failed packets, finite nonnegative outputs, predicted
versus direct detector integer error<=3590 (one per original detector path, the
same earlier rule), increasing wavelength ratio. No clinical slope certificate.
No changed geometry/source/seed or tolerance after execution.

### Curved uncertainty mechanism gates before execution

Before a curved-fixture chain, retain nominal inputs[60,15,1] and half-bands
[2,0.75,0.1] from the frozen box chain. Measure27 grid points/isolated spans and
finite-difference sensitivities twice on the saved curved detector paths. Three
observation gates unchanged: tables repeat exactly, finite positive ratios,
Hb/fraction relative sensitivity difference<=1e-7. Composition is a measured
outcome, not forced by adjusting bands. Excludes anatomical and clinical errors.


Measured before chain: nominal ratio1.2399585366270993; derivatives with respect
to oxygen input/Hb/fraction0.0043763409670/0.0122624885165/0.183937327747.
Composition excess-0.0000704806778966. Hb and fraction remain unidentifiable
separately through this signal; no clinical uncertainty interpretation.

Separate curved chain uses frozen propagate/leg functions, seed20260912,
256 draws, bands[2,0.75,0.1], two complete runs. Eight gates unchanged: exact
chain repeats, finite physical output, zero-band null, zero-blood null, oxygen
direction, sample composition, grid composition and samples within measured grid.
No fitted bands or enlarged tolerance. Geometry/extinction/photon uncertainty
still excluded; the clinical calibration failure remains independent.


Curved direct transport evidence: reports/nirs_curved_direct_transport_l4.json.
Ratios1.09147217731317/1.23995852864639/1.29717758872134/1.32205330026006.
All eight wavelength/oxygen cases match predicted detector integer energy exactly
within unchanged3590 bound. Detected packet count can decrease from3590 because
rounded-zero weights remove packets; this is retained and does not imply fixed
live-path identity across different absorption settings. Each complete repeat at
fixed settings remains exact. No packet-energy residual or failed boundary event.

Conditional chain evidence: reports/nirs_curved_chain_v1.json. Fixed256 draws,
seed/bands unchanged; central95% interval[1.2151813391955226,1.2639812440041762].
Grid span0.0726743472168303<=0.07274482789472692. This conditional synthetic
interval excludes geometry, extinction, finite-photon and clinical uncertainty.
The OWN-GATE-FAIL calibration observer still has0 matching anchors. Neither
curved synthetic geometry nor these passing engineering gates establish measured
head/forearm anatomy or pulse-ratio clinical slope equivalence.

### Label-interface mechanism gates before execution

Before any general voxel-label mesh conversion, measure small exact synthetic
volumes: solid2x2x2 block, two adjacent layers, enclosed internal region, and two
voxels touching only on an edge. Count differing face-neighbor label pairs with
explicit exterior0 padding, per-region voxel counts and expected integer volume.
Two observation gates: complete tables repeat and every nonzero region has a
positive boundary-face count. No mesh/manifold or anatomical claim yet. Edge
contacts must remain explicit; no welding or tolerance choice before the table.


Measured before design:

| Fixture | Region voxel counts | Unique interface squares | Internal shared squares |
|---|---|---:|---:|
| block | 1:8 | 24 | 0 |
| layers | 1:4, 2:4 | 28 | 4 |
| cavity | 1:26, 2:1 | 60 | 6 |
| edge_contact | 1:2 | 12 | 0 |

### Shared label-face conversion preregistration

A separate label_interfaces_v1 module will emit every differing face-neighbor
interface once, using exact integer corner welding and positive-axis normals.
The label on the negative-axis side is back; the other is front. Region-outward
faces use back faces directly and front faces reversed. No smoothing or repair.
Before execution, gates are: full arrays and tables bit-identical in two runs;
triangle counts exactly twice the measured squares; signed regional volumes
exactly equal integer voxel counts at pitch1; every regional edge has two incident
triangles; paired directed regional edges cancel exactly. All four fixtures must
pass for an unrestricted manifold claim. Edge contact is deliberately included;
any failed gate is retained with its count, without changing the fixture.


### Internal transparent-region transport gates before execution

Use only the measured closed3x3x3 cavity control, internal label2 at[1,1,1].
Source[1.25,1.375,0.0001] avoids mesh diagonals, direction+z, initial label1.
Two fixed arms: every region transparent and non-scattering; or exterior tissue
absorption0.2 with internal label2 transparent. All indices1, reflection disabled.
Each arm uses two full1e6 runs. Ideal regional lengths1.9999 and1; the frozen
backend retracts1e-5 at each handled internal crossing (two crossings here), so
zero path-length error is not claimed. Fixed absolute regional error bound1e-4.
Gates: two complete outputs identical; exact packet energy; zero leak/cap/residual;
all photons exit at z3; both regions visited by every packet; maximum regional
path error<=1e-4; attenuation predicted from observed paths within1 integer packet
unit; transparent arm has exactly zero absorption. No statistical or anatomical
claim and no change to the existing boundary retraction.


Measured representation mechanism before the follow-up observer:

| Absorption used in prediction | Regional path length | Rounded escaped packet | Difference from transport |
|---|---:|---:|---:|
| double0.2 | 1.9998900000025261 | 719766504 | 5 |
| double(float32(0.2))=0.20000000298023224 | 1.9998900000025261 | 719766499 | 0 |

The transport interface explicitly casts optical properties to float32. A new
photon_internal_air_contract_probe keeps the same eight gates, including the
unchanged1 integer bound, using that represented input in the prediction. A ninth
gate requires all five output hashes for each arm to match the preceding run.
No transport change or tolerance adjustment; the original double-input failure
remains. Two complete1e6 runs per arm are required again before a passing claim.


The corrected observer reproduces every preceding transport output hash, with
zero leaks, residuals or caps. Transparent arm absorbs0 and escapes1073741824000000;
absorbing arm absorbs353975325000000 and escapes719766499000000. Every packet
visits both regions and exits at z3. Regional path errors are9.999997473864042e-6
and1.0000000000065512e-5, inside the original1e-4 bound. These nonzero retraction
biases and the earlier5>1 double-input prediction failure remain explicit.

### Tetrahedral interface mechanism preregistration

Measure before designing a tetrahedral material-boundary adapter: one tetrahedron,
two face-adjacent tetrahedra with equal and unequal labels, reversed input order,
duplicate overlapping cells, and three cells sharing one face. Record signed
six-volumes, unique/exterior faces, material interfaces, maximum face incidence
and opposite-vertex sidedness. Observation gates only: complete table repeat and
finite volumes. Invalid fixtures stay in the table; no geometry repair is assumed.


| Fixture | Signed six-volumes | Exterior faces | Material shared faces | Maximum incidence | Same-side pairs |
|---|---|---:|---:|---:|---:|
| single | 1 | 4 | 0 | 1 | 0 |
| same | 1,1 | 6 | 0 | 2 | 0 |
| layers | 1,1 | 6 | 1 | 2 | 0 |
| inverted | -1,-1 | 6 | 1 | 2 | 0 |
| duplicate | 1,1 | 0 | 4 | 2 | 4 |
| triple | 1,1,2 | 9 | 0 | 3 | 0 |

### Tetrahedral boundary adapter preregistration

A new tetra_interfaces_v1 accepts finite nodes, integral4-node cells and positive
integral material labels. Exterior is0. Emit one triangle per exterior or differing
material face, orienting from back to front by the opposite vertex, independently
of input cell ordering. Reject zero-volume cells, face incidence>2, and shared
faces whose opposite vertices lie on the same side. No geometric repair.
Six gates before execution: all arrays/tables repeat; expected triangle counts
4/6/7/7 for single/same/layers/inverted; regional signed six-volumes exactly match
the sum of absolute cell six-volumes; every region edge has2 incidents and directed
balance0; layers/inverted full output arrays identical; duplicate and triple
controls both rejected. This does not detect arbitrary interpenetration between
cells that do not share vertex indices, or certify every possible input mesh.


The adapter preserves input coordinates and emits each material boundary once.
It removes only same-material shared faces; no vertices are smoothed or moved.
These exact synthetic results do not certify external anatomical inputs or
intersection freedom between cells with unrelated vertex indices.

### Ocular transport mechanism preregistration

Start with an explicitly synthetic four-region optical stack (cornea, aqueous,
lens, vitreous) and a planar retinal scoring exit. This is a transport/refraction
control, not a measured eye or clinical retinal-dose model. A published schematic
model motivates the region ordering, not the chosen numerical parameters:
https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0109373 .
Published model tables differ in conic/index conventions; no claim of reproducing
one of those prescriptions is made. Geometry/optical constants below are the
complete synthetic definition; wavelength dependence must be explicit, not inferred.

Before a new mesh/transport adapter, measure three polar tessellations of four
curved interfaces plus retinal plane: radial/angular8/32,16/64,32/128. Radius3,
vertex z=[0,0.6,3.6,7.6,23.6], curvatures=[1/8,1/7,1/10,-1/6,0], spherical conics0.
Report triangles, regional watertightness/volumes, maximum facet sag and normal
errors against analytic surfaces. Observation gates: complete table repeat,
finite quantities and closed consistently oriented positive-volume regions.
Geometry acceptance for subsequent design is fixed sag<=0.002 and normal angle
<=0.02 radians; failed resolutions remain in the table.

Independently compare float32 transport-style Snell arithmetic with the unchanged
external patent-lens tracer on flat interfaces at incidence0/10/30/50 degrees,
using the same represented indices1.38/1.34/1.42/1.34. Record directional error,
Fresnel reflection probability and critical-angle/refraction status. No Monte Carlo
or retinal irradiance claim before actual propagation. External tracer path is
provided by the caller; only its content hash is recorded, no reference edits.


| Radial/angular | Triangles | Maximum sag error | Maximum normal angle | Fixed geometry bounds |
|---|---:|---:|---:|---|
| 8/32 | 2656 | 0.010641837277220567 | 0.02322184445333638 | FAIL |
| 16/64 | 10432 | 0.002758670862621493 | 0.01088402745110873 | FAIL |
| 32/128 | 41344 | 0.0007019793901950067 | 0.005258810098375222 | PASS |

All15 transmitted plane-refraction controls agree with the unchanged external
tracer within maximum5.98900040671424e-8 per direction component; one50-degree
high-to-low-index control is total internal reflection in both. The16/64 sag
failure is not excused by its passing angular bound. At32/128 a finite facet-normal
bias remains, so no stochastic agreement with smooth optics follows automatically.

### Analytic ocular boundary candidate preregistration

Use a separate spherical-cap/cylindrical-wall intersection backend for this exact
synthetic prescription, preserving the existing photon propagation/energy logic
beside the immutable tissue-path baseline. This removes the measured facet bias
without claiming arbitrary anatomy support. Add a transparent launch region from
z=-1 to anterior cornea so incident photons cross that refracting surface. Region
order: outside0, launch air1, cornea2, aqueous3, lens4, vitreous5. The final plane
z=23.6 is an ideal terminating retinal receiver (no receiver back-reflection),
not an anatomical retina. Radius3 bounds the side wall. Internal retraction1e-5
is retained and its directional/position bias is not hidden.

Two explicit synthetic wavelength arms532 and650 use declared optical rows;
indices and absorption are illustrative, not inferred patient properties. First
verify a bounded diagnostic before1e6 packets per arm. Gates fixed beforehand:
integer packet energy exact and balance>=0.99; zero leaks/residuals/caps; two full
runs bit-identical; all finite nonnegative scoring; zero-reflection retinal rays
agree with the unchanged patent-lens tracer within2e-6 direction components and
3e-4 position; per-launch direct-transmission counts within6 binomial standard
deviations plus1 of Fresnel prediction. The stochastic test covers reflection
sampling, while the separate deterministic bounds cover represented arithmetic
and retained retraction. No mesh tolerance or baseline solver is changed.

Initial bounded ocular diagnostic:10000 photons per wavelength/leg,6/6 fixed
gates. At532/650, retinal packets9716/9727, zero leaks/residuals/caps, exact packet
energy. Maximum direct retinal-position error2.453277971950496e-7 and direction
error7.086032516312457e-8. Largest per-launch Fresnel count deviation/bound is
0.4382126500245677. Both complete wavelength outputs repeat exactly. Proceed to
1e6 per wavelength/leg with identical geometry, optical rows, seed and gates;
no new tolerance, disabled interface, or smaller accepted sample count.


| Wavelength | Retinal packets | Absorbed packet units | Escaped packet units | Retinal map overflow units | Max position error | Max direction error | Fresnel error/bound |
|---|---:|---:|---:|---:|---:|---:|---:|
| 532 | 972655 | 47501765173862 | 1026240058826138 | 37115983084 | 1.9985093928018283e-07 | 7.086032516312457e-08 | 0.39461335936823866 |
| 650 | 973485 | 28397003551856 | 1045344820448144 | 25521725326 | 2.453277971950496e-07 | 5.8208959319294706e-08 | 0.3880684328380319 |

Both wavelength arms repeat exactly over1e6 packets. Each balances
1073741824000000 integer input units. Some reflected paths land outside the fixed
[-0.25,0.25] retinal map; their nonzero energy remains in the explicit overflow
bucket, never dropped or hidden by changing the scoring region. This is a
29-launch synthetic pupil quadrature with illustrative spectral coefficients and
an ideal receiver. No anatomical eye, clinical retinal dose, diffraction or
volumetric-scattering validation follows from ballistic Fresnel agreement.

### Retinal irradiance normalization preregistration

Convert committed integer maps to normalized incident-power density by dividing
by photons*PACKET and exact bin area. Keep the overflow fraction separate. Three
gates before running: bin-area integral plus overflow reproduces total retinal
energy fraction within1e-12; two complete normalized maps byte-identical; invalid
negative counts/non-increasing edges/energy beyond launch total rejected. A zero
map control must integrate to zero. No source spectrum or physical power is
invented; each wavelength is independently normalized to unit incident power.


| Wavelength | Receiver fraction | Explicit overflow fraction | Integral residual |
|---|---:|---:|---:|
| 532 | 0.9284932341826605 | 0.00003456695292517543 | 0 |
| 650 | 0.9470811107766163 | 0.000023768958939239383 | 0 |

Reproduce normalization with `python probes/optics/ocular_receiver_probe.py`.
The helper returns owned density/area arrays and separate overflow/total fractions;
it does not write artifacts. Its script entry point runs the probe and writes
only the explicitly named reports. Input energy uses exact integer summation.
Item15's synthetic ballistic transport/receiver gates are complete; anatomical
prescriptions, volumetric scattering and clinical dose remain outside this evidence.

### Data-defined NIRS chain mechanism preregistration

Before designing a declarative runner, observe both unchanged path libraries and
chains twice. Record path shapes, full propagation hashes, ratio spans and seam
shapes at the nominal input. Fixed observer gates: complete leg repeat equality,
finite physical outputs and existing nulls, exact identity between manual
Hill-to-optical seam composition and frozen propagation. Retain any difference
from archived chain reports as a separate negative; no rebaselining silently.


| Frozen chain | Path shape | Nominal ratio | Joint span | Sum isolated spans | Propagation seam |
|---|---|---:|---:|---:|---|
| Box | 4119 x 2 | 1.215323496343131 | 0.055127841665889354 | 0.06987483101741154 | 1 x 11, exact |
| Curved | 3590 x 2 | 1.2399585366270993 | 0.05616480545854996 | 0.07217381811721912 | 1 x 11, exact |

### Declarative NIRS runner preregistration

The measured seam is scalar oxygen -> Hill saturation; scalar blood fraction
scale -> two regional fractions; detector paths + saturation + hemoglobin +
fractions -> optical response; response and original inputs ->11-column record.
Implement a separate ordered node graph with a fixed module registry and explicit
input/output seam references, with no executable code/imports in specifications.
Preserve frozen functions, random ordering, reductions and record layout. A
parameterized observer retains existing null/composition gates. The two existing
box/curved specifications use their unchanged bands. A third specification alone
halves all three input bands on curved paths; it is another conditional synthetic
scenario, not a new physiology law.

Fixed gates: complete records and observer result byte-identical to each frozen
chain; all three full repeated observations exact; all eight original gates for
each spec; third narrower bands reduce measured span; unknown modules, missing
seams, forward references and malformed bands rejected. No tolerance introduced
for reconstruction; no clinical calibration or general-purpose compiler claim.

Initial declarative candidate:5/5 registered gates. Box/curved full records and
serialized observer results are exactly equal to the unchanged chains; all three
specs pass8/8 original gates twice. Third joint span0.028074704181443977 versus
curved0.05616480545854996. Four malformed specifications rejected. This initial
probe checks reference existence, not incompatible seam types. Before delivery,
add static seam-type validation and deterministic standalone entrypoint emission;
require wrong-type seams rejected and two generated source byte strings identical,
with all prior numerical gates unchanged. No new numerical tolerance.


| Specification | Joint span | Sum isolated spans | Original gates | Generated child repeats |
|---|---:|---:|---|---|
| nirs_box.json | 0.055127841665889354 | 0.06987483101741154 | 8/8 | exact |
| nirs_curved.json | 0.05616480545854996 | 0.07217381811721912 | 8/8 | exact |
| nirs_narrow.json | 0.028074704181443977 | See complete report | 8/8 | exact |

Run `PYTHONPATH=src python probes/optics/nirs_spec_probe.py` with the two
committed path libraries. `Chain(spec).propagate(paths, inputs)` returns records;
`Chain(spec).observe(paths)` returns the conditional uncertainty/null report;
`generate(spec)` returns deterministic standalone Python source. These APIs write
nothing. The generated script takes one explicitly supplied NPZ path argument,
reads `paths_0` and prints its result; the probe writes the named reports and
uses temporary generated children. JSON specifications cannot introduce imports
or executable code. The fixed registry has four module types, not arbitrary
program synthesis. All original numerical functions and chains remain unchanged.

### Photon cross-backend mechanism preregistration

With the ocular transport source frozen at86dc97e, capture10000 packets per
wavelength, two full runs on each of three remote CUDA architectures. Reuse all
six original transport gates, exact inputs and external reference hash. Before
any arithmetic redesign, report each full array's differing element count and
maximum absolute difference, including integer ledgers. Cross-backend gate is
strict byte identity, no epsilon. Capture validity additionally requires exact
source/reference hashes and each backend's repeat/physical gates. A missing
architecture is missing evidence, not a passing row. This bounded matrix does
not replace the existing million-packet energy evidence or certify anatomical
transport; local hardware is outside this reservation and remains unmeasured.


| Pair | Changed elements | Maximum absolute difference | Two full repeats |
|---|---:|---:|---|
| L4 / A10 | 0 | 0 | exact on each |
| L4 / H100 | 0 | 0 | exact on each |

All six original transport gates pass per backend at10000 packets per wavelength
and leg. Source/reference hashes agree; the capture inventory names match the
requested architectures. Four auditor controls separately catch one-ULP change,
signed-zero byte change and shape mismatch, while identical arrays pass. Run
`python probes/optics/ocular_cross_backend_probe.py` to re-audit committed
captures. The complete table includes all arrays, not only retinal summaries.
Local hardware and the rest of the determinism-suite families are unmeasured by
this module; no "all backends" or nightly-scheduler completion claim.

### Scattering photon cross-backend preregistration

The preceding cross-backend table covers ballistic Fresnel transport only. Before
any scattering arithmetic changes, capture frozen tissue_photon_mc_v6 on the
unchanged60-cube coefficients at10000 photons per arm, reflection off/on, twice
per remote architecture. This is a determinism observer, not a rerun of the
million-packet MCX statistical acceptance. Fixed capture gates: exact integer
energy, zero leaks/residual/caps, complete within-backend repeats and nonnegative
finite scoring. Cross-backend gate: full bytes identical, zero tolerance. Report
every differing array with counts and maximum difference if it fails. No change
to the baseline, physical coefficients, retraction, window, RNG seed or math flags.


| Pair | Changed elements | Maximum absolute difference | Within-backend repeats |
|---|---:|---:|---|
| L4 / A10 | 0 | 0 | exact |
| L4 / H100 | 0 | 0 | exact |

Each reflection arm has10000 packets per leg. Measured event totals822834 with
reflection off and1354942 with reflection on; boundary hits10000 and20178.
Every integer ledger balances10737418240000 input units, with zero failed packets.
Run `python probes/optics/photon_cross_backend_probe.py` to re-audit committed
captures. The common runtime was Warp1.13.0, NumPy2.5.3, SciPy1.18.1,
Trimesh5.1.0 and CUDA12.9.1 container runtime. This establishes bounded synthetic
scattering determinism on these three architectures, not arbitrary inputs or
local hardware. The baseline source remains unchanged.

### Hundred-million packet scale preregistration

| Measured reflection-on control | Value | Scale implication |
|---|---:|---|
| Packets | 10000 | Next requested report scale100000000 |
| Events | 1354942 | Linear estimate13549420000; not an observed large-run count |
| Terminal storage bytes per packet | 32 | 3200000000 bytes for one large ledger |
| Packet energy units | 1073741824 | Total107374182400000000 fits signed64 |
| Counter representation | signed64 | No counter-width change needed |

Use the unchanged reflection-on cube and source, one100000000-packet launch per
leg, two legs, on remote H100 with a bounded600-second child. Keep two complete
ledger hashes (streamed, no quantization), both absorption arrays and all counts.
Do not tile identical seeds into repeated small batches. Gates fixed beforehand:
exact integer energy, zero leaks/residual/caps, full output byte hashes identical,
finite nonnegative arrays, original MCX total-absorption relative error<=0.02 and
5x5x5-binned fluence L1<=0.05 against both committed reflection-on reference fields.
No tolerance changes or anatomical claims. This closes only the specified large
synthetic cube report if it passes; a timeout/resource failure stays recorded.

The committed `reports/ocular_receiver.png` visualizes the normalized receiver
arrays on one shared logarithmic scale. White bins contain zero scored energy;
the explicitly reported overflow is outside the pictured fixed region. Discrete
spots reflect the29-position synthetic pupil quadrature, not diffraction or an
anatomical retinal irradiance estimate. The image introduces no new numerical
claim beyond the receiver table above.


| Quantity | First run | Second run | Fixed gate |
|---|---:|---:|---|
| Absorbed integer energy | 29252410329520630 | 29252410329520630 | Accounting exact |
| Escaped integer energy | 78121772070479370 | 78121772070479370 | Accounting exact |
| Residual / leaked / capped | 0 / 0 / 0 | 0 / 0 / 0 | All zero |
| Events | 13862141954 | 13862141954 | Full repeat exact |
| Boundary hits | 200981651 | 200981651 | Full repeat exact |
| Maximum MCX absorption error | 0.0029679609100792875 | 0.0029679609100792875 | <=0.02 |
| Maximum MCX binned fluence L1 | 0.006544806559352599 | 0.006544806559352599 | <=0.05 |

Each run launches107374182400000000 integer energy units. Complete absorption,
terminal and counter byte hashes are identical. Both absorption/counter arrays
are retained in `reports/photon_scale_arrays.npz`; the two3.2GB terminal ledgers
were hashed in full on the worker, not truncated, and are not retained locally.
The immutable source SHA and MCX reference archive SHA are in the report.
Run `PYTHONPATH=src python probes/optics/photon_scale_probe.py` on an explicitly
reserved capable CUDA worker; no local throughput or cross-backend100M claim.
Zero observed failures here is a measured sample result, not proof of zero failure
probability for arbitrary geometry. Original anatomical topology/calibration
failures remain unchanged and unpromoted.

### Reusable remote matrix runner preregistration

The two measured photon families each pass strict cross-backend4/4, but captures
were orchestrated by a separate session wrapper. Provide a public, bounded Modal
entrypoint that uploads only these frozen modules/probes and the explicitly
supplied hash-checked refraction oracle. It runs L4, A10 and H100 sequentially,
each family in separate children, then emits one full-array matrix certificate.
No scheduler or local CUDA context is created. Fixed runner gates: three actual
architectures present, pinned input/source identities, all original family gates,
complete same-backend repeated arrays exact, cross-backend full arrays exact.
A third-party caller supplies the output directory; refuse an existing directory
instead of overwriting archived evidence. Two complete emitted matrices must be
identical after excluding no fields (omit volatile times from the certificate).

Reusable matrix integrity gate, before acceptance: bind each packet-array hash to
its family report, require the exact nonempty array inventory and four distinct
arm/leg records per family, reject a modified packet array and a missing array.
Matrix generation must not pass vacuously on empty archives. This adds an artifact
integrity gate without changing any physical or byte-equality tolerance.

Initial public-runner startup failed before capture entry:0 complete captures and
0 numerical children. Repository-root lookup at module import assumed a nested
source tree, while the remote import used a flat module path. Failed application
was stopped; move this client-only lookup into the local entrypoint and set an
explicit startup bound. Original zero-sample evidence remains in
`reports/photon_suite/startup_failure.json`. Numerical kernels/probes and all
acceptance tolerances remain unchanged; retry uses a distinct output directory.


| Family | L4/A10 max difference | L4/H100 max difference | Full repeated outputs |
|---|---:|---:|---|
| Ocular ballistic/reflection | 0 | 0 | exact |
| Tissue scattering/reflection | 0 | 0 | exact |

The public runner's fresh captures reproduce the zero-difference result. Matrix
integrity requires exact nonempty inventory, each packet-array hash bound to its
report and all four distinct arm/leg records. The startup failure above remains
separate; it contributes zero samples to this successful matrix.

Invoke with an explicit oracle and a NEW output directory:
`modal run probes/optics/photon_modal_suite.py --lens-reference <oracle-file> --output-dir <new-directory>`.
The supplied oracle must match the frozen hash in the script. Only the four
listed public source files and that explicit oracle are uploaded. Captures run
sequentially on remote backends; the client never creates a CUDA context. The
runner writes each capture plus two complete deterministic matrix JSON files,
refuses an existing output directory, and installs no scheduler. The original
local hardware and remaining suite families are not certified by this command.
Run `python probes/optics/photon_suite_integrity_probe.py` for CPU-only artifact
controls. The optional Modal client package is required for these runner imports.

### Tetrahedral partition topology mechanism preregistration

The anatomical-input boundary audit remains failed on regional nonmanifold edges.
Before any new transport or repair, distinguish the tetrahedral volume partition
from its per-material surface topology. Measure whole-domain boundary edge counts,
face incidence, connected cell components, signed volumes and regional edge
incidence histograms twice. Synthetic control: four tetrahedra around a shared
edge fill a closed octahedron; alternating labels make two material boundaries
nonmanifold without changing the tetrahedral partition. Fixed observer gates:
complete table repeat, whole-domain edge incidence exactly2 and signed volume
agreement<=1e-10 on valid controls, alternating-label regional defect detected.
No vertex splitting, material relabeling, smoothing or acceptance-gate change.
A conforming volume-partition observation does not certify the existing
surface-based transport on a failed material boundary.


| Control | Cells | Outer faces | Cell components | Outer edges with2 incidents | Material edges with4 incidents | Volume error |
|---|---:|---:|---:|---:|---:|---:|
| One material | 4 | 8 | 1 | 12 | 0 | 0 |
| Alternating materials | 4 | 8 | 1 | 12 | 2 | 0 |

Both volumes are1.3333333333333333. No index splitting or smoothing hides the
four-incident material edges. A cell-neighbor transport representation would
require its own measured traversal/physical gates; the existing boundary
representation remains unaccepted on this alternating-label control.

### Cell-ray traversal mechanism preregistration

Before a neighbor-walk design, use an exhaustive barycentric interval oracle on
all four unchanged synthetic octahedron cells. Measure512 fixed-seed directions
from a strict cell-interior source: number of positive-length cell intervals,
minimum segment length, interval gap/overlap and analytic octahedron exit error.
Keep a source exactly on the shared edge as an explicit ambiguous-source control;
never repair it with an undisclosed displacement. Fixed oracle gates: full table
and array repeats exact; interval coverage/analytic exit error<=1e-12; edge source
belongs to more than one closed cell and is refused as a unique interior launch.
This measures straight-ray geometry only, no scattering/refraction certificate.


| Rays | Mean cell intervals | Maximum intervals | Minimum length | Maximum gap/overlap | Analytic exit error | Cells at ambiguous edge source |
|---:|---:|---:|---:|---:|---:|---:|
| 512 | 1.650390625 | 3 | 0.0039033028447270635 | 0 | 2.220446049250313e-16 | 4 |

### Separate cell-neighbor ray walk preregistration

Implement a prepared cell adjacency traversal beside the exhaustive oracle and
surface converter. Use strict interior source selection, global-ray barycentric
exit parameters, explicit neighbor lookup and first exterior exit. Reject exact
edge/vertex exit ties instead of choosing an arbitrary material normal. Reject
nonmanifold outer boundaries; material edge junctions may exist inside a valid
volume partition. Nonadjacent intersecting cells remain outside the conforming
input contract and are not silently repaired.

Five fixed gates: all512 cell sequences match exhaustive intervals and parameters
within1e-12; complete arrays repeat byte-identically; input-cell vertex inversion
retains cell sequences/parameters within1e-12; ambiguous source and exact edge-hit
rays rejected; duplicate cells and malformed rays rejected. No scattering,
refraction, clinical or invalid-anatomical-input acceptance follows from this
straight-ray traversal. Keep all prior material-boundary failures unchanged.


| Cell order | Rays per run | Sequence differences | Maximum oracle parameter error | Full repeats |
|---|---:|---:|---:|---|
| Original | 512 | 0 | 4.440892098500626e-16 | exact |
| Inverted | 512 | 0 | 4.440892098500626e-16 | exact |

Inversion preserves the numerical contract within the fixed1e-12 bound but does
not preserve every parameter byte; its distinct hash is retained. This is not an
orientation-independent byte certificate. Source-on-edge, exact edge-exit,
zero/NaN direction and duplicate-cell controls all raise explicit errors. The
prepared geometry is read-only, and `trace(origin,direction)` returns an owned
array of cell index, entry distance and exit distance through the first exterior
boundary. Run `PYTHONPATH=src python probes/optics/tetra_ray_walk_probe.py`.

### Structured cell-network observer preregistration

Extend the frozen neighbor walker to a6x6x6 grid with six conforming tetrahedra
per cube (1296 cells),512 unchanged seeded directions and a strict interior
source. No walker edits. Compare each full cell sequence and entry/exit parameters
against exhaustive all-cell intervals, and compare final exit against the analytic
box. Fixed gates: sequences exact/parameters<=1e-12, analytic exit<=1e-12,
complete two-run trace hashes exact, closed valid partition. Record explicit
refusal counts as failures rather than omitting rays or reducing the grid.

First larger-grid observer exited1 during JSON serialization after two numerical
observations;0 numerical reports were exported, so no result was accepted. A
NumPy boolean in a gate was not JSON serializable. Convert only report gate flags
to builtin booleans, retain the failure record, and repeat the full observer with
the walker and numerical bounds unchanged.


This aggregate-only external summary contains no raw geometry or launch landmarks.
The input also fails at the whole-domain boundary, so its rejection cannot be
explained solely by internal material labels. No components were removed and no
acceptance tolerance changed. The new closed-partition ray walker does not accept
this geometry merely because its summed volume is accurate.


| Cells | Rays | Total intervals | Mean / max intervals | Refusals | Max parameter error | Max analytic exit error |
|---:|---:|---:|---|---:|---:|---:|
| 1296 | 512 | 2735 | 5.341796875 / 23 | 0 | 4.440892098500626e-15 | 1.7763568394002505e-15 |

Both complete observations and all saved arrays repeat exactly. The initial
report-only serialization error remains separate; no rays were omitted from the
successful repeat. The original neighbor walker was not changed.

### Separate CUDA cell walk preregistration

The measured maximum23 intervals motivates an explicit32-slot default output on
this fixture, with an independently tested two-slot overflow control. Port only
the straight-ray neighbor traversal to a separate Warp module; retain prepared
CPU geometry and strict source validation. Gates fixed before execution:512-ray
cell sequences exactly match the frozen CPU walker, parameters differ<=1e-12,
all full GPU arrays repeat byte-identically, no ordinary-ray error status, and
explicit edge-hit/capacity controls report failure rather than silently truncating.
No scatter/refraction or speed claim; no local GPU execution. The original CPU
walker, exhaustive oracle, geometry converter and photon transport remain frozen.


| Case | Rays | Parameter difference vs CPU | Error status |
|---|---:|---:|---|
| Normal capacity32 | 512 | 0 | 0 for all |
| Capacity2 | 512 | Not an accepted full path | 275 capacity failures |
| Exact edge-hit | 1 | Not an accepted path | Ambiguity status2 |

All cells, intervals, counts and status arrays repeat byte-identically. Failed
rays retain an explicitly counted prefix; it must not be treated as complete.
`trace(prepared, origin, directions, capacity=32)` writes no report artifacts;
the Warp runtime manages its normal compilation cache. The script entry point
runs the fixed probe on an explicitly reserved CUDA worker. No speed claim or
invalid-anatomical-input acceptance is inferred from these geometry results.

### Cell attenuation mechanism preregistration

Before a GPU energy-deposition stage, measure sequential Beer-Lambert attenuation
on the complete frozen cell intervals for transparent, heterogeneous and strongly
absorbing synthetic material rows. Record optical-depth range, zero-weight rays,
integer absorbed/escaped totals and complete array hashes. Fixed observer gates:
exact integer accounting with PACKET=2^30 per ray, nonnegative scoring and full
repeat hashes. No scattering/refraction or dose calibration is included.


| Material arm | Optical-depth range | Zero escaped weights | Absorbed integer units | Escaped integer units |
|---|---|---:|---:|---:|
| Transparent | 0 to0 | 0 | 0 | 549755813888 |
| Heterogeneous | 0.0017052748995819277 to0.21728975839142928 | 0 | 17129928324 | 532625885564 |
| Strong | 1.7052748995819278 to130.06569789860583 | 107 | 527206061720 | 22549752168 |

Strong attenuation quantizes107 escaped packet weights to zero; retain this
measured censoring rather than claiming an uncensored reusable exit library.

### Separate cell attenuation CUDA preregistration

Add a separate scoring kernel over accepted complete cell traces, with float64
continuous attenuation and int64 per-cell dose/terminal energy. Exact accounting
uses unchanged PACKET=2^30. Five gates before execution: complete integer dose and
terminal arrays exactly equal the CPU sequential reference for all three arms;
all full repeats exact; global integer energy exact; transparent null is zero
dose/full escape; incomplete two-slot traces deposit no partial dose and retain
one full residual packet for each explicitly failed ray. No path completion,
scattering, Fresnel or clinical-dose claim is added by scoring known intervals.


| Arm | Absorbed units | Escaped units | Residual units | Full CPU comparison |
|---|---:|---:|---:|---|
| Transparent | 0 | 549755813888 | 0 | exact |
| Heterogeneous | 17129928324 | 532625885564 | 0 | exact |
| Strong | 527206061720 | 22549752168 | 0 | exact |
| Explicitly capped | 811221310 | 253665590978 | 295279001600 | exact |

Every row balances549755813888 input units in both runs. Strong attenuation still
has107 zero escaped weights. No failed path's prefix contributes dose. The API
`attenuate(prepared, traversal, mu)` validates complete interval ledgers, returns
owned cell-dose/terminal arrays, and performs no report writes; Warp's normal
runtime cache remains external to those outputs. This is ballistic attenuation
on supplied paths, not an anatomical photon solver or a clinical dose estimate.

### Packet-rounding boundary mechanism preregistration

The512-path agreement does not prove general CPU/GPU equality near integer
rounding boundaries. Freeze the attenuation kernel and construct2048 fixed-seed
half-integer target weights, then use represented lengths -log(target/PACKET) and
their immediate lower/upper float64 neighbors (6144 isolated cell intervals).
Before GPU comparison, measure the CPU continuous-weight offsets from those
half-integers, integer output weights and complete repeated arrays. Fixed mechanism
gates: full repeats, exact integer energy, maximum offset<=1e-6 packet units.
These are isolated scoring inputs, not a whole geometry/transport experiment.
Then require strict full CPU/GPU integer-array equality with no tolerated tick
error; also retain same-GPU full repeats and exact accounting if that strict
cross-reference gate fails. Do not change packet resolution or round-to-nearest
policy to force the gate.


| Intervals | Below half | Exactly half in CPU arithmetic | Above half | Maximum offset in packet units |
|---:|---:|---:|---:|---:|
| 6144 | 1327 | 3529 | 1288 | 1.1920928955078125e-7 |

Integer totals3234854203864 absorbed and3362215562792 escaped balance the
6597069766656 input units. This is the table before a GPU boundary comparison;
the subsequent GPU comparison is recorded below. No broad CPU/GPU integer-equivalence claim follows
from the three passing CPU observer gates. The CPU table predates the GPU comparison; all existing kernels remain unchanged.

### Packet-rounding GPU boundary result


| Reference | Absorbed units | Escaped units | Residual | Differing rays versus CPU |
|---|---:|---:|---:|---:|
| Frozen CPU | 3234854203864 | 3362215562792 | 0 | 0 |
| GPU, each of two runs | 3234854203893 | 3362215562763 | 0 | 125 |

Strict full-array reference equality FAILS; full GPU repetition and exact energy
PASS. Net absorbed difference is29 packet units; maximum per-ray difference is1.
The previously passing512-path attenuation fixture retains its bounded result,
but does not establish universal CPU/GPU integer equivalence. Input lengths,
packet resolution, arithmetic and rounding policy were not changed. No timing
claim. A portable rounding-boundary contract remains open.

### Public-atlas path-library mechanism, before oxygen/detector design

Existing public-source diagnostics supply the prerequisite table:10000 packets,
5403027133501 absorbed and5334391106499 escaped units, zero leaks/caps/residual;
regional visited counts10000/2345/1524/1463/1228/9. The voxel interfaces retain56450
summed regional nonmanifold incidences. The distinct tetra atlas also fails outer
closure (171 edges), so neither geometry is promoted by transport accounting.

Freeze the existing path observer and repeat10000 packets twice with original
optical properties and twice with only absorption set to zero. This explicit
proposal-distribution control measures censoring, surviving detector populations
and six-region path support before any oxygen model or placement optimisation.
Record full output hashes, zero escaped weights, regional path summaries and
source-distance shell populations with fixed edges0/5/10/20/40/80/infinity.
Numerical gates fixed before execution: complete full-array repeats within each
arm; exact terminal energy; zero leaks/residual/caps; finite nonnegative paths;
all zero-absorption packets exit with full weight. Failure remains failure.
These are mechanism gates, not the final anatomy/calibration gate. The published
oximetry-band gate remains unresolved with zero matching calibration anchors;
a static760/850 ratio is not assumed equivalent to pulse AC/DC calibration.


| Arm, each repeated twice | Absorbed | Escaped | Residual | Leaks | Caps | Positive exits |
|---|---:|---:|---:|---:|---:|---:|
| Original | 5403027133501 | 5334391106499 | 0 | 0 | 0 | 9339 |
| Zero absorption | 0 | 10533407293440 | 204010946560 | 0 | 190 | 9810 |

| Source-distance shell | Original exits | Zero-absorption exits |
|---|---:|---:|
| [0,5) | 4658 | 4663 |
| [5,10) | 2210 | 2216 |
| [10,20) | 1490 | 1518 |
| [20,40) | 673 | 776 |
| [40,80) | 287 | 520 |
| [80,infinity) | 21 | 117 |

All four ledgers exactly balance10737418240000 input units. However, the proposal
has1.9% residual energy: counting that as escaped would hide the190 capped paths.
Full repetition, energy accounting and finite nonnegative paths pass; zero failed
paths and an uncensored proposal fail. No event cap, geometry or tolerance changed.
Raw public-atlas captures are caller-owned outside this repository. The pinned
atlas hash is in the report; no raw anatomy or source coordinates are emitted.

### Atlas oxygen diagnostic preregistration

The table above motivates an explicit unresolved-tail bound, not a larger event
cap. For a zero-absorption path with nonnegative target absorption, its eventual
weight is bounded above by exp(-partial optical depth). Assign every capped path's
remaining upper weight to each detector when computing conservative signal
intervals. This bounds omitted continuation for these captured histories only;
it excludes finite-photon sampling, geometry and physiology uncertainty.

Use a separate six-region illustrative oxygen model: blood fractions
[.02,.01,0,.04,.04,0], background absorption[.002,.004,.0004,.004,.004,0], Hb15,
760/850 extinction[[586,1548.52],[1058,691.32]] from
https://omlc.org/spectra/hemoglobin/summary.html. These fractions/backgrounds are
explicit assumptions, not fitted or validated atlas physiology. Keep scattering,
index and geometry fixed. Saturations .7/.8/.9/1, fixed collection shell[5,40).
Gates: full repeated weight/results bytes; finite positive signal intervals;
strict monotonic separated ratio intervals; analytic conditional log-ratio slope
vs central finite differences at .7/.8/.9 with relative error<=1e-7; matching
published calibration anchor present. The last gate remains false with0 anchors.
No accepted clinical or pulse-oximetry inference is made by this diagnostic.


| Saturation | Static760/850 ratio | Maximum unfinished intensity contribution |
|---:|---:|---:|
| .7 | 1.0379590590055212 | 1.953894756173434e-17 |
| .8 | 1.1120436680165557 | 3.185647562082733e-16 |
| .9 | 1.1937021362653342 | 5.264593561061942e-15 |
| 1 | 1.2841763725533222 | 8.858594404568358e-14 |

Report: reports/atlas_oxygen_bounds.json. Maximum analytic slope relative error
4.675738719157381e-10<=1e-7. The continuation bound follows from nonnegative
absorption along each unfinished history; numerical endpoints use ordinary
float64 and are not outward-rounded machine interval certificates. They exclude
sampling, input topology and optical-property uncertainty. Geometry is unchanged,
the original190 proposal caps remain failures, and no clinical band is certified.
Item25 therefore remains unaccepted despite monotone conditional signal and exact
original transport energy/leak observations. Detector information can next be
examined only within this explicitly assumed optical model.

### Atlas detector-support mechanism preregistration

Before placement design, measure12 candidate centers at equal azimuths and target
radius20 relative to the source in the input xy frame. Snap each target to the
nearest escaped training-half point (even packet indices), with first-index ties;
this yields a represented surface exit, not a validated anatomical landmark.
Detector acceptance is Euclidean radius5 about that center. Report even/odd packet
support and pairwise overlap, without moving centers after inspecting held-out
odd indices. Both capture repeats must give exact centers/masks/counts; every
candidate must have at least20 detections in each half and no two acceptance masks
may overlap. These support gates precede Fisher-information optimisation. Raw
centers remain caller-owned outside this repository; public reports use candidate
IDs and source-relative target angles only. No graph code changes.


| Candidate ID | Training detections | Held-out detections |
|---|---:|---:|
| 0 | 32 | 29 |
| 1 | 31 | 25 |
| 2 | 23 | 24 |
| 3 | 32 | 44 |
| 4 | 30 | 34 |
| 5 | 37 | 29 |
| 6 | 30 | 21 |
| 7 | 42 | 31 |
| 8 | 55 | 53 |
| 9 | 86 | 73 |
| 10 | 35 | 30 |
| 11 | 19 | 18 |

Report reports/atlas_detector_support.json. All12 centers are distinct. The full
candidate grid is not a valid independent-detector ensemble; its failures remain.

### Constrained detector information preregistration

The measured overlap motivates evaluating only simultaneous triples with center
separation strictly greater than10, so radius5 acceptance regions cannot overlap.
Use candidates with training support>=20 only; never use held-out counts for
selection. The fixed uniform reference is IDs[0,4,8]. Keep all12 candidate results,
including excluded11, in the mechanism report. No radius/center or threshold changes.
For the illustrative oxygen model at saturation.8, compute the Poisson-whitened
Jacobian of intensities with columns saturation and one shared log-gain nuisance;
parameter scales are1 and1, equal launched exposure in both wavelengths. Score the
smallest singular value of the6x2 Jacobian for three detectors. Select maximum
training score with lexicographic ties, then evaluate frozen choice on odd-index
held-out photons. Noise/exposure are assumptions, not measured detector noise.
Gates: full repeat bytes; selected/reference physical separation>10; selected and
reference support>=20 in both halves; positive finite full-rank scores; held-out
sigma_min gain strictly>1 versus the fixed uniform reference. Report all scores
and failed gates without changing candidates after held-out inspection. This is
conditional design evidence, not anatomical or clinical acceptance and not a
certificate of Monte Carlo sampling uncertainty.


| Triple | Training sigma_min | Held-out sigma_min | Held-out detected packets |
|---|---:|---:|---|
| Fixed uniform0/4/8 | 4.564514342213444 | 4.4183192187962765 | 29/34/53 |
| Training-selected5/7/9 | 5.265807908541792 | 4.916942056439588 | 29/31/73 |

Selected target azimuths150/210/270 degrees; minimum separation across selected
and reference pairs16.493781393984754>10. Training gain1.15364034675993; held-out
gain1.1128535112452005. All122 training scores retained in
reports/atlas_detector_information.json. Independent central-difference check of
the12 held-out whitened Jacobians gives maximum relative error
7.642126777727886e-11<=1e-7 (reports/atlas_detector_jacobian_check.json).
No held-out samples choose centers, eligibility or the winning triple. The support
observer retains its1/3 full-grid failure; conditional feasible selection does not
make that grid pass. The two repeats use identical captured histories, with even
and odd packet indices separating training from evaluation; an independent-seed
replication and uncertainty estimate remain outstanding. Raw positions stay with
the caller's public-atlas capture. This is an optical/noise-assumption-dependent
information comparison, not a validated physical detector recommendation.

### Independent-seed detector replication preregistration

Prior table: fixed IDs5/7/9 versus uniform0/4/8; held-out sigma_min4.9169420564
versus4.4183192188, gain1.1128535112. Freeze all centers from the original even
packet indices, selected IDs, acceptance radius5, optical/noise assumptions and
parameter scales. A separate capture uses seed20260913 (original20260912),10000
zero-absorption packets twice, without changing the50000-event cap. No new center
selection or optimisation on the new photons. Retain leak/cap/residual counts.
Replication gates before execution: all new full transport arrays repeat exactly;
zero leaks and exact energy ledger; at least20 detections per selected/reference
center; positive finite full-rank information; independent-seed sigma_min gain>1.
These conditional information gates do not waive the separate uncensored-path,
anatomical topology or clinical calibration failures. Report capped-path residual
and unfinished target-weight upper contributions explicitly.

### Packet arithmetic precision mechanism preregistration

The frozen CPU/GPU table differs on125/6144 intervals, maximum1 unit. Neither
backend is thereby proven correctly rounded. Compare both against a separate
Decimal exp observer with exact float64-to-decimal input conversion at80 and112
digits, each run twice. Keep length and packet scale unchanged; form continuous
weight PACKET*exp(-length), round half-up to the integer ledger. Gates: complete
repeated integer arrays;80/112-digit integer equality; exact packet accounting;
all observed differences from the frozen CPU/GPU ledgers at most1. Report counts
for each backend and the nearest half-integer margin; this finite-precision
stability check is not a formal arbitrary-input correct-rounding proof.


| Reference on6144 represented intervals | Escaped integer units | Absorbed integer units | Differences from Decimal |
|---|---:|---:|---:|
| Decimal80 and112 digits, two runs each | 3362215561016 | 3234854205640 | 0 |
| Frozen CPU | 3362215562792 | 3234854203864 | 1776 |
| Frozen GPU | 3362215562763 | 3234854203893 | 1747 |

Nearest high-precision weight to a half-integer is2.0110271567820346e-12.
Report reports/packet_rounding_decimal.json and complete repeated integer arrays.
No exact half-tie occurs in this finite input set at these precisions. Stable
high precision does not make the existing CPU reference disappear: its strict
GPU-equality failure remains125 rays. CPU agreement and correct rounding of the
represented mathematical expression are distinct requirements. Arithmetic-site
localisation is next, before any replacement design.

### Packet arithmetic-site localisation preregistration

Keep the same6144 inputs. At112 digits, compute exact-input exp, then separately
observe: rounding that exp to binary64 before power-of-two packet scaling;
rounding the high-precision continuous weight directly to binary64; and using
the existing CPU exp. Apply the unchanged floor(weight+.5) to each. Compare full
arrays twice, require equal integer arrays for the two high-precision-to-binary64
routes, and report differences of each route from both Decimal integers and CPU.
These gates localise represented intermediate rounding versus exp approximation;
no CUDA instruction attribution follows without observing that instruction.


| Route | Continuous differences versus CPU | Integer differences versus CPU | Integer differences versus Decimal |
|---|---:|---:|---:|
| High-precision exp rounded to binary64, then scaled | 2 | 0 | 1776 |
| High-precision weight rounded to binary64 | 2 | 0 | 1776 |
| Original CPU exp route | 0 | 0 | 1776 |

Report reports/packet_rounding_site.json; full intermediate and integer arrays
repeat exactly. Replacing CPU exp with a correctly rounded binary64 exp would
not fix any of the1776 integer disagreements on this input set: intermediate
binary64 rounding followed by integer rounding already reproduces all of them.
The two CPU continuous differences do not alter final integers. The observed
GPU-versus-CPU125-case difference still requires a GPU intermediate capture;
no particular GPU instruction or compiler operation has been blamed.

### GPU packet-intermediate observer preregistration

The CPU table localises1776 high-precision integer differences to intermediate
binary64 rounding. Freeze all prior inputs/arrays/kernels. A separate GPU kernel
records exp(-length), PACKET*exp(-length), weight+.5 and floor-to-int64 for6144
isolated intervals. Gates: full intermediate/integer repeats; exact final equality
to the previously captured GPU terminal weights; exact power-of-two scaling of
recorded exp; every CPU/GPU integer disagreement has a differing continuous weight.
Record differing continuous-value counts and maximum positive-float64 ULP distance
without asserting an instruction-level cause or a universal exp error bound.


| New-seed triple | sigma_min | Detections |
|---|---:|---|
| Frozen uniform0/4/8 | 6.601549285674637 | 50/89/116 |
| Frozen selected5/7/9 | 7.139809384260563 | 67/81/141 |

Capture report reports/atlas_independent_capture.json; decision report
reports/atlas_detector_replication.json. Seed20260913 differs from20260912 used
for placement; two10000-packet runs, all five full output hashes exact. Escaped
10525891100672 plus residual211527139328 equals10737418240000; leaks0, caps197.
The uncensored-proposal failure remains1.97% residual. The illustrative target
unfinished intensity contribution is at most4.129998466834137e-16 per detector.
Independent-seed gain8.1535% is smaller than original held-out11.2854%; no new
selection, tolerance or center changes. Sampling uncertainty remains unmeasured.

### Detector sampling uncertainty preregistration

Measure uncertainty in the observed independent-seed8.1535% gain before calling
it robust. Freeze all data, centers, six detector IDs, noise model and matrix
scales. Resample complete packet histories jointly for all detectors/wavelengths
using2000 ordinary multinomial bootstrap draws of10000 histories, fixed seed
20260914, repeated twice. Report full gain arrays,2.5/50/97.5 percentiles and
fraction<=1. Gates: full gain-array identity; finite positive scores;2.5-percentile
gain>1. This last gate may fail and will not be replaced by the mean or median.
This is an empirical bootstrap diagnostic, not a coverage-certified interval or
anatomical/physiological uncertainty model. Paired resampling preserves shared
packet contributions across detector channels.


Bootstrap gain percentiles2.5/50/97.5 are0.9879559049466795 /
1.0848807501074456 /1.1863267678398044. Fraction not above1 is0.041.
The observed independent point gain8.1535% does not pass the preregistered
uncertainty gate. Identical zero-contribution histories were aggregated into one
multinomial category; all532 active histories retain their joint channel values.
No cutoff, selection, radius or confidence percentile was changed. Reports
atlas_detector_bootstrap.json/.npz retain both complete score/gain arrays.

GPU arithmetic report packet_gpu_intermediates.json and complete arrays show
power-of-two scaling exact. Every integer disagreement has a differing
continuous weight; the separate observer exactly reproduces all frozen GPU
terminal weights. This localises the disagreement before integer conversion in
this arithmetic observer. It does not prove a particular instruction/compiler
cause or a universal transcendental error bound. CPU intermediate-precision
negative1776 and original strict CPU/GPU failure125 remain separate.

### Larger fixed detector replication preregistration

The10000-packet bootstrap fails its lower-percentile gate. Freeze centers, chosen
IDs and every model/threshold again; use one preregistered larger independent
seed20260915 at100000 packets twice on a remote accelerator. Evaluate the same
fixed replication and2000-draw paired bootstrap observers without retuning.
Lower2.5-percentile gain must still exceed1; no optional within-run early stopping
or replacement by a median. Keep prior failed sample results. This sample-size
extension is not an anytime-valid coverage certificate. Record leaks/caps/residual
and all full-array repeat hashes; no local GPU or throughput claim.

### Pulsatile observable mechanism preregistration

Existing static atlas ratio changes1.037959->1.284176 over saturation.7->1,
but matching pulse calibration anchors remain0. DOI10.1117/1.JBO.28.11.115002,
Eq.1-3, defines pulse optical-density change as log(diastolic/systolic intensity)
and identifies wavelength-dependent pathlength as part of its saturation
relationship. This does not supply a calibration band for the present model.

Before inverse-model design, measure a separate illustrative1% relative blood
fraction pulse, leaving background, scattering, geometry and capture unchanged.
Use the original complete/partial proposal paths and shell[5,40), saturations
.5/.6/.7/.8/.9/1. Compute pulse optical-density ratios and effective
blood-fraction-weighted lengths from direct exponential reweighting. Freeze the
pathlength ratio at saturation.8 only, then apply the two-extinction algebraic
inverse at every saturation. Report absolute saturation bias and unclipped
estimates (even outside[0,1]). Gates: full arrays/results repeat; positive finite
pulse optical densities; algebraic inverse exact at the.8 anchor within1e-12;
constant-ratio absolute saturation error<=.01 over all six cases. No tolerance or
pulse size will be adjusted after observing the result. This is a measured model
approximation, not clinical validation or a reproduction of the source experiment.


| True saturation | Weighted path ratio | Frozen-ratio estimate | Absolute error |
|---|---:|---:|---:|
| .5 | 1.0760871925589113 | .6151436951312054 | .11514369513120537 |
| .6 | 1.0234240991782382 | .674741066794849 | .07474106679484904 |
| .7 | .9714173965965845 | .736134739872306 | .03613473987230609 |
| .8 | .919936017797717 | .8000000000000002 | 1.1102230246251565e-16 |
| .9 | .868839096622828 | .8671620675617504 | .03283793243824962 |
| 1 | .8179671622026652 | .9386602583842492 | .06133974161575084 |

Report atlas_pulse_mechanism.json/.npz. Positive pulse signals and the single
algebraic anchor pass, while the fixed ratio fails the unchanged.01 error bound.
This is a consequence measured in the assumed model, not clinical validation.

### Path-dependent pulse inverse preregistration

The table motivates a separate bracketed scalar inverse that recomputes pulse
optical density from supplied paths at each trial saturation, instead of freezing
the path ratio. Keep the old approximation untouched. Train only on the original
proposal's detected paths; generate targets from independent-seed detected paths.
Use fixed shell[5,40), blood/optical assumptions and1% modulation. Five held-out
saturations.55/.65/.75/.85/.95; solve bracket[.5,1] with48 bisection steps. No clamp
or refit when the observed ratio is outside the endpoint range: refuse explicitly.
Gates: full predictions/records exact twice; same-library self consistency<=1e-10;
independent-library absolute saturation error<=.01; invalid/unbracketed controls
refused. Bracketing finds a root of the supplied model, not uniqueness or clinical
truth. Keep frozen-ratio estimates beside the new ones on the same targets.

The completed100000-packet seed20260915 capture is also a second held-out input
for the pulse inverse. Before evaluating it, keep the original4510-path training
library, five saturation targets, fixed-ratio comparator and.01 error bound. No
model update from the larger capture. Report this separately from seed20260913.

### Larger fixed replication result

The unchanged detector observers pass5/5 replication and3/3 bootstrap gates at
100000 packets, seed20260915 on H100. Both full transport captures are exact;
this is not a cross-device equality claim because seed and count also changed.
Selected sigma_min22.49426613604595 versus uniform20.111309468281373 gives
1.1184883894071078 gain. Bootstrap2.5/50/97.5 percentiles are
1.0857909664142673 /1.1177191683701746 /1.1510359447434173;0/2000 gains<=1.
All2000 gains and both repeats retained in atlas_detector_large_bootstrap.json/.npz.
The10k failure remains above; neither percentile nor detector selection changed.
Capture atlas_large_capture.json records1926 capped paths, residual2068026753024,
escaped105306155646976, leaks0 and exact107374182400000 ledger. Capped proposal,
atlas topology and clinical-calibration failures are not promoted by these gates.


| Independent library | Detected histories | Maximum fixed-ratio error | Maximum new inverse error |
|---|---:|---:|---:|
| Seed20260913,10000 packets | 4540 | .09395694224518636 | .002750026448719689 |
| Seed20260915,100000 packets | 45211 | .09463937431060943 | .00019492034174484196 |

Original4510 training paths and all optical coefficients remain frozen. Reports
atlas_pulse_inverse.json and atlas_pulse_inverse_large.json include all five
held-out targets, both repeats and unclipped estimates. Selftest forward error
2.545741395465484e-13, inverse error1.1102230246251565e-14, four malformed/refusal
controls; two complete selftests saved in pulse_path_inverse_selftest.json.
Initial selftest serialization failed after numerical evaluation with0 exported
reports (NumPy boolean); report-scalar conversion only fixed it, and both full
selftests reran. Receipt pulse_inverse_selftest_serialization_failure.json remains.
The prepared inverse owns read-only arrays and writes no artifacts. It brackets
a root of its supplied pulse model; it does not establish uniqueness, anatomy,
calibration or robustness to wrong optical coefficients.

### Pulse inverse optical-assumption stress preregistration

The matched independent100k library gives maximum inverse error.00019492034,
while the constant-ratio comparator gives.09463937431. Next measure whether this
agreement depends on assumed coefficients. Freeze the inverse trained at Hb15,
1% pulse and original background. Independently generate targets on the100k
library with one changed assumption at a time: Hb12/18, relative pulse.005/.02,
background scale.8/1.2, plus matched control. These are explicit numerical stress
values, not clinically validated ranges. Same five saturations and error<=.01;
record unbracketed observations as refusals, never clamp. Gates before execution:
full records repeat; matched control<=.01; all stress cases accepted and<=.01.
A failed last gate names the unsupported assumption rather than relaxing it.


| Changed target assumption | Maximum saturation error |
|---|---:|
| Matched | .00019492034174484196 |
| Hb12 | .011751110172818935 |
| Hb18 | .009062870354318298 |
| Pulse.005 | .00024276107296139493 |
| Pulse.02 | .00021336956013739972 |
| Background scale.8 | .005536694032364897 |
| Background scale1.2 | .00506906893786252 |

No unbracketed targets in these35 cases. Full table/repeats in
atlas_pulse_assumption.json. No target coefficient was fitted into the inverse;
the matched-model result remains conditional on its assumed optical inputs.

### Unknown-hemoglobin detector mechanism preregistration

The Hb stress failure motivates measuring a third information parameter before
redesigning detector selection. Keep the existing six centers and triples5/7/9
versus0/4/8 fixed on the100k independent capture. Append log-Hb to saturation and
shared log-gain in the Poisson-whitened Jacobian; all parameter scales remain1.
Compare sigma_min within this same three-parameter model, not against the numeric
magnitude of a two-column matrix. Gates: full Jacobian/SVD outputs repeat; analytic
log-Hb column agrees with central finite differences at step1e-5 within relative
error1e-7; finite positive full-rank triple matrices; fixed-selected/reference
sigma_min gain>1. Record individual-detector dimensional nullity (2observations,
3parameters) and full triple singular values. No center selection using this data.

### Direct atlas spectral-transport preregistration

The unknown-Hb detector mechanism above is deferred without execution while
direct atlas transport is verified. Packet-rounding investigation is closed.
Existing table: original atlas10000 packets, zero leaks/residual/caps and exact
repeat; zero-absorption proposal190caps; static reweighted ratios1.037959->1.284176
under explicit six-region optical assumptions. This motivates direct spectral
transport against the same geometry/seed rather than additional reweighting alone.

Freeze mesh converter/path transport/source/scattering/index; set absorption only
to the declared six-region oxygen model at saturations.7/.8/.9/1 and760/850. Each
of8 arms launches10000 packets twice, seed20260912. Use the exact represented
float32 optical properties in the separate path-weight predictor, retaining the
original float64 oxygen diagnostic unchanged. Collection shell[5,40) is fixed.
Numerical gates: all complete direct arrays repeat per arm; exact energy ledger;
zero leaks/residual/caps; direct collection intensity vs zero-absorption path
prediction differs by at most1/PACKET plus the explicit capped-tail bound;
four direct760/850 ratios strictly increase with saturation. The additional
matching-published-calibration gate remains false with0 matching anchors.
Record all full hashes, per-arm integer energy and detector counts. No event cap,
geometry, optical coefficient or acceptance threshold changes after observation.
This tests the assumed spectral transport, not validated atlas physiology.

### Published pulse reference input/algebra preregistration

Public source DOI10.1184/R1/24530353.v1 contains seven small files. All downloaded
bytes match the repository MD5 metadata; source stays outside staging. Dataset
SHA25649ae7d8e41c5e96464575a7d5a578bc08e03933c7d11815e8a4e2fb56b6242c6
has an8x5x10x2 intensity array; extinction SHA256
5abb6f2906a45ec1dc77f6e20bd8893847197293564bc230d3dadcaf58b9920a.
It is a Monte Carlo reference, not the paper's clinical recording dataset.

Before fitting or importing a calibration, inventory full intensity/optical-depth
arrays and test algebraic consistency of the demo's conventional expression on
known homogeneous mixtures at saturations.2/.4/.6/.8/1 and equal pathlengths.
Use760/840 and the source extinction conversion2.3026/1000 unchanged. Compare the
literal demo expression separately with an independently rearranged two-mixture
equation. Gates: full arrays/table repeat; all source intensities/pulse densities
finite positive; independently derived inverse error<=1e-12; literal demo inverse
error<=1e-12. A failure of the latter concerns this expression/control only, not
the separate self-calibrated algorithm or the clinical study. Source hashes pin
which expression was read; no upstream files are edited or republished.


At true saturations.2/.4/.6/.8/1, the literal demo expression returns
-.17982038020214547 /.07410003258828915 /.32802044537872355 /
.581940858169158 /.8358612709595925. An independent rearrangement returns the
inputs with maximum error2.220446049250313e-16. The expression uses the oxygenated
second-wavelength coefficient in its numerator; solving the declared mixture
ratio uses the deoxygenated coefficient there. This finding is limited to the
pinned demo's conventional-expression control and does not test its separate
self-calibrated estimator. Neither upstream code nor source data were changed.
The supplied reference is simulated760/840-compatible data, not matching clinical
anchors for the atlas760/850 observable. Report pulse_reference_mechanism.json.

### Published spectral-grid mechanism preregistration

Next evaluate the separate normalized-path estimator's actual numerical surface
before designing an adapter. Inventory all five simulated saturations at each of
ten source-detector distances (50 cases), source grid0:.01:1, all eight wavelengths.
Record the literal measured-path expression's nonfinite entries at zero saturation,
its finite-grid minimizer, and the algebraic mixture-form continuous limit at zero.
Keep the published concentration/scattering/extinction conversion unchanged.
Gates: full objective arrays repeat; mixture/source expressions agree for all
positive-grid points within1e-12; literal objective finite on the whole grid;
source-demo distance30 has maximum saturation-grid index error<=1. Other distances
are reported in full and are not silently promoted. Do not replace the literal
zero entry in the baseline; the limiting expression is a separately named arm.

### Direct atlas spectral transport result


| Saturation | Direct760/850 intensity ratio |
|---|---|
| 0.7 | 1.037959069441 |
| 0.8 | 1.112043678170 |
| 0.9 | 1.193702147425 |
| 1.0 | 1.284176382642 |

Maximum direct-versus-path intensity difference3.63659e-12 is below the frozen
1/PACKET plus capped-tail bound (approximately9.31323e-10). The direct transport
and proposal library used different cloud backends; this is bounded intensity
agreement, not cross-device full-array identity. All direct repeat pairs are exact.
The optical assumptions and nonmanifold atlas findings remain unvalidated;
monotonicity does not establish a published clinical calibration band. Per-arm
integer ledgers and all full-array hashes are in atlas_spectral_transport.json.
No local accelerator, timing claim, or baseline change was involved.

### Published spectral-grid mechanism result


The full50x101 objectives for both literal and continuous-limit expressions repeat
bit-for-bit, including nonfinite entries. Finite-only literal minimizers equal
continuous-limit minimizers in all50 cases. This explicitly separate diagnostic
does not assert how the source runtime handles NaN in its minimum operation.

| Source distance | Maximum saturation-grid index error |
|---|---|
| 5 | 3 |
| 10 | 1 |
| 15 | 1 |
| 20 | 1 |
| 25 | 0 |
| 30 | 0 |
| 35 | 0 |
| 40 | 2 |
| 45 | 0 |
| 50 | 3 |

Each index step is0.01 saturation. The distance30 control satisfies the unchanged
one-index gate, while short/far-distance errors reach3indices and are retained.
Source concentration, scattering and coefficient conversion are unchanged. The
reference contains simulated inputs, so matching clinical anchors remain0.
Full objective arrays and hashes: pulse_spectral_grid_mechanism.npz/.json.
Next: item26 unknown-concentration information audit under the already frozen
three-parameter preregistration; no further packet-rounding investigation.

### Unknown-hemoglobin detector result

The deferred mechanism resumed only after the direct atlas spectral result.
No centers were reselected and all three parameter scales remain one.


| Frozen triple | Singular values, descending | Detected packet counts |
|---|---|---|
| Selected5/7/9 | 55.12012485037246,23.437362872947187,4.0713906612328445 | 645,845,1557 |
| Reference0/4/8 | 47.077579632607154,21.410809668911163,2.6960210488318763 | 552,655,1129 |

The log-Hb finite-difference relative error is1.2012340586258946e-10 against the
unchanged1e-7 gate. Both triple matrices have rank3; all six individual detectors
have nullity1 (two observations for three parameters). Both complete100k capture
legs produce identical Jacobians and singular values, archived in
reports/atlas_detector_hemoglobin.npz with hashes/rows in the matching JSON.
The1926 capped proposals remain excluded detector observations, and no clinical,
anatomical, sampling-uncertainty or absolute parameter-precision certificate is
implied. This three-parameter gain is not compared numerically to two-column
singular values. The fixed-center reference and all earlier baselines are unchanged.

Next own mechanism: measure paired sampling variability of this three-parameter
gain before redesigning placement. Reuse the established2000 paired multinomial
packet resamples and seed20260914, preserve joint channels and aggregate only
identically zero contributions. Gates before execution: full score arrays repeat,
all scores finite positive, empirical2.5th percentile gain>1. This is an empirical
resampling diagnostic, not a formal coverage certificate or detector reselection.

### Three-parameter paired sampling result


| Gain percentile | Value |
|---|---|
| 2.5 | 1.221659866545761 |
| 50 | 1.5088149501864867 |
| 97.5 | 1.8864274198621611 |

Zero of2000 gains are<=1;5261 of100000 histories contribute to the six fixed
channels. Identical zero histories alone are aggregated. The1926 capped paths
and optical assumptions remain limitations; these empirical percentiles are not
certified coverage. Reports atlas_detector_hemoglobin_bootstrap.json/.npz preserve
both complete score sets and resampling-count hashes. Frozen baselines unchanged.

### Fixed-detector optical operating-point preregistration

The local three-parameter gain1.51014795 and empirical lower1.22165987 motivate
an operating-point mechanism before redesign. Keep the same six centers and
selected/reference triples, unit parameter scales, fixed background and scattering.
Evaluate the Cartesian grid saturation(.7,.8,.9,1) and hemoglobin(12,15,18) on the
existing independent100k captures, twice. No reselection or new photons.
Gates fixed before execution: full Jacobians/singular values repeat; every triple
matrix finite and rank3; minimum selected/reference sigma_min gain>1 across all12
points. Record every gain and singular value, including failures. This finite
assumed-optics grid is not a continuous robustness or physiological certificate.

### Fixed-detector operating-grid result


| Saturation | Hemoglobin | Selected/reference gain |
|---|---|---|
| 0.7 | 12.0 | 1.475374333755 |
| 0.7 | 15.0 | 1.516513690953 |
| 0.7 | 18.0 | 1.556119547501 |
| 0.8 | 12.0 | 1.470818946889 |
| 0.8 | 15.0 | 1.510147950439 |
| 0.8 | 18.0 | 1.548159432483 |
| 0.9 | 12.0 | 1.465376529018 |
| 0.9 | 15.0 | 1.502164702815 |
| 0.9 | 18.0 | 1.537866300738 |
| 1.0 | 12.0 | 1.459397898735 |
| 1.0 | 15.0 | 1.492892978798 |
| 1.0 | 18.0 | 1.525496805898 |

No detector centers were fitted or selected using these observations. Full spectra
and exact repeated arrays are preserved in atlas_detector_operating_grid.json/.npz.
The finite grid does not bound between-grid behavior or remove capped histories,
assumed scattering/background, unknown anatomy or missing clinical anchors.
Next mechanism: measure sensitivity to independent detector gain parameters;
shared-gain information does not establish identifiability with per-detector gains.

### Independent detector-gain mechanism preregistration

Existing12-point shared-gain table has minimum gain1.4593978987347527. Before any
placement redesign, reuse its complete Jacobians; split the shared log-gain column
into three disjoint detector columns. Parameter order is saturation, log-Hb,
log-gain0, log-gain1, log-gain2, each scale1. Six measurements yield five columns.
Measure full singular spectra and saturation sensitivity after orthogonal removal
of all nuisance columns, for selected and reference triples at all12 points.
Gates fixed before execution: complete arrays repeat; all full matrices finite and
rank5; all selected/reference profiled saturation-sensitivity ratios>1; adding
nuisance columns cannot increase sensitivity beyond relative1e-12. Record both
shared and independent-gain sensitivities. Include a six-independent-channel-gain
control: its nuisance space spans all six measurements, so residual sensitivity
must be<=1e-12 times the original saturation-column norm. This diagnostic does
not certify clinical precision, continuous-grid robustness or calibration.

### Independent detector-gain result


| Grid index | Selected/reference profiled saturation sensitivity |
|---|---|
| 0 | 0.318846664516 |
| 1 | 0.368699956575 |
| 2 | 0.397145060457 |
| 3 | 0.309786739109 |
| 4 | 0.362576355750 |
| 5 | 0.392408587016 |
| 6 | 0.300512073015 |
| 7 | 0.355391371462 |
| 8 | 0.387892975244 |
| 9 | 0.292264693383 |
| 10 | 0.347012540420 |
| 11 | 0.383254267282 |

Indices follow the saturation-major/Hb-minor operating-grid table above. All
five-column matrices have full rank; adding independent gain nuisance directions
reduces sensitivity as required. The per-channel-gain control removes saturation
sensitivity entirely within the fixed relative1e-12 bound. Thus the earlier
shared-gain placement advantage cannot be promoted to independently calibrated
detectors. Baselines and thresholds are unchanged. Full arrays and numeric spectra
are preserved in atlas_detector_gain_nuisance.json/.npz.

### Nuisance-aware selection preregistration

The observed all12-point reversal motivates a separate selection objective.
Use ONLY original even-packet training observations and existing candidate centers,
training count>=20 and pairwise separation>10. Enumerate the same eligible triples;
maximize saturation sensitivity after profiling log-Hb and three detector gains
at saturation.8/Hb15. Lexicographic tie rule; no held-out reselection.
Then evaluate the chosen triple and fixed reference0/4/8 on both existing100k
independent captures. Gates before execution: full training scores and held-out
Jacobians repeat; positive finite full-rank matrices; original support/separation
constraints; held-out profiled sensitivity gain>1. Record IDs and gain even if
failure. This is conditional on wavelength-shared per-detector gains; independent
per-channel gains remain unidentifiable. No clinical or sampling certificate.

### Nuisance-aware selection result


| Quantity | Value |
|---|---|
| Feasible training triples | 122 |
| Selected training counts | 23,30,35 |
| Selected evaluation counts | 512,447,621 |
| Minimum selected separation | 32.327428387693104 |
| Selected profiled sensitivity | 0.06781794141710024 |
| Reference profiled sensitivity | 0.2220963980761792 |

All rank/support/separation gates pass; the gain>1 gate fails unchanged. Centers
and scores use only original even training packets. The evaluation library was
already inspected for earlier mechanisms, so this is a reused diagnostic set,
not an untouched confirmatory validation set. No placement is promoted. Complete
arrays and all122 training scores are in atlas_detector_nuisance_selection.npz,
with identities/rows in the matching JSON. The old shared-gain advantage remains
conditional on that different model. Next: quantify training score instability
before any further placement design or fresh validation run.

### Detector score stability and contrast mechanism preregistration

Prior training-only winner2/6/10 loses on reused independent100k data with gain
0.3053536302459018. Measure before choosing another design: retain all122 frozen
eligible triples and centers; compare original training and full independent
score rankings. Partition the independent100k capture by packet index modulo10,
without reshuffling or dropping histories. Record every triple score and each
block winner; winners are diagnostics, never promoted placements.
For both complete capture legs require exact full matrices/scores/rankings, finite
nonnegative scores, and equivalence of the five-column nuisance projection to an
independent three-by-two contrast reduction within relative1e-10 (absolute1e-12
floor). In the reduction, eliminate each detector gain with its normalized
orthogonal two-channel contrast, then remove the remaining log-Hb direction.
Also test the prior fixed2/6/10 gain>1 in every block; retain each failure.
Record full-data rank correlation, winner frequency, all fixed-placement block
gains, and saturation/log-Hb contrast angle. No data-derived threshold changes,
new selection proposal, untouched-validation claim or sampling certificate.

### Frozen candidate score stability result


| Quantity | Value |
|---|---|
| Training/full score rank correlation | -0.15003254896388551 |
| Selected contrast angle sine | 0.003695060996178306 |
| Reference contrast angle sine | 0.010197048452380678 |
| Maximum independent contrast difference | 8.79296635503124e-14 |

Block gains: 0.923143472581, 0.618688091260, 1.022627995301, 0.487928249084, 0.538315959035, 1.182331681553, 0.515059850970, 4.618774889530, 0.933000194217, 1.491235604477.
Full score/rank/Jacobian/contrast arrays repeat exactly. Five-column projection
and independent contrast reduction agree under the frozen bound. The very small
contrast angles show nearly parallel saturation and log-Hb directions in this
model. Rank instability is measured, not attributed exclusively to Monte Carlo
noise. Full-data winner1/4/8 is diagnostic only and is not promoted. Reports
atlas_detector_score_stability.json/.npz retain all122 scores per data partition.

### Third-wavelength information mechanism preregistration

Pinned public extinction source5abb6f2906a45ec1dc77f6e20bd8893847197293564bc230d3dadcaf58b9920a
contains760:(586,1548.52),800:(816,761.72),850:(1058,691.32). The endpoint
coefficients match the frozen two-wavelength model exactly. Near-collinear
saturation/Hb contrast sines0.003695/0.010197 motivate measuring the third channel
before proposing different hardware or placement.
Keep all paths, scattering, background, blood fractions and centers fixed; evaluate
760/800/850 against760/850 at saturation.8/Hb15 for fixed triples5/7/9,2/6/10,0/4/8.
Match total incident exposure: per-channel weight2/3 for three wavelengths versus1
for two. Profile unknown log-Hb plus one wavelength-shared gain per detector.
Run both complete100k captures. Gates fixed before execution: full Jacobian/SVD/
profile arrays repeat; all matrices finite and rank5; every three/two profiled
sensitivity ratio>1; duplicated760/850/760/850 at half exposure preserves two-band
information within relative1e-10 (absolute1e-12 floor). Record all scores and ratios.
This is a path-reweighted mechanism under wavelength-independent scattering, not
an experimentally validated spectral optimum or extra-photon performance claim.

### Equal-exposure third-wavelength result


| Fixed triple | Two-band sensitivity | Three-band sensitivity | Ratio |
|---|---|---|---|
| [5, 7, 9] | 0.080526902640 | 13.732236553804 | 170.529799404663 |
| [2, 6, 10] | 0.067817941417 | 10.760981006786 | 158.674545141418 |
| [0, 4, 8] | 0.222096398076 | 12.680887197867 | 57.096320821537 |

The half-exposure duplicated-channel control changes sensitivity by at most
4.2368886177257536e-14. These are local sensitivity ratios under the stated
Poisson/optical/calibration model, not throughput or observed measurement-accuracy
improvements. Wavelength-shared detector gain remains an assumption; arbitrary
per-channel gains still remove identifiability. No new photons were simulated.
All spectra/Jacobians/projections: atlas_detector_third_wavelength.json/.npz.

### Three-band training selection preregistration

The measured third-band sensitivity motivates a separate equal-exposure selector.
Keep original even training packets, same122eligible triples and centers, training
support>=20 and pairwise separation>10. Use760/800/850 at per-channel exposure2/3,
saturation.8/Hb15, profile log-Hb and wavelength-shared per-detector gains.
Maximize profiled saturation sensitivity with lexicographic ties. Evaluate chosen
triple versus fixed0/4/8 on both existing100k independent capture legs. Gates:
complete scores/Jacobians repeat; finite positive full-rank matrices; support and
separation constraints; evaluation gain>1. Retain every failure; reused evaluation
is not untouched confirmation, and detector gain must not be changed by wavelength.

### Three-band selection result and fresh-seed preregistration


Training support30/42/86, evaluation655/845/1557, separation17.71093099033924;
all122 scores and complete matrices retained in atlas_detector_three_band_selection
JSON/NPZ. This reused-library result does not validate a spectral hardware choice.

Freeze centers and selected4/7/9 now; no further optimization before fresh seed
20260916,100000 zero-absorption paths twice using unchanged transport and50000event
cap. Preserve full repeat hashes, exact energy ledger, zero leaks; retain capped
counts/residual explicitly without requiring or claiming zero cap. On the new
capture use the same three-band exposure and nuisance model, reference0/4/8.
Acceptance gates before observation: full repeated inputs/Jacobians exact, finite
rank5 matrices, frozen support/separation, fresh-seed sensitivity gain>1. No
reselection or optional repeat until pass. Matching clinical anchors remain0.

### Analytic spectral identifiability control preregistration

The measured three-band gain57-171 could reflect removal of a two-band ambiguity.
Test this independently of atlas paths using homogeneous Beer-Lambert signals,
three distinct lengths10/20/30, unknown per-detector wavelength-shared log gains,
known background. Unknown oxygenated/deoxygenated concentrations and three gains
form a linear log-signal model. Use the pinned760/800/850 coefficients unchanged.
Before running: require three-band design rank5, two-band design rank4; exact
synthetic saturation/Hb reconstruction errors<=1e-10/1e-8 on saturation.7/.8/.9
and Hb12/15/18; full arrays repeat. Construct a nonzero two-band null perturbation
(delta oxygenated concentration.1, deoxygenated delta chosen from endpoint
extinction differences, detector gains compensated at760). Require two-band
log-signal difference<=1e-12 and third-band difference>1e-6. This is an analytic
identifiability control, not empirical anatomy or photon-noise validation.

### Homogeneous spectral control result


The explicit nonzero concentration/gain perturbation changes two-band log signals
by1.3877787807814457e-17 but third-band signals by0.021765924663144892. Full design,
forward signal, recovered coefficient and perturbation arrays repeat exactly in
spectral_identifiability_control.json/.npz. This verifies an analytic ambiguity,
not accuracy in noisy or clinical data.

### Reusable profiled-information design preregistration

The independent contrast observer agrees with fixed full-rank QR projections to
8.79296635503124e-14. Existing experiments assume full-rank nuisance columns;
a reusable seam must also handle redundant nuisance columns without projecting
away extra information. Add a new module beside these frozen observers.
API receives already-whitened interest and nuisance Jacobians with identical
observation rows; no implicit noise, calibration or parameter scaling. Use SVD
of nuisance columns and the explicit numerical-rank cutoff machine-epsilon times
max(matrix shape) times its largest singular value. Project only resolved nuisance
directions. Return projected Jacobian, complete singular spectrum padded for
underdetermined interest columns, Gram matrix and explicit numerical ranks.
No artifact writes in the API, no caller-input mutation or inference of physiology.

Gates before execution: homogeneous known spectra and redundant-nuisance controls
agree within1e-12; full-observation nuisance removes information within1e-12;
rank-deficient interests report rank loss; caller arrays unchanged; six invalid
shape/nonfinite cases refuse; complete selftest arrays repeat. Integration gate:
all saved independent-gain atlas QR profiles agree in norm within relative1e-10
with absolute1e-12 floor. Original observers and rank thresholds stay unchanged.

### Reusable profiled-information result


Run `PYTHONPATH=src python probes/optics/profiled_information_probe.py
reports/atlas_detector_gain_nuisance.npz --output reports/profiled_information_probe.json`.
API `profile_information(interest, nuisance)` returns owned read-only projected,
Gram and singular-value arrays plus numerical ranks/cutoffs; inputs must already
be whitened and parameter-scaled. No implicit calibration or artifact writes.
Raw projected values are retained even when rank is numerically zero. The six
malformed-input controls refuse before decomposition. Numerical tests do not
establish interval rank or clinical uncertainty. Full integration outputs are
preserved in profiled_information_probe.json/.npz. Existing observers unchanged.

### Common spectral-gain nuisance preregistration

Third-band gains assume relative wavelength response known. Before interpreting
these gains, append two common unknown relative wavelength gains to the existing
three-detector gain model (760 reference, free800/850 responses), with log-Hb
unchanged. Use the saved three-band9x5 matrices for all three fixed triples;
construct each additional column from the corresponding per-channel square-root
signal. All scales remain1. Compare saturation sensitivity using the new reusable
projection. Gates before execution: full arrays repeat; finite full rank7;
sensitivity cannot increase beyond relative1e-12; prior5/7/9 profile gain versus
reference0/4/8 must remain>1. Record all sensitivity retention ratios even if the
placement gate fails. No spectral recalibration, additional photons or new layout.

### Common spectral-gain result


| Triple | Sensitivity before | After | Retained fraction |
|---|---|---|---|
| 5/7/9 | 13.732236553803837 | 2.161971800011409 | 0.15743770445117636 |
| 2/6/10 | 10.760981006785501 | 0.45964648331177704 | 0.04271418033559765 |
| 0/4/8 | 12.680887197867023 | 1.5153011452196774 | 0.11949488403891467 |

The relative placement comparison passes while absolute sensitivity loses most
of its magnitude. It does not remove the need for calibration. This observer
uses the earlier fixed triples; it does not reselect or reinterpret the separately
frozen4/7/9 confirmation. Full arrays: atlas_spectral_gain_nuisance.json/.npz.

### Fresh three-band paired-bootstrap preregistration

Before reading the fresh-seed placement result, freeze an additional sampling
diagnostic: selected4/7/9 versus0/4/8;2000paired multinomial packet resamples,
seed20260914, two complete repeats on the new20260916 library. Same three wavelengths,
equal exposure and per-detector wavelength-shared gains. Retain all joint channel
contributions and aggregate only zero histories. Use algebraic per-detector gain
projection followed by log-Hb projection; check its unresampled score against the
frozen QR implementation within relative1e-10, absolute1e-12 floor. Gates: full
score arrays repeat; finite positive scores; algebraic agreement; empirical lower
2.5th percentile gain>1. No optional reselection or repetition until pass, no
formal bootstrap coverage or clinical statement. Preserve capped proposals.

### Frozen three-band fresh-seed confirmation


| Fresh-seed quantity | Value |
|---|---|
| Selected sensitivity | 13.73785641155489 |
| Reference sensitivity | 12.46319317345603 |
| Selected counts | 654,799,1539 |
| Reference counts | 571,654,1089 |
| Bootstrap2.5/50/97.5 percentiles | 1.0730823174431443,1.1022673420252738,1.1332502000647395 |
| Escaped packet energy | 105332999192576 |
| Residual packet energy | 2041183207424 |
| Capped histories | 1901 |
| Leaked energy | 0 |

Seed20260916 used100000 packets twice, unchanged transport/cap and zero absorption
proposal. Exact energy ledger includes residual; it is not an uncensored or99%
escaped-path acceptance. Every full transport-array hash repeats; receipt
atlas_confirm_capture.json binds inputs and raw-array hashes. Both confirmation
and bootstrap reports retain full matrix/score arrays. Algebraic bootstrap score
versus frozen QR differs at most1.0658141036401503e-14. Shared detector4 is preserved
in both triples by paired resampling, never treated as independent observations.
This frozen fresh-seed comparison supports the conditional three-band placement;
1901caps, optical assumptions, spectral response and missing clinical calibration
remain limitations. No event cap, tolerance, detector choice or baseline changed.

### Reusable grouped detector selection preregistration

Fresh fixed-placement gain1.10227421 supports extracting the measured selection
method into a caller-owned Jacobian API. Accept per-group whitened interest,
shared nuisance and group-local nuisance arrays, plus explicit feasible candidate
index tuples. Enumerate canonical candidate order; compose shared columns and
block-diagonal local columns, call the separate SVD profiler, maximize minimum
profiled singular value with lexicographic ties. No implicit geometry, noise,
calibration, feasibility relaxation, file writes or graph-repository dependency.
Refuse empty/duplicate/malformed candidates and experiments with no resolved
interest parameters. Gates before execution: known analytic selection/tie controls,
input nonmutation and invalid-input refusal, full arrays repeat; integration with
all122 stored three-band candidates keeps4/7/9 and matches frozen QR scores within
relative1e-10, absolute1e-12 floor. Confirm selected/reference fresh sensitivity
norms under the same bound. Preserve original examples as frozen baselines.

### Reusable grouped selector result


Fresh selected/reference sensitivity differences are1.7763568394002505e-15/0
against their frozen QR profiles. Full scores, ranks and fresh replay outputs are
in grouped_detector_selection_probe.json/.npz. Call `select_groups(interest,
shared_nuisance, local_nuisance, candidates)` with explicit feasible tuples and
already-whitened/scaled group Jacobians. Returned indices use canonical
lexicographic ties, scores/ranks are owned read-only arrays. Numerical rank loss
is scored zero; if no candidate resolves all interest columns, the API refuses.
This does not silently infer detector geometry, support, calibration or noise.

### Verification inventory mechanism and preregistration

Initial table scan found362 status declarations,227 VERIFIED-FRESH and4 unresolved
VERIFIED-FRESH paths. All four were historical pre-relocation optical paths with
explicit RENAMES mappings and separate current status rows. Six old-path table
entries (including two negative variants) are now historical prose retaining their
original statuses/numbers; no source or historical result changed. They must not
be treated as live verification promises or redirected to differently scoped gates.

There was no Makefile. Add a scoped `make verify-detectors` for two fresh CPU
recipes: profiled-information and grouped-selector integration, covering their
four current module rows. Bind each recipe to source, input artifact, decisive
RUNNING row and expected complete report hashes. Run in distinct temporary output
directories; preserve source artifacts. Gates before execution: both full fresh
reports match expected bytes and repeat; all internal gates pass; corrupted source,
row or evidence bindings refuse; full `make verify` refuses if any VERIFIED-FRESH
row lacks an explicit recipe. Report coverage, never claim whole-repo completion
from this slice. No CUDA context, remote job, new scheduler or private input needed.

### Explicit verifier result and legacy-cell preflight


Reports verification_detectors.json retain complete coverage and recipe hashes.
An explicit output report is marked incomplete before preflight, preventing stale
success after a later failure. Empty manifests and changed source/row/input/report
bindings refuse. This verification slice covers4current declarations, not225.

The legacy cell test checks subprocess exit0 and nonempty stdout. Read-only
preflight found the two aging cells write explicit required-gate dictionaries but
do not convert false required gates into nonzero process exit. Their original
required/bonus distinction remains authoritative; exit0 alone is insufficient.
Next explicit recipe observes both unchanged cells in isolated output directories,
twice. Gates before execution: complete output JSON bytes repeat; every declared
required gate and overall-required flag strictly true; valid nonempty boolean
gate schema; missing/false/string-valued required gates refuse in controls. Keep
bonus failures and the explicit abstention in the source record as such. Hash
complete records; do not republish citation prose or imply fresh literature review.

### Typed legacy-cell verification result


All three optional bonus gates pass in this replay; the separate deliberate
bonus-failure control is retained without changing required-gate semantics.
Existing abstention remains in the byte-identical full source output; no citation
or abstention is newly verified as medical evidence. Complete record hashes and
explicit counts are in aging_cells_probe.json. The recipe binds both original
cell sources and current status rows, and this observer's source/status row.

Current explicit coverage is7/226 fresh declarations in3recipes. Run
`make verify-detectors`, `make verify-cells` or `make verify-selected`; these are
scoped targets. Full `make verify` still exits2 with219unmapped declarations and
zero numerical children. Reports verification_selected_v2.json preserve the
complete current inventory and repeated recipe output hashes; the earlier4/225
coverage report remains historical evidence. No cell numerical source changed.

### Legacy cell schema inventory and numerical replay preregistration

Static AST inventory of182 current fresh cell sources finds4 with both the
explicit gates dictionary and required-overall marker,72 with only a gates key,
and106 with neither marker. These are source-schema counts, not execution results.
The two aging sources are already covered. The remaining glucose-meal and
four-state iron sources explicitly require ALL19/7 boolean gates respectively.
This differs from the aging required-prefix/bonus convention and gets a separate
observer; no permissive schema fallback is introduced.

Run each unchanged numerical source twice in isolated output directories on a
bounded remote CPU worker, with pinned source identities. Gates before execution:
all19/7 declared gates strictly true, required-overall flag true, exact gate count
and key inventory, full result-record bytes repeat, all child runs complete within
the fixed200-second child cap. False/missing/string-valued gates must refuse in
controls. Keep reported negative values if any gate fails; no source or numerical
tolerance changes. Raw citation prose stays out of new public reports; full record
hashes bind the replay. No fresh literature/clinical or GPU performance claim.

### Numerical required-cell replay result


The legacy glucose fresh status was based on process-level evidence and is now
corrected to OWN-GATE-FAIL. Source hash1c96a5970119e2ae2bc59b26f547b950ca425c35f102d476717de8ced97b1500
is unchanged. Full3069-byte output hashd76c278c39ab1d17cf82d1b23980885feec77004d529ace6522cd45b4fec85cb
repeats; the iron2470-byte output also repeats and its7/7status remains unchanged.
The failed glucose gates are five suppression points and two intervention timing/
endpoint controls. Next observer captures their computed values versus unchanged
recorded targets and5% bands, twice, without changing the numerical source.
Require full-record repeat and exact source identity; retain the original all-gates
failure rather than accepting exit0. Full gate vectors: required_cells_probe.json.
No clinical/literature claim follows from this software replay.

### Frozen glucose target mismatch values


| Control | Computed | Recorded target | Relative error |
|---|---|---|---|
| suppression_40 | 33.86057685641762 | 30.6 | 0.10655479922933403 |
| suppression_50 | 40.84520231861217 | 37.2 | 0.09798930964011211 |
| suppression_60 | 47.58136927051477 | 43.6 | 0.09131580895676075 |
| suppression_75 | 57.299892784142 | 52.8 | 0.08522524212390163 |
| suppression_100 | 72.72334174289064 | 67.4 | 0.07898133149689367 |
| endpoint | 230.23203113736403 | 214.7 | 0.0723429489397487 |
| crossing | 190.0 | 276.0 | 0.3115942028985507 |

Both full output hashes still equal the prior frozen replay. The five suppression
errors are7.898-10.655%, endpoint7.234%, crossing31.159%; all exceed the unchanged
5% source gates. This is a mismatch against the source's recorded targets, not an
independent clinical diagnosis. Source, target constants and integration settings
are unchanged. Report glucose_gate_values_probe.json; no retry-until-pass.

### Explicit all-gates cell adapter preregistration

Read-only inspection of the capillary and iron cells shows explicit all-gates
contracts:14 under `overall_pass` and7 under `required_gates_overall_pass`.
Unlike aging, neither excludes bonus gates. Add a pinned adapter for these two
contracts only; do not infer arbitrary legacy gate semantics from exit codes.
Run both complete source records twice with exact key inventories, strict boolean
values and original overall flags. Require all declared gates, exit0, full record
byte repeat and refusal of malformed/missing/false flags; keep any failure and its
values. Sources/thresholds unchanged, no claim of fresh literature verification.

### Explicit all-gates replay result


Both cells keep their original numerical sources/thresholds and source-defined
all-gates scopes. The source-key inventory is explicit, including different
overall flag names; no inference from process exit alone. Full output hashes and
gate vectors: explicit_all_gates_probe.json. Neither record's citation text is
newly validated. This separate passing recipe does not hide the glucose7/19
failure from the previous combined observer.

### Isolated CPU verification preregistration

Four explicit recipes now cover10/226 fresh declarations locally. Before calling
this slice reproducible in an isolated environment, ship only their pinned public
sources, inputs, expected reports, manifest and status document to a bounded remote
CPU worker. Run the complete selected verifier twice, each with its own fresh
numerical children and temporary outputs. Require both complete outer reports
byte-identical, all four expected report hashes unchanged, and all seven refusal
controls passing. Missing dependencies or environment-dependent numbers are
failures to retain, not grounds to weaken comparison. No private inputs, CUDA
context or implicit validation of the remaining216 declarations.

### Isolated CPU verification negative and dispatch preregistration

The isolated selected verifier repeats bit-for-bit but only1/4 complete expected
report hashes matches. All original child numerical gates pass. Fresh archive
verification_modal_0.json equals verification_modal_1.json; no thresholds were
changed. Full captured child reports locate profiler maximum difference
1.1379786002407855e-13 locally versus1.4151180227628402e-13 remotely; grouped maximum
score difference8.881784197001252e-16 versus5.329070518200751e-15. The cell recipe's
complete record hashes differ, while all14/7 booleans remain true. These are
cross-environment byte failures, not the glucose7/19 numerical failure.

Local inspection reports NumPy2.4.1, Python3.13.11, OpenBLAS0.3.30 dispatch Haswell,
NumPy baselineX86_V2/foundX86_V3. The remote image used Python3.12 with the same
NumPy/SciPy package versions; backend dispatch has not yet been isolated.
Before claiming a cause or changing any numerical source, run four process-local
remote arms: default, OPENBLAS_CORETYPE=Haswell only, disable NumPy X86_V4/
AVX512_ICL/AVX512_SPR only, and both. Same sources and complete expected reports;
two complete selected verifier runs per arm. Record actual CPU-library dispatch.
Gate: complete repeat per arm and all4expected report hashes exact for the combined
profile. Individual arms localize the effect; an unchanged host is not evidence of
a causal dispatch effect. No host/GPU setting, source tolerance or reference edit.

### CPU dispatch experiment result and limitation

All four process-local arms pass4/4 expected recipe report hashes twice on the
newly observed worker. Its DEFAULT runtime already reports the same Haswell/
X86_V3 dispatch as the local environment. Thus the original failure was NOT
reproduced in this experiment; the controls do not establish a repair or sole
cause. Python3.12.10 on this worker matches the local Python3.13.11 reports, so
the interpreter-version difference alone is insufficient to explain the earlier
mismatch. Full records: cpu_dispatch_mechanism.json; earlier failed isolated
verification_modal_0/1.json remain unchanged. No universal CPU portability claim.

The new default worker's entire selected-verifier report equals the local
verification_selected_v3.json. Current scoped evidence is4recipes/10declarations;
216of226fresh declarations remain unmapped. Full verification remains
OWN-GATE-FAIL and its earlier cross-environment byte failures are not erased.
Numerical sources, thresholds and expected reports remain unchanged. Optional
explicit artifact export now preserves each fresh child JSON/NPZ for diagnosis;
caller chooses a new directory, and an existing directory refuses overwrite.

### Pipeline-only cell replay preregistration

Read-only inspection identifies explicit all-gates scopes for geometric cardiac
output16 gates and minimal glucose-insulin9 gates. Both separate model plausibility
and unexecuted coupling gaps from their pipeline gates; those exclusions remain
unchanged and are not clinical validation. The output filenames differ, so a new
pinned adapter declares both exact paths rather than guessing a result filename.
Run each unchanged source twice in isolated remote CPU output directories.
Require exact gate key/count/type, source-declared overall flag, zero exit and
full output-record byte repeat; preserve any failure. Four malformed-gate controls
must refuse. No source parameter/threshold changes or private fixture inputs.

### Pipeline-only cell replay result


### Diagnostic relocation gates before edits

| Measured mechanism | Observation | Required relocation contract |
|---|---|---|
| Diagnostic discovery | Optical examples contain mutually importing observers and their launch helpers | Move the complete optical diagnostic family together; keep sibling imports valid |
| Cell observers | Five explicit verification observers use repository-relative roots | Preserve directory depth and complete numerical reports |
| Frozen numerical kernels | Script entry points refer to diagnostic paths | Change entry-point paths only; AST numerical definitions must remain identical |

Relocate diagnostics to probes, preserve report archives and all numerical gates.
Before committing require zero duplicate/dead current status paths, no forbidden
text in the requested source/test/document trees, identical numerical definitions
for moved kernels and engine entry-point edits, and two selected recipe runs with
all unchanged complete expected report hashes. Source and status pins may change
only to bind the explicit path edits; expected numerical report pins stay frozen.
Full verification must continue refusing unmapped declarations.

### Relocation verification result

All92diagnostic modules now live under probes. Main Status contains273engine
declarations, followed by92distinct diagnostic declarations; all365paths exist.
Frozen numerical definitions are unchanged, including moved failure kernels.
Two archived full cross-device reports replay exactly through the explicit
entry-point path translation; no GPU capture or archive rewrite was performed.
Selected recipes now cover13/227fresh declarations and all5expected complete
reports match. Full coverage still refuses214unmapped declarations. Earlier
cross-environment byte failures and the glucose7/19failure remain unchanged.

### Direct third-band transport preregistration

| Measured mechanism | Existing number | Next fixed comparison |
|---|---|---|
| Three-band fixed-placement fresh-seed gain | 1.1022742101770213; paired lower1.0730823174431443 | Direct800 transport against the same represented-coefficient path estimator |
| Existing direct760/850 transport | Eight arms, maximum intensity discrepancy3.63659e-12 | Same1/PACKET plus capped-proposal tail bound |
| Remaining spectral assumption | Scattering and refractive indices fixed across bands | No wavelength-dependent scattering or clinical claim |

Use the pinned800 extinction pair816/761.72, saturation0.8, concentration15,
float32represented regional absorption, unchanged source/mesh/transport/scattering.
Two10000-packet runs per seed20260912/20260913 on each H100/L4; match each seed
with its already captured zero-absorption proposal. Freeze five gates: full
per-arm repeat, exact integer energy, zero leaks/caps/residual, direct intensity
difference within1/2^30 plus unfinished-proposal bound, full cross-device arrays
identical. Preserve every failing gate; no kernel or numerical tolerance change.
All raw public-atlas captures stay outside this repository; ship aggregate numbers
and complete-array hashes only. Distinct result paths preserve earlier captures.

### Fresh-placement spectral nuisance preregistration

| Mechanism | Existing observation | Unresolved comparison |
|---|---|---|
| Independent detector gains | Frozen4/7/9 fresh-seed gain1.1022742101770213 | Unknown common relative800/850 response adds two nuisance columns |
| Earlier spectral nuisance on different placement | Sensitivity retention0.04271418 to0.1574377 | Does the currently selected placement still beat0/4/8 on fresh data? |

Use only the archived fresh-seed selected/reference9x5Jacobian pairs. Append
common relative800/850log-gain columns from the existing whitened signal, keeping
760response as the gauge. Fixed gates: complete repeated extended matrices and
profiles; all full matrices rank7 and finite; profiled sensitivity cannot increase
beyond1+1e-12; selected/reference gain strictly above1. No reselection, new
calibration datum, score scaling or changed reference. Save a loss with its number.

### Fresh spectral-nuisance result and sampling preregistration

| Fixed placement | Profile before | Profile after two common spectral gains | Retained fraction |
|---|---|---|---|
| Selected4/7/9 | 13.737856411554889 | 1.7136677017012683 | 0.12474054542162089 |
| Reference0/4/8 | 12.46319317345603 | 1.1635104559104081 | 0.09335572671604254 |

All4gates pass; gain1.4728425455879384. Loss of absolute sensitivity remains.
Next, freeze2000paired empirical multinomial draws, seed20260914, using the same
fresh100000packet library and fixed placements. Shared detector4 is resampled
jointly. Add both common spectral-gain nuisance columns before projection. Gates:
full repeated score arrays; finite positive scores; batchedSVD projection agrees
with independent frozenQR within max(1e-12,1e-10*score); lower2.5percentile gain>1.
Rank loss must refuse. Keep all1901capped paths in the original history count.
No calibration, continuous coverage or physiological uncertainty certificate.

### Spectral-nuisance sampling result

The two2000-draw arrays are byte-identical. Gain percentiles are
1.1343164615513959 / 1.4608122698298334 / 1.9368715620349748;
2/2000draws are not above1 and remain in the result. All4fixed gates pass.
BatchedSVD versus independentQR differs by at most2.6201263381153694e-14.
The100000histories include1901capped proposals; empirical resampling does not
establish coverage or correct unknown tissue properties. This does not restore
the absolute sensitivity lost to uncalibrated spectral response.

### Spectral observer recipe binding before replay

Two additional explicit detector recipes bind the existing common-response
observer and the fresh-placement observer to their source hashes, archived input
hashes, current status rows and complete expected reports. Each recipe gets two
fresh executions. Preserve all expected report hashes, finite/rank/monotonic/gain
gates and refusal controls. Coverage counts increase only for these two modules;
imported dependencies do not count as independently verified declarations.

### Strict endocrine output replay preregistration

| Source gate mechanism | Literal count | Required overall record key |
|---|---|---|
| hpg_male_axis | 16 | overall_pass_strict_all |
| ovarian_cycle_axis | 17 | overall_pass_strict_all |

Both source scripts call main without propagating gate failure to process exit.
A new typed adapter must therefore check every literal gate and the explicit
strict overall boolean, not infer success from exit0. Run each pinned unchanged
source twice on bounded remote CPU,200seconds per child. Require full record
byte identity, exact gate inventory/type/count, all required booleans, and four
malformed-record refusals. Existing biological/model exclusions remain unchanged;
this is numerical replay, not literature or clinical revalidation. Timeouts and
numerical failures remain explicit negatives, with no widened bound or tolerance.

### Direct third-band two-seed result

| Seed | Intensity on both backends | Path difference | Frozen bound |
|---|---|---|---|
| 20260912 | 0.19991624446930364 | 8.540390616929017e-14 | 9.31322721557791e-10 |
| 20260913 | 0.20018051357734948 | 7.605027718682322e-13 | 9.313227644531166e-10 |

All5gates pass. Eight10000packet launches, four per backend, have exact integer
energy and zero leaks/caps/residual. All40downloaded arrays (465864960bytes total)
match their recorded hashes and match fully within seed across repeats and
H100/L4, including shapes and dtypes. No throughput claim. Spectral scattering,
anatomy validity and clinical calibration remain unaccepted.

### Strict endocrine replay result and bounded recipe design

| Unchanged source | Strict gates | Complete bytes per run | Repeat |
|---|---|---|---|
| hpg_male_axis | 16/16 | 27796 | Exact |
| ovarian_cycle_axis | 17/17 | 136857 | Exact |

The adapter passes4/4, including four refusal controls. Its declared budget is
four sequential numerical children at200seconds each. Existing verifier recipes
retain their60second budget; this new recipe explicitly declares840seconds,
covering its predeclared child bounds and overhead. This is a resource limit,
not a numerical or performance acceptance gate. Reject noninteger, boolean,
nonpositive and above900second recipe limits before any child. Add an explicit
single-recipe selector so this larger replay does not rerun unrelated recipes.
Two fresh recipe executions must match the unchanged complete expected report.
Full coverage remains incomplete; imported dependencies do not count.

### Reusable numerical gate-record contract preregistration

| Measured mechanism | Observation | Required design |
|---|---|---|
| Existing typed adapters | Three assess functions have identical ASTs | One new library beside the frozen adapters |
| Known numerical failure | Glucose12/19required gates, process exit0 | A valid failure record must return accepted=false |
| Optional checks | Aging records distinguish required gates from bonuses | Explicit declared required subset; preserve bonus failures |

The new API accepts a complete declared gate inventory, its nonempty required
subset, and an explicit overall key. Require literal booleans and exact inventory;
reject missing/extra/type-invalid records, malformed declarations and inconsistent
overall booleans. Return an immutable verdict; do not mutate inputs or launch/write
anything implicitly. Before acceptance require known pass/fail/optional controls,
full repeat, invalid-input refusals, and parity with saved typed-adapter outcomes
including the seven-gate glucose failure. Do not replace frozen numerical sources
or existing adapters, broaden their validation scope, or infer required keys.

### Numerical gate-record contract result

The library passes6/6 controls, including seven malformed records and six
malformed declarations. The independent archived-vector observer passes4/4 on
16records per complete repeat, including the unchanged12/19glucose failure.
Required and optional failures remain distinguishable. Ordinal replay keys do
not revalidate the original source key inventory; that evidence stays with each
frozen typed adapter. No numerical cell or existing adapter was edited.

A new explicit recipe binds this library and observer to the four frozen adapter
reports and the complete expected parity report. Require two fresh executions,
all original controls and exact report hash; no added biological coverage.

### Declarative cell adapter preregistration

| Measured duplication | Existing evidence | Next seam |
|---|---|---|
| Three identical assessment functions | 16archived verdicts agree with the new contract | Explicit data recipe instead of another copied adapter |
| Pipeline cell outputs | Source hashes, gate-key hashes, counts and filenames already fixed | Declare only those fields; never guess required gates or output names |

Add a separate declarative all-required adapter for the existing16/9pipeline pair.
Require each exact source hash; a unique literal gate dictionary selected by the
explicit variable; exact count and existing key hash; explicit overall key and
output filename. Preserve the prior complete observer report format/hash, two
fresh full records per cell and four malformed-record controls. Refuse malformed
specs before children, paths outside the repository, duplicated cells, ambiguous
source dictionaries and changed pins. This first spec has all gates required;
optional-gate cells are outside its scope. No existing adapter or cell is changed.

### Declarative adapter and strict recipe results

The strict endocrine verifier ran twice on isolated CPU; both full outer reports
are byte-identical and its two fresh child report hashes match the original
expected report on each invocation. Other declarations are not certified by that single-recipe command; its report
names the selected ID.

The separate declarative pipeline adapter reproduces the exact original complete
pipeline_cells_probe.json bytes with all4gates passing. Its controls-only command
rejects10malformed specs, including ambiguous dictionaries, before any numerical
child. The explicit data spec contains only frozen source/key hashes, counts,
variable names, overall keys and filenames. The frozen bespoke adapter remains
unchanged; a new recipe binds the generic adapter and data spec to the same
expected numerical report, counting only the new adapter declaration.

### Next declarative cell mechanism and preregistration

| Inspected strict-overall schemas | Measured outcome | Mapping decision |
|---|---|---|
| Six additional sources | Four incompatible with one literal boolean dictionary: two dictionaries, dynamic nested gates, unpacked optional gates, embedded overall gate | Leave unmapped; no schema relaxation |
| melanin_photoprotection | One literal8gate dictionary and separate strict overall | Explicit frozen data recipe |
| nfkb_signaling_dynamics | One literal20gate dictionary and separate strict overall | Explicit frozen data recipe; optional coupling inputs absent |

Run the compatible pair twice on bounded remote CPU with the unchanged generic
adapter. Freeze200seconds per numerical child; timeout remains a failed result.
Require complete record repeat, all declared gates and overall values, zero exit,
and the existing four invalid-record controls. The absent optional coupling
inputs stay absent and their disclosed scope is not promoted. Do not infer any
new biological or clinical validity from numeric gate replay. Ten malformed-spec
controls must pass before launching numerical children.

### Additional declared-cell result

The data-only8/20gate recipe passes4/4; full24554/23369byte records repeat
exactly, and10malformed-spec controls refuse. Optional coupling inputs remain
absent. No new adapter implementation was needed. Bind only these two newly
executed cell declarations in the next explicit recipe; the generic adapter and
contract library are dependencies, not additional coverage. Preserve the current
complete expected report and four200second child bounds in its840second budget.

### Capped-history derivative mechanism preregistration

| Measured mechanism | Existing number | Unresolved effect |
|---|---|---|
| Fresh zero-absorption proposal | 1901capped histories of100000 | Omitted future detector signal and optical derivatives |
| Direct800 signal-tail accounting | Bound included in direct transport test | Derivative tails were not separately measured |

At the fixed0.8/15three-band point and frozen detector centers, retain all capped
histories. For accumulated optical depth t, future depth is at least t under
nonnegative absorption. Signal remainder is bounded by exp(-t); a parameter
with abs(dmu)<=k*mu has derivative remainder at most k*sup(u*exp(-u),u>=t),
where the supremum is1/e for t<=1 and t*exp(-t) otherwise. Optical trajectory
independence and fixed scattering remain assumptions. Compute the worst per-
detector remainder against every selected/reference signal and derivative.
Freeze gates: full repeated tables; pinned100000history source and1901caps;
finite positive observations; sampled analytic-envelope controls within1e-14;
all relative remainder estimates<=1e-6. This is a floating-point evaluation of
an analytic envelope, not an outward-rounded interval or clinical certificate.

### Capped-history mechanism result and profile-bound design

| Band | Minimum capped optical depth | Unnormalized signal remainder |
|---|---|---|
| 760 | 26.820835466983624 | 1.7488361970192683e-11 |
| 800 | 27.554067289433025 | 8.009137177757614e-12 |
| 850 | 32.495124975848775 | 4.3101844934653997e-14 |

All5mechanism gates pass. Worst relative signal/derivative remainder is
5.352299007440456e-12, with1901capped histories retained and both full tables exact.

Next add a separate reusable scalar profiled-information perturbation bound.
Given interest column s, full-column-rank nuisance B, and componentwise absolute
radii, use Frobenius norm bounds es/eB. If eB<sigma_min(B), a projector-distance
bound is eB/(sigma_min(B)-eB), capped at1; thus profile norm changes by no more
than es+projector_bound*norm(s). Refuse rank loss or an unresolved denominator.
Require analytic zero-radius control, exhaustive64corner perturbations of a
small independent projection example, immutable/nonmutating repeat and invalid
input refusals. Propagate the measured tail envelopes through signal square-root
normalization with stable rational differences, adding the measured difference
between reconstructed and frozen Jacobians to the radius. Require alignment
<=1e-10, full repeated bounds, and worst-case selected/reference gain>1. These
are evaluated perturbation formulas, not outward-rounded interval certificates.

### Capped-history profile-bound result

| Placement | Nominal profile norm | Evaluated lower | Evaluated upper |
|---|---|---|---|
| Selected4/7/9 | 1.7136677017012683 | 1.7136676934205193 | 1.7136677099820172 |
| Reference0/4/8 | 1.1635104559104081 | 1.1635104420615707 | 1.1635104697592455 |

The evaluated worst-case gain is1.4728425209401963>1. All4integration gates
pass; reconstructed/frozen Jacobian difference<=3.552713678800501e-15 is added
to the radius. The reusable bound passes5/5, including64independent perturbation
corners and6invalid-input refusals. All bound/radius arrays repeat exactly.
Nuisance uncertainty, fixed-optics assumptions and empirical sampling uncertainty
remain distinct; this does not provide an outward-rounded interval certificate.

Bind the bound library and observer to this complete expected report and its
archived tail/Jacobian inputs in one explicit verification recipe; two fresh
executions must keep every gate and report hash unchanged.

### Embedded-overall contract mechanism and preregistration

| Inspected mechanism | Measured schema | Required interpretation |
|---|---|---|
| deglutition_swallowing.compute_gates | 18literal subscript assignments, then one embedded overall boolean | Check all18 independently; do not count overall as a numerical gate |
| Output validity | Separate nan_inf_free_check | Require literal true in addition to the numerical gates |
| Script exit | main returns without gate-dependent exit | Exit0 alone cannot establish success |

Add a separate pinned observer; the existing flat-dictionary adapter stays
unchanged. Extract only literal assignments in the named function after source
hash validation. Require exact18key inventory plus the explicit embedded overall
key, literal booleans, agreement with all18required values, and the separate
finite-output check. Two bounded CPU runs,200seconds each, must produce complete
byte-identical records. Require seven malformed-record refusals and one valid
negative control. No source, tolerance, output data or biological scope changes.

### Nested gate-record repeat mechanism preregistration

| Source mechanism | Measured schema | Frozen check |
|---|---|---|
| spermatogenesis_sertoli | 21explicit gate objects with literal pass fields | Require all21booleans plus separate overall |
| Two output filenames | Results and evidence aliases | Require byte identity within each run |
| Generated metadata | _meta.generated_utc uses current time | Full cross-run bytes remain a gate; do not strip metadata to pass |

Run the pinned source twice with200seconds per child. Require all21nested pass
fields, strict overall consistency, source finite flag, both output aliases exact,
full record repeat and malformed-nested-record controls. Record changed JSON leaf
paths to diagnose failures. Any comparison excluding the known timestamp is
only explanatory evidence; it cannot satisfy the full-record repeat gate.
No clock patch, timestamp freeze, numerical-source edit or tolerance change.

### Embedded-overall replay result and recipe binding

The observer passes5/5. All18required booleans and the separate finite flag pass;
complete34770byte records are identical. Seven malformed controls refuse and the
valid negative remains unaccepted. Bind the source and observer to this complete
expected report in an explicit recipe,440seconds for the two200second children
plus overhead. Existing recipe limits remain unchanged. Two fresh recipe outputs
must match the complete report hash; no timestamp stripping or tolerance change.

### Nested full-record negative

The source passes21/21numerical gates, finite checking and both output aliases,
but full19928byte records differ. Exactly one JSON leaf changed:
/_meta/generated_utc. The observer is5/6 and the existing source declaration is
now OWN-GATE-FAIL for full-record repeat, retaining its numerical pass. No clock,
source, report field or acceptance gate was changed. This demotion reduces the
fresh-declaration denominator by one; it is not additional verification coverage.

## Source wording revision

A prose revision changes 204 source files. Identifier, operator and numeric tokens
and literal dictionary/index keys are unchanged. The prose-pattern search retains
only the executable `is_watertight` property and the existing descriptive coronary
gate key; both are preserved identifiers. Local framework/chain/photon tests pass
15/15; the separately checked coronary cell passes 1/1.

The graph test now creates its persistent lock before both complete byte snapshots.
This fixes a fresh-checkout setup failure without excluding any comparison. Two
clean temporary checkouts each pass three tests and skip the optional external
validator. The first cloud snapshot also omitted the existing JSONL ledger; that
transport omission is corrected in the complete sharded test input.

The initial sequential cloud run was stopped after 1393.53 seconds when the
corrected snapshot was already running. It recorded 166 passes, four graph setup
errors and two unchanged per-cell timeouts: arterial pressure sweep at 120 seconds
and gastrointestinal slow waves at 400 seconds. The arterial case passes locally
in 58.16 seconds under the original timeout. These failed observations are retained. The corrected full collection was then
executed in eight isolated cloud shards: 243 passed, three timed out and one skipped
(the external graph validator was not installed). All 247 collected test identifiers
were covered exactly once. Timeouts were the arterial sweep at 120 seconds,
gastrointestinal slow waves at 400 seconds and spinal CPG at 120 seconds. Local
isolated reruns pass all three in 58.16, 298.38 and 69.60 seconds respectively, with
the original bounds. The final fold-comment revision separately passes all nine
framework tests. `reports/wording_tests.json` records the environment versions and
counts. This is not a claim of a green full suite in a single environment.

Existing numerical receipts are retained as recorded. At that checkpoint, nine explicit recipe source
pins no longer matched after text changes: aging_cells, explicit_all_gates,
pipeline_cells, strict_endocrine, declared_pipeline, additional_cells,
embedded_overall, optional_coupling_measurement and coupled_activation. These
historical bindings were retained, and their archived full-output hashes do not
constitute a fresh verification of the edited text. Full `make verify` remains
unaccepted. The photon segmentation implementation and its measured source hashes
are unchanged by this revision.

Current source/receipt repair is documented in `docs/wording_verification_refresh.md`;
its new receipts supersede the earlier nine stale bindings without overwriting old reports.

## Diagnostic probes

These declarations retain their measured scope and negative results. Engine modules are listed in the main Status table.

| Module | Status | Evidence |
|---|---|---|
| `probes/optics/mcx_reference_probe.py` | OWN-GATE-FAIL | L4 external MCX0.7.1 cube60 finite/config gates PASS, whole-field repeat FAIL: maximum flux difference10224; absorbed packet energy identical177974.2482680604 of1e6. |
| `probes/optics/photon_cube_probe.py` | OWN-GATE-FAIL | Two/six gates PASS (full arrays identical, finite tally); four physical/reference gates FAIL with residual543971957541510 integer energy. |
| `probes/optics/hemoglobin_mechanism.py` | SYNTHETIC-ONLY | Three arithmetic gates PASS, two tables exact; dmu760/dS=-0.51541493109503 and dmu850/dS=0.19635160509280386 per mm at15g/dL. |
| `probes/optics/photon_failure_sites.py` | VERIFIED-FRESH | Two observation gates PASS;307890 grid exits and270136 segment-limit failures, zero region or ray-distance failures. |
| `probes/optics/photon_cube_v2_probe.py` | OWN-GATE-FAIL | Four/six unchanged gates PASS; one leaked packet retains480925106 integer energy, no cap packet. |
| `probes/optics/photon_last_failure.py` | VERIFIED-FRESH | Both observation gates PASS; packet431164 stops at zero sampled ray distance inside the cube, with no boundary hit. |
| `probes/optics/photon_cube_v3_probe.py` | VERIFIED-FRESH | Maximum MCX absorption error0.006667085187647551<=0.02 and binned-field L1 error0.01119298472586431<=0.05; two complete repeats exact. |
| `probes/optics/mcx_boundary_references.py` | OWN-GATE-FAIL | L4 external finite/diffusion gates PASS; exact external field-repeat gate FAIL for both cases. Maximum slab/diffusion bin error0.005022675449151793<=0.05. |
| `probes/optics/photon_boundary_probe.py` | OWN-GATE-FAIL | Two complete L4 runs per case; original cube 5/5, refractive cube 4/6, slab 6/8 gates PASS; failures retained above. |
| `probes/optics/tissue_photon_failure_probe.py` | CUDA-ONLY | Byte-identical relocation of the prior nonperturbing diagnostic; not rerun at the new path. |
| `probes/optics/tissue_photon_last_failure.py` | CUDA-ONLY | Byte-identical relocation of the prior nonperturbing diagnostic; not rerun at the new path. |
| `probes/optics/photon_arithmetic_mechanism.py` | SYNTHETIC-ONLY | Two exact CPU tables and finite-threshold gates PASS; photon_arithmetic_mechanism.json. Arithmetic observations, not yet transport attribution. |
| `probes/optics/tissue_photon_boundary_failure.py` | VERIFIED-FRESH | All frozen canonical absorption/terminal/counter hashes exact, both complete repeats identical for both cases; four observation gates PASS. |
| `probes/optics/photon_boundary_failure_probe.py` | VERIFIED-FRESH | photon_boundary_failure_l4.json; failure classification above, no transport change. |
| `probes/optics/photon_ray_replay.py` | SYNTHETIC-ONLY | Two identical CPU replay tables, all six double hits present; photon_ray_replay.json. |
| `probes/optics/photon_boundary_v4_probe.py` | OWN-GATE-FAIL | Two1e6 runs for each fixture, unchanged reference/profile/accounting gates. |
| `probes/optics/photon_boundary_v5_probe.py` | OWN-GATE-FAIL | photon_boundary_v5_l4.json; cube4/6, refractive2/6, slab6/8; no promotion. |
| `probes/optics/photon_vertex_accessor_probe.py` | SYNTHETIC-ONLY | Both CPU gates PASS; all36 extra-indirection corners wrong; photon_vertex_accessor.json. |
| `probes/optics/photon_boundary_v6_probe.py` | VERIFIED-FRESH | photon_boundary_v6_l4.json records all fixed reference, diffusion, accounting and repeat gates. |
| `probes/optics/layer_geometry_mechanism.py` | SYNTHETIC-ONLY | All3 CPU measurement gates PASS, two exact tables; layer_geometry_mechanism.json. |
| `probes/optics/mcx_layer_reference.py` | OWN-GATE-FAIL | Four/five external gates PASS; full-field repeat FAIL with max1.9669532775878906e-5; detector4217 and fraction0.00281081862729486 twice. |
| `probes/optics/photon_layer_probe.py` | VERIFIED-FRESH | All7 unchanged gates PASS; two1e6 full outputs identical; photon_layer_l4.json. |
| `probes/optics/photon_detector_probe.py` | VERIFIED-FRESH | Detector4119, fraction0.0027638959100153297; fixed external fraction/path5-SE gates PASS; energy reconstruction0. |
| `probes/optics/photon_path_reuse_mechanism.py` | SYNTHETIC-ONLY | All3 CPU bound/repeat gates PASS; photon_path_reuse.json. No baseline absorption truncation inside the fixed window. |
| `probes/optics/nirs_uncertainty_mechanism.py` | SYNTHETIC-ONLY | Three observation gates PASS; measured composition margin5.572968314870863e-5; Hb/fraction sensitivity difference0. |
| `probes/optics/nirs_direct_transport_probe.py` | VERIFIED-FRESH | All6 gates PASS,16 full1e6 transport runs; prediction integer error0 for all8 parameter pairs, full repeats exact. |
| `probes/optics/nirs_observable_contract.py` | OWN-GATE-FAIL | Three/four gates PASS, matching calibration anchors0 FAIL;12 rows/complete weights repeat exactly. |
| `probes/optics/curved_layer_mechanism.py` | SYNTHETIC-ONLY | Three/three observer gates PASS twice;32 segments fails fixed0.002 volume bound,64/128 pass. |
| `probes/optics/photon_curved_layer_probe.py` | SYNTHETIC-ONLY | Seven/seven gates PASS, two1e6 runs exact; leaks/caps/residual0, absorbed386616805627734 plus escaped687125018372266 equals launched integer energy. |
| `probes/optics/photon_curved_detector_probe.py` | SYNTHETIC-ONLY | Five/five gates PASS;3590 detected paths, fraction0.0021331821667263284, canonical/full hashes repeat; reconstructed energy error0. |
| `probes/optics/nirs_curved_uncertainty.py` | SYNTHETIC-ONLY | Three/three observer gates PASS;27-point grid repeated, joint span0.0726743472168303<=0.07274482789472692 isolated sum, Hb/fraction sensitivity difference0. |
| `probes/optics/nirs_curved_direct_probe.py` | SYNTHETIC-ONLY | Six/six gates PASS;16 full1e6 transports, all8 detector integer prediction errors0, exact full repeats, zero leak/cap/residual. |
| probes/optics/label_interface_mechanism.py | SYNTHETIC-ONLY | 2/2 observation gates; full table repeats; reports/label_interface_mechanism.json. |
| probes/optics/label_interface_mesh_probe.py | OWN-GATE-FAIL | reports/label_interface_mesh_probe.json; block48, layers56, cavity120, edge_contact24 triangles. No repair or unrestricted anatomy certificate. |
| probes/optics/photon_internal_air_probe.py | OWN-GATE-FAIL | 7/8 gates; absorption prediction error5 exceeds1, reports/photon_internal_air_l4.json; energy/leaks/repeats pass. |
| probes/optics/photon_internal_air_contract_probe.py | SYNTHETIC-ONLY | 9/9 gates, reports/photon_internal_air_contract_l4.json; 4 million packets total, both arms repeated exactly; represented-input prediction error0. |
| probes/optics/tetra_interface_mechanism.py | SYNTHETIC-ONLY | 2/2 observation gates; reports/tetra_interface_mechanism.json, exact repeated table. |
| probes/optics/tetra_interface_probe.py | SYNTHETIC-ONLY | Two full repeated tables/mesh hashes exact; volumes exact, regional edge defects0, duplicate/triple both rejected. |
| probes/optics/ocular_geometry_mechanism.py | SYNTHETIC-ONLY | 3/3 observer gates, full tables repeat, all regions closed; reports/ocular_geometry_mechanism.json. Geometry resolution failures retained below. |
| probes/optics/ocular_transport_probe.py | SYNTHETIC-ONLY | All six complete packet arrays and retinal histograms repeat; direct refraction and six-sigma Fresnel gates pass. |
| probes/optics/ocular_receiver_probe.py | SYNTHETIC-ONLY | Two independent transport maps per wavelength produce byte-identical normalized maps; zero-map control integrates to zero. |
| probes/optics/nirs_spec_mechanism.py | SYNTHETIC-ONLY | 3/3 mechanism gates; both complete repeated legs also exactly match archived evidence; reports/nirs_spec_mechanism.json. |
| probes/optics/nirs_spec_probe.py | SYNTHETIC-ONLY | Three specs each8/8 inherited gates twice; five invalid schema/seam controls rejected, third spec-only scenario reduces span. |
| probes/optics/ocular_cross_backend_probe.py | SYNTHETIC-ONLY | 4/4 strict gates across L4/A10/H100, all complete arrays identical; reports/ocular_cross_backend.json and reports/ocular_cross capture inputs. |
| probes/optics/photon_cross_capture.py | SYNTHETIC-ONLY | 4/4 capture gates per architecture; reflection off/on full scattering ledgers repeat exactly; reports/photon_cross. |
| probes/optics/photon_cross_backend_probe.py | SYNTHETIC-ONLY | 4/4 strict audit gates; complete L4/A10/H100 absorption, terminal and counter arrays identical; reports/photon_cross_backend.json. |
| probes/optics/photon_scale_probe.py | SYNTHETIC-ONLY | 6/6 unchanged gates at two100000000-packet launches; full output hashes exact, zero leaks/residual/caps; reports/photon_scale.json. |
| probes/optics/photon_modal_suite.py | SYNTHETIC-ONLY | Corrected public runner6/6, complete matrix bytes identical, two families x three remote architectures x two full repeats; reports/photon_suite_v1/matrix.json. |
| probes/optics/photon_suite_integrity_probe.py | SYNTHETIC-ONLY | 3/3 integrity-control gates twice; both valid captures accepted, all six modified/missing/duplicated controls rejected; reports/photon_suite_integrity.json. |
| probes/optics/tetra_partition_mechanism.py | SYNTHETIC-ONLY | 3/3 observer gates, full tables exact; closed four-cell controls retain zero volume error while alternating labels expose two material-edge defects; reports/tetra_partition_mechanism.json. |
| probes/optics/tetra_ray_mechanism.py | SYNTHETIC-ONLY | 3/3 exhaustive geometry gates,512 rays/full interval arrays repeat; reports/tetra_ray_mechanism.json. |
| probes/optics/tetra_ray_walk_probe.py | SYNTHETIC-ONLY | Two full repeats for original/inverted cells; all five ambiguous/invalid controls rejected. |
| probes/optics/tetra_ray_grid_probe.py | SYNTHETIC-ONLY | 4/4 on1296 cells and512 rays:2735 intervals, zero refusals/sequence errors, maximum4.440892098500626e-15 parameter error; reports/tetra_ray_grid.json. |
| probes/optics/tetra_ray_cuda_probe.py | SYNTHETIC-ONLY | Normal32-slot output has0 errors; two-slot control explicitly rejects275 rays with status4; exact edge-hit reports status2, controls repeat. |
| probes/optics/tetra_attenuation_mechanism.py | SYNTHETIC-ONLY | 3/3 seam observer gates, all full dose/terminal arrays repeat; reports/tetra_attenuation_mechanism.json. |
| probes/optics/tetra_attenuation_probe.py | SYNTHETIC-ONLY | Transparent, heterogeneous, strong and incomplete-path controls all pass;275 failed paths keep295279001600 residual units and deposit no partial dose. |
| probes/optics/packet_rounding_mechanism.py | SYNTHETIC-ONLY | 3/3 CPU mechanism gates,6144 represented intervals/full arrays repeat; reports/packet_rounding_mechanism.json. |
| probes/optics/packet_rounding_gpu_probe.py | OWN-GATE-FAIL | 2/3 frozen gates;125 of6144 rays differ from CPU by at most1 integer unit, both full GPU repeats exact; reports/packet_rounding_gpu.json. |
| probes/optics/atlas_path_mechanism.py | OWN-GATE-FAIL | 3/5 mechanism gates; both full-array repeats exact, original661 zero-weight terminations and proposal190 capped paths; reports/atlas_path_mechanism.json. |
| probes/optics/atlas_oxygen_bounds.py | OWN-GATE-FAIL | 4/5 gates;4510 selected exits,190 unfinished histories, four saturation ratios strictly separated and full weights exact twice; matching calibration anchors0. |
| probes/optics/atlas_detector_support.py | OWN-GATE-FAIL | 1/3 gates: repeated centers/masks exact, candidate11 only19/18 samples vs20, seven overlapping pairs with maximum10 shared packets. |
| probes/optics/atlas_detector_information.py | VERIFIED-FRESH | 5/5 conditional design gates; IDs5/7/9 selected from122 feasible triples, held-out sigma_min gain1.1128535112452005, full Jacobians/scores repeat exactly. |
| probes/optics/packet_rounding_decimal_mechanism.py | SYNTHETIC-ONLY | 4/4 precision-observer gates; all four80/112-digit integer arrays exact, CPU1776 and GPU1747 one-unit differences from stable Decimal result. |
| probes/optics/packet_rounding_site_mechanism.py | SYNTHETIC-ONLY | 3/3 gates; two binary64 intermediate routes exact and both reproduce all CPU integers, including1776 differences from high precision. |
| probes/optics/atlas_detector_replication.py | VERIFIED-FRESH | 5/5 conditional replication gates; independent-seed gain1.0815354207464527 at frozen centers/IDs, full transport/Jacobians repeat exactly. |
| probes/optics/atlas_detector_bootstrap.py | VERIFIED-FRESH | 3/3 at100000 packets, lower percentile1.0857909664142673; original10000-packet2/3 failure0.9879559049466795 retained. |
| probes/optics/packet_gpu_intermediates.py | SYNTHETIC-ONLY | 4/4 gates; frozen GPU integers exact,304 CPU/GPU continuous differences at maximum1 ULP produce125 integer differences; all arrays repeat. |
| probes/optics/atlas_pulse_mechanism.py | OWN-GATE-FAIL | 3/4 gates; constant path-ratio inverse maximum saturation error0.11514369513120537>.01; complete pulse weights repeat exactly. |
| probes/optics/atlas_pulse_inverse_probe.py | SYNTHETIC-ONLY | 4/4 on each independent library; largest error.002750026448719689 at10k and.00019492034174484196 at100k, both<=.01, full predictions repeat. |
| probes/optics/atlas_pulse_assumption_probe.py | OWN-GATE-FAIL | 2/3 gates;35 cases repeat exactly, Hb12 target against fixedHb15 model gives maximum error.011751110172818935>.01. |
| probes/optics/pulse_reference_mechanism.py | OWN-GATE-FAIL | 3/4 gates;800 positive intensities and400 positive pulse densities; literal conventional expression error up to.3798203802021455 on homogeneous controls. |
| probes/optics/atlas_spectral_transport_probe.py | OWN-GATE-FAIL | 5/6 gates; eight arms x10000 packets x2, complete arrays exact, energy100%, zero leaks/residual/caps; matching clinical anchors0. |
| probes/optics/pulse_spectral_grid_mechanism.py | OWN-GATE-FAIL | 3/4 gates;50 literal zero-grid NaNs retained, positive-grid expression difference1.3322676295501878e-15, source-demo distance30 index error0. |
| probes/optics/atlas_detector_hemoglobin.py | VERIFIED-FRESH | 4/4 conditional information gates, selected/reference minimum singular-value gain1.5101479504386228; full Jacobian/SVD repeats exact. |
| probes/optics/atlas_detector_hemoglobin_bootstrap.py | VERIFIED-FRESH | 3/3 conditional bootstrap gates;2000 paired draws twice, lower gain1.221659866545761, full score arrays exact. |
| probes/optics/atlas_detector_operating_grid.py | VERIFIED-FRESH | 3/3 conditional finite-grid gates; all12 gains>1, all triple matrices rank3, complete Jacobian/SVD repeats exact. |
| probes/optics/atlas_detector_gain_nuisance.py | OWN-GATE-FAIL | 4/5 gates; all12 profiled saturation gains below1, range0.29226469338334465-0.3971450604567481; complete repeated matrices/spectra/projections exact. |
| probes/optics/atlas_detector_nuisance_selection.py | OWN-GATE-FAIL | 3/4 gates; training-selected2/6/10 gains only0.3053536302459018 in the existing independent-seed library; all complete repeated scores/Jacobians exact. |
| probes/optics/atlas_detector_score_stability.py | OWN-GATE-FAIL | 3/4 gates; fixed placement loses6/10blocks,10different diagnostic winners; training winner ranks107/122 in full reused library. |
| probes/optics/atlas_detector_third_wavelength.py | VERIFIED-FRESH | 4/4 conditional information gates; three/two sensitivity ratios57.096-170.530 at equal total exposure; full arrays repeat exactly. |
| probes/optics/atlas_detector_three_band_selection.py | VERIFIED-FRESH | 4/4 conditional gates; training-selected4/7/9 yields sensitivity13.778019308030023 versus12.680887197867023, gain1.0865185608107564; complete repeats exact. |
| probes/optics/spectral_identifiability_control.py | SYNTHETIC-ONLY | 4/4 gates; two-band rank4 vs three-band rank5; maximum saturation error3.3306690738754696e-16/Hb error2.6645352591003757e-14 across9 noiseless controls. |
| probes/optics/profiled_information_probe.py | VERIFIED-FRESH | 3/3 integration gates,24 archived atlas matrices twice exact; largest norm difference from frozen QR1.1379786002407855e-13. |
| probes/optics/atlas_spectral_gain_nuisance.py | VERIFIED-FRESH | 4/4 conditional gates; full rank7 and exact repeated arrays, fixed5/7/9 gain1.4267604870700352 after profiling common spectral response. |
| probes/optics/atlas_detector_three_band_confirmation.py | VERIFIED-FRESH | 5/5 conditional gates; frozen4/7/9 fresh-seed gain1.1022742101770213, full transport/Jacobian repeats exact and receipt-bound. |
| probes/optics/atlas_three_band_bootstrap.py | VERIFIED-FRESH | 4/4 empirical gates;2000paired draws twice exact, lower gain1.0730823174431443, zero draws<=1. |
| probes/optics/grouped_detector_selection_probe.py | VERIFIED-FRESH | 4/4 replay gates; all122 scores repeat, same4/7/9 selection; maximum frozen-score difference8.881784197001252e-16. |
| probes/verification/aging_cells_probe.py | VERIFIED-FRESH | Current wording-bound replay 4/4 gates; two complete receipts exact, existing numerical gate values unchanged. reports/wording_refresh/aging_cells.json; archived scope in docs/wording_verification_refresh.md. |
| probes/verification/required_cells_probe.py | OWN-GATE-FAIL | 3/4 observer gates; glucose12/19 versus iron7/7required gates, both child pairs exit0 and complete output bytes repeat. |
| probes/verification/glucose_gate_values_probe.py | OWN-GATE-FAIL | 2/3 gates; seven original5% target checks fail identically, maximum relative error0.3115942028985507. |
| probes/verification/explicit_all_gates_probe.py | VERIFIED-FRESH | Current wording-bound replay 4/4 gates; two complete receipts exact, existing numerical gate values unchanged. reports/wording_refresh/explicit_all_gates.json; archived scope in docs/wording_verification_refresh.md. |
| probes/verification/pipeline_cells_probe.py | VERIFIED-FRESH | Current wording-bound replay 4/4 gates; two complete receipts exact, existing numerical gate values unchanged. reports/wording_refresh/pipeline_cells.json; archived scope in docs/wording_verification_refresh.md. |
| probes/optics/atlas_fresh_spectral_gain_probe.py | VERIFIED-FRESH | 4/4 on frozen fresh placement, gain1.4728425455879384; selected absolute sensitivity retention0.12474054542162089, full extended matrices repeat exactly. |
| probes/optics/atlas_spectral_nuisance_bootstrap.py | VERIFIED-FRESH | 4/4,2000paired draws twice exact; lower1.1343164615513959,2/2000gain<=1 retained; SVD/QR difference2.6201263381153694e-14. |
| probes/optics/atlas_third_band_transport_probe.py | VERIFIED-FRESH | 5/5;800band two seeds x two repeats on H100/L4,40full arrays exact within seed; maximum direct path difference7.605027718682322e-13 below9.313227644531166e-10. |
| probes/verification/strict_endocrine_probe.py | VERIFIED-FRESH | Current wording-bound replay 4/4 gates; two complete receipts exact, existing numerical gate values unchanged. reports/wording_refresh/strict_endocrine.json; archived scope in docs/wording_verification_refresh.md. |
| probes/verification/numerical_gate_contract_probe.py | VERIFIED-FRESH | 4/4,16archived gate vectors per repeat exactly agree with typed adapters;12/19negative preserved, no original-key or biological revalidation. |
| probes/verification/declared_cell_recipe.py | VERIFIED-FRESH | Current wording-bound replay 4/4 gates; two complete receipts exact, existing numerical gate values unchanged. reports/wording_refresh/declared_pipeline.json; archived scope in docs/wording_verification_refresh.md. |
| probes/optics/atlas_capped_history_mechanism.py | VERIFIED-FRESH | 5/5,1901capped of100000 retained; maximum relative signal/derivative remainder5.352299007440456e-12, full tables exact twice under fixed optics. |
| probes/optics/atlas_capped_profile_bound_probe.py | VERIFIED-FRESH | 4/4, evaluated gain lower1.4728425209401963 after capped-history radii; Jacobian alignment3.552713678800501e-15, full arrays repeat. |
| probes/verification/embedded_overall_probe.py | VERIFIED-FRESH | Current wording-bound replay 5/5 gates; two complete receipts exact, existing numerical gate values unchanged. reports/wording_refresh/embedded_overall.json; archived scope in docs/wording_verification_refresh.md. |
| probes/verification/nested_gate_repeat_probe.py | OWN-GATE-FAIL | 5/6;21/21nested gates pass, two19928byte records differ only at one timestamp leaf; aliases exact,7malformed controls pass. |
| probes/verification/optional_coupling_measurement.py | VERIFIED-FRESH | Current wording-bound replay 5/5 gates; two complete receipts exact, existing numerical gate values unchanged. reports/wording_refresh/optional_coupling_measurement.json; archived scope in docs/wording_verification_refresh.md. |
| probes/verification/coupling_schema_control.py | OWN-GATE-FAIL | 5/6; two synthetic arms falsely pass with0/3 and1/3comparisons; ten full source records repeat within each arm. |
| probes/verification/coupling_evidence_probe.py | SYNTHETIC-ONLY | 7/7; ten archived verdicts exact twice, four source-PASS legs correctly ABSTAIN; source numerical failure unchanged. |
| probes/verification/coupled_activation_measurement.py | VERIFIED-FRESH | Current wording-bound replay 4/4 gates; two complete receipts exact, existing numerical gate values unchanged. reports/wording_refresh/coupled_activation.json; archived scope in docs/wording_verification_refresh.md. |
| probes/verification/coupled_activation_contract.py | VERIFIED-FRESH | Current wording-bound replay 4/4 gates; two complete receipts exact, existing numerical gate values unchanged. reports/wording_refresh/coupled_activation.json; archived scope in docs/wording_verification_refresh.md. |
| probes/optics/photon_query_capture.py | VERIFIED-FRESH | 4/4prefix and complete capture gates;187882full-event queries, instrumented transport arrays unchanged, captures exact twice. |
| probes/optics/photon_segmented_probe.py | VERIFIED-FRESH | 4/4query gates;0/187882face mismatches twice,1.620-1.664ms vs2352.590-2395.178ms; input arrays unchanged. |
| probes/optics/photon_segmented_transport_probe.py | VERIFIED-FRESH | 5/5per fixed L4/H100sample; full5arrays exact vs baseline and repeats, zero leaks/caps, exact energy; no large-count or clinical claim. |
| probes/optics/tetra_ray_batch_probe.py | SYNTHETIC-ONLY | 6/6 gates;full record/offset bytes exact twice,oracle error4.440892098500626e-16<=1e-12,32/32 controls. |

| src/bodytwin/chains/nirs_reduced_chain_v1.py | SYNTHETIC-ONLY | 8/8 original gates and 5/5 comparison gates on two fixtures; complete records and arrays exact twice, 14 reduction controls. Modal CPU flat 0.262834/0.261987 to 0.210012/0.208811 s; curved 0.238953/0.239906 to 0.192976/0.193998 s. Extra transient 65904/57440 bytes. Unique-input batching has0bit mismatches but no strict flat speed win (in-place15.257/17.006 vs15.225/52.562ms); docs/nirs_exp_phase.md. Post-profile exp63.48/64.85ms;255 duplicate calls isolated to zero-band control47.934/48.123ms, no cache candidate or control bypass. reports/nirs_reduction.json, reports/nirs_duplicate_phase.json, docs/nirs_reduction.md. |
| src/bodytwin/geometry/optics/tissue_photon_launch64_v1.py | VERIFIED-FRESH | 5/5 unchanged gates at1024/8192 on L4/H100, full5arrays exact versus prepared baseline and repeats; zero leaks/caps. L4 API412.805/411.262 to322.675/321.052ms at1024 and417.088/416.776 to376.222/377.347ms at8192. H100235.466/252.957 to233.408/216.259ms and239.644/334.878 to224.120/224.843ms; variable timing and narrow2.058ms margin retained. docs/photon_launch.md; reports/photon_launch64_* JSON. |
| src/bodytwin/geometry/optics/tissue_photon_window2_v1.py | VERIFIED-FRESH | 5/5 unchanged transport gates all4count/device cases; full5arrays exact against launch64 and repeats,10boundary/tie/fallback controls CPU/GPU. At1024 L4319.578/318.315 to268.563/267.903ms,H100218.059/234.086 to164.843/181.271ms; at8192 L4374.043/374.105 to311.544/315.047ms,H100225.505/224.018 to184.077/183.920ms. Zero leaks/caps; initial upload failure0timed samples retained. reports/photon_window2_* JSON, docs/photon_window.md. |
| src/bodytwin/geometry/optics/tissue_photon_cell_v1.py | VERIFIED-FRESH | Certified integer-plane meshes only; Final focused CPU photon tests5/5, five sibling rows unique, sweep source pins current and targeted scrub0; reports/final_compute_check.json.5/5 original gates all4count/device cases, full5arrays exact against window2 and repeats,10CPU/GPUboundary controls plus unaligned refusal. At1024 L4251.698/250.170 to183.527/185.334ms,H100174.572/174.430 to127.669/117.086ms; at8192 L4298.359/299.624 to239.262/241.491ms,H100191.743/176.021 to121.820/115.402ms. Construction penalty64.337–241.597ms; no cold speed claim. Final sweep10000/100000/1000000 passes18/18 call checks, full reference/repeat arrays exact: L4 about38k/56k/60k photons/s,H10078k/334k/408k. L4 last-decade gain1.068x;H1001.221x, final H100 saturation unresolved;100million/s not reached. docs/photon_throughput_sweep.md, reports/photon_throughput_* JSON. docs/photon_cell_certificate.md, reports/photon_cell_* JSON. |
