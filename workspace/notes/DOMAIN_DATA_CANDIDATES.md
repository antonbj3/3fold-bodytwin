# Domain data candidates in bodytwin/data

**PENDING_INDEPENDENT_REVIEW**

Source (read only): `source_repository/data/`. Nothing was written there.

## Sammanfattning

164 directories traversed (of 505 entries total; 341 entries are top-level files) · 108 directories with numeric domain data · 13 already covered by a cell · 36 directories excluded as accounting/orchestration/non-domain · 15 directories without parsable numeric file (only mesh/video/binary) · 2 empty · 337 top-level files not opened.

### How the filter worked
1. All 164 directories were listed with file count and size. For each directory, the largest parsable JSON/JSONL-file under 8 MB (depth ≤ 2) was read and the number of numeric leaves was counted. 108 catalogs bar ≥ 10 numeric leaves under keyname with physical quantity; 3 bar < 10 (`cardiac_output_resting_spread_gate` 8, `msk_models` 9, `necroptosis` 5).
2. Excluded as accounting/orchestration/non-domain data data (36 directories): `acquired`, `anatomy`, `anchors`, `beslut`, `cache`, `claims`, `context_twin`, `camera_models`, `dataset_profiles`, `design_state`, `external`, `federations`, `fleet_metrics`, `fold_drafts`, `gpu_results`, `hpa_gtex_ingest`, `model_registry`, `oversikt`, `papers`, `prefetch`, `promotion_loop`, `provenance`, `raw_fetch`, `real_external`, `recept`, `ritningar_dragsystem_v1`, `ritningar_lindare_v1`, `ritningar_v1`, `scene_models`, `suffstat_cache`, `twin`, `twin_corpus`, `video_index_tables`, `video_segments`, `video_sources`, `vlm_dataset`. Base: content is run log, manifest, cache, provenance, dataset profiles, video index or CAD/drawing geometry — numerical but not physiological quantities. `ritningar_*` wears true measurements in mm but are workshop drawings, not domain data.
3. 15 directories had no parsable numeric JSON at depth ≤ 2 (mesh, npz, video or empty textures only): `bom`, `cad_meshes`, `cert_dataset`, `ensemble_legs`, `evidence`, `forensics`, `growth_legs`, `material_gt_ambientcg`, `material_leg_frames`, `polarization_leg_cache`, `prosth_passive_fit_shell_short_demo`, `repro_mjwarp_bug`, `surface_to_deep`, `video_frames`, `youtube`. The contents of them are **not opened** beyond filename and size — no claims about what they're missing.
4. 2 directories are empty: `altitude_acclimatization`, `aqueous_humor_iop` (0 files at depth ≤ 2). The latter is the only aqueous humor/IOP-the entry in the entire source that I found by name, and it's empty.
5. Of the 341 top-level files, 4 (`empir_photoled_spds.json`, `THROUGH_LENS_V2.json`, `BRDF_LIBRARY_V0.json`, plus a name grab over the small ones) was opened. 217 of them have names in UPPERCASE-form (`ARROW_LEDGER.jsonl`, `BUILD_QUEUE.jsonl`, `DISPATCH_SURFACE`, `DPO_PAIRS`, `CHEW_SCHEDULER_LOG.jsonl`, `*_manifest.json`, …) and 75 matches accounting/queue patterns directly. These are **not opened** and are treated as orchestration by name.

### Already covered by a cell (13)
`complement_cascade` → IMMUNITY/complement_cascade.py + complement_ap_bistability_eculizumab.py · `lymphatic_return` → IMMUNITY/lymphatic_return.py · `tcell_activation_exhaustion` → IMMUNITY/tcell_activation_exhaustion.py · `pd1exh_1784818729` → same cell (partial) · `autophagy_mitophagy` → MITOSTRESS/autophagy_mitophagy.py · `mitochondrial_oxphos` → MITOSTRESS/mitochondrial_oxphos.py + three oxphos variants · `blood_brain_barrier` → Q154 · `csf_davson_icp_flux` → Q156 · `capillary_starling` → Q020 (vasculature–interstitium–lymph) · `nephron_autoregulation` → Q019 · `renal_filtration` → Q140 · `glucose_meal` and `glucose_meal_trajectory` → SOLBENCH/native_glucose.py.

Deviation from the background task: `complement_cascade` and `capillary_starling` were specified as not covered, but `IMMUNITY/complement_cascade.py` and `Q020_model.py` are in the workspace cell list. The coverage is per domain, not per datafile — no cell references any path below `bodytwin/data/` in its `DAT_DATA_SOURCES.json`/`SOURCE_MANIFEST.json`, so the datafiles are still unlinked. `corneal_transparency`, `bcell_affinity_maturation` and `tissuetwin_recovery` have no corresponding cell.

## Top 25, ranked by tracks (a) → (e)

| # | Track | Path | Byte | Covered by Cell |
|---|---|---|---|---|
| 1 | (a) eye and optics | `tissue_lit_refs/measurements.nir_photobiomodulation_optics.jsonl` | 250358 | no |
| 2 | (a) eye and optics | `corneal_transparency/corneal_transparency_results.json` | 14467 | no |
| 3 | (a) eye and optics | `lens_accommodation/lens_accommodation_results.json` | 3570 | no |
| 4 | (a) eye and optics | `empir_photoled_spds.json` | 9300900 | no |
| 5 | (b) tissue behavior | `ipf_anchors/lung_fibroblast_stiffness_rawdata_2.json` | 244885 | no |
| 6 | (b) tissue behavior | `tissue_lit_refs/measurements.bending_vs_stretch_mode.jsonl` | 90042 | no |
| 7 | (b) tissue behavior | `tissuetwin_recovery/recovery_evidence.json` | 34258 | no |
| 8 | (b) tissue behavior | `erythrocyte_membrane/erythrocyte_membrane_results.json` | 16467 | no |
| 9 | (b) tissue behavior | `calcium_pth_vitd/calcium_pth_vitd_results.json` | 40448 | no |
| 10 | (b) tissue behavior | `msk_smoketest/moment_arm_validation/moment_arm_validation_evidence.json` | 283742 | no |
| 11 | (b) tissue behavior | `pulmonary_surfactant/pulmonary_surfactant_results.json` | 15749 | no |
| 12 | (c) immune system | `bcell_affinity_maturation/bcell_affinity_maturation_results.json` | 13568 | no |
| 13 | (c) immune system | `complement_cascade/complement_cascade_results.json` | 20592 | yes — Yes. IMMUNITY/complement_cascade.py and IMMUNITY/complement_a… |
| 14 | (c) immune system | `renal_c3g_joint_lesion_probability/results.json` | 3985 | no |
| 15 | (c) immune system | `hematopoiesis/hematopoiesis_results.json` | 35618 | no |
| 16 | (c) immune system | `iron_hepcidin/iron_hepcidin_results.json` | 27924 | no |
| 17 | (c) immune system | `erythropoiesis/erythropoiesis_results.json` | 17664 | no |
| 18 | (d) surgery and medical robots | `socket_contact_winkler/evidence.json` | 47667 | no |
| 19 | (d) surgery and medical robots | `socket_contact_edge_singularity/evidence.json` | 14095 | no |
| 20 | (d) surgery and medical robots | `prosth_real_task_grasp_sweep/evidence.json` | 24911 | no |
| 21 | (d) surgery and medical robots | `socket_two_patch_free_body_resolution/evidence.json` | 11388 | no |
| 22 | (d) surgery and medical robots | `prosth_socket_wholebody_thermal_transient_v1/trajectories.json` | 462826 | no |
| 23 | (e) other physiology | `deglutition_swallowing/deglutition_swallowing_results.json` | 34835 | no |
| 24 | (e) other physiology | `mucociliary_clearance/mucociliary_clearance_results.json` | 12809 | no |
| 25 | (e) other physiology | `fluid_compartments/fluid_compartments_results.json` | 25919 | no |

## Per kandidat: ordagranna tal ur filen

### 1. `source_repository/data/tissue_lit_refs/measurements.nir_photobiomodulation_optics.jsonl`
- Size: 250358 byte · Track: (a) eye and optics
- Covered by cell: no
- Numbers read from the file (verbatim, no rounding):
  - `mu_a (whole blood, HbO2 100%, 630 nm, Hct 0.45, C_Hb 150 g/L)` = 3.2665 cm^-1
  - `si (samma post, SI)` = 326.6458 1/m
  - `molar_extinction_HbO2` = 610.0 cm^-1/M
  - `molar_extinction_Hb` = 5148.8 cm^-1/M
  - `mu_s_reduced (whole blood)` = 270.0 cm^-1
  - `mu_eff (dermis, in vivo)` = 5.4748 cm^-1
  - `mu_s_reduced_powerlaw (skin)` = 48.9 cm^-1
  - `blood_volume_fraction (dermis)` = 0.2 %
- What a cell can calculate: Monte-Carlo or diffusion calculation of light transport in tissue: penetration depth, fluence versus depth and wavelength, and which wavelength maximizes dose in a given layer. The entries carry prop/unit/si/mods (wavelength_nm, SO2, hematocrit) plus DOI- or PMC-key, so a cell can be directly parameterized per wavelength.
- Note: Primary sources with verification fields (VERIFIED_PRIMARY_TABLE); The omlc table is specified to have 376 rows over 250–1000 nm.

### 2. `source_repository/data/corneal_transparency/corneal_transparency_results.json`
- Size: 14467 bytes · Tracks: (a) eye and optics
- Covered by cell: no
- Numbers read from the file (literally, no rounding):
  - `params.a_nm (fibril radius)` = 14.0 nm
  - `params.phi_areal` = 0.28 dimensionless
  - `params.d_hex_nm` = 50.3916571460783 nm
  - `params.n_fibril` = 1.411 brytningsindex
  - `params.n_matrix` = 1.365 brytningsindex
  - `params.dn` = 0.04600000000000004 brytningsindex
  - `params.L_nm (tjocklek)` = 500000.0 nm
  - `C_PREF` = 16.35109626093563 (enhet ej angiven i filen)
  - `validation.large_disorder_agreement.median_relerr_pct` = 7.9271096618596975 %
- What a cell can calculate: Transmission T(lambda) and attenuation tau_per_um for cornea at varying fibril order (sigma 0.00–0.06) over 400–700 nm. Can be converted to how much disorder is required for a given visible turbidity, and to edema sensitivity via phi_areal.
- Note: File has tau_per_um and T vectors per sigma level over 31 wavelengths.

### 3. `source_repository/data/lens_accommodation/lens_accommodation_results.json`
- Size: 3570 byte · Trace: (a) eye and optics
- Covered by cell: no
- Numbers read from the file (literally, no rounding):
  - `calibration.n_eq` = 1.4085022063402688 refractive index
  - `calibration.P_lens_relaxed` = 19.10999999999997 D
  - `accommodated_fixed_index.P` = 26.614173473168197 D
  - `accommodated_fixed_index.delta` = 7.504173473168226 D
  - `index_rise_required[0].target_amplitude_D` = 10.0 D
  - `index_rise_required[0].delta_n_relative_pct` = 0.49407233146678925 %
  - `stiffness_gradient_toy.average.k_per_year` = 0.09545699348069321 1/ar
  - `hofstetter_vs_measured[0].measured_D (age 9.24, n=5444)` = 14.44 D
- What a cell can calculate: Amplitude of accommodation in diopters from lens geometry (R1/R2) and refractive index; age-dependent loss via E/E14-curve, and how large an index gradient is required for a desired amplitude.
- Note: Also contains rigid_lens_falsifier.delta_D = 0.0 as null control.

### 4. `source_repository/data/empir_photoled_spds.json`
- Size: 9300900 byte · Trace: (a) eye and optics
- Covered by cell: no
- Numbers read from the file (literally, no rounding):
  - `wavelength_grid.step_nm` = 1 nm
  - `spds[0].meta.cct_min_duv_dataset` = 1705.35231381655 K
  - `spds[0].cct_agreement.diff_K` = 0.029950384083804238 K
  - `spds[1].meta.cct_min_duv_dataset` = 1713.38753774762 K
  - `spds[1].cct_agreement.diff_K` = 0.02993916675291075 K
- What a cell can calculate: Illumination side of an optical cell: spectral power distribution on 1 nm rasters as input to transmission and absorption in cornea, lens and skin, for example retinal dose per LED-spectrum.
- Note: Technical lighting and camera optics rather than ocular physiology; therefore ranks last in track (a). Only the first two spd records opened.

### 5. `source_repository/data/ipf_anchors/lung_fibroblast_stiffness_rawdata_2.json`
- Size: 244885 prey · Trace: (b) tissue behavior
- Covered by cell: no
- Numbers read from the file (literally, no rounding):
  - `fits.wi38_1|ACTA2|WT|excl_rigid|linear.E50` = 61.1344 kPa
  - `... .n (Hill-exponent)` = 1.782 dimensionslos
  - `... .r2` = -0.3993 dimensionslos
  - `... .N` = 61 datapunkter
  - `... .stiffness_levels_kPa` = [0.5, 32.0] kPa (forsta..sista av 4)
  - `... .leave_one_stiffness_out_E50_range_kPa` = [45.683, 61.134] kPa
  - `fits.wi38_1|ACTA2|WT|incl_rigid|linear.E50` = 4415823.4096 kPa
- `... incl_rigid stiffness_levels_kPa` = [0.5, 2300000.0] kPa (first..last of 5)
- What a cell can calculate: Dose-response curve for fibroblast mechanotransduction against substrate stiffness: E50, ie the stiffness where half ACTA2-the answer is reached, and the Hill exponent, including how much E50 moves when the stiff glass level is included or excluded. A cell can figure out what level of stiffness in a tissue model tips over to the fibrotic phenotype.
- Note: The directory has 60 files, including collagen_mass_ratio.json, fvc_decline.json and focus_fibroblast_density.json. Only this file and the file list opened.

### 6. `source_repository/data/tissue_lit_refs/measurements.bending_vs_stretch_mode.jsonl`
- Size: 90042 prey · Trace: (b) tissue behavior
- Covered by cell: no
- Numbers read from the file (literally, no rounding):
- `R_voidfloor_pointwise_linear_compression_inert` = 0.25 ratio (bending/stretching with matched top stretching)
  - `R_pointwise_linear_compression_symmetric` = 0.5 ratio
  - `R_pointwise_threshold_at_half_peak` = 0.125 ratio
  - `R_pointwise_quadratic_tension_only` = 0.16666666666666666 ratio
  - `R_peak_driven_network_integrates` = 1.0 ratio
  - `collagen_increase_cyclic_flexure_vs_STATIC (ovint SMC pa PGA/PLLA)` = 63.0 %
  - `effective_stiffness_increase_flex_arm` = 429.0 %
  - `effective_stiffness_increase_static_arm` = 351.0 %
- What a cell can calculate: Thickness-integrated growth law: how much collagen synthesis and effective stiffness that bending provides compared to pure stretching at the same peak strain, i.e. which mode of loading builds tissue. The entries carry both derived ratios and measured percentages using PMID.

### 7. `source_repository/data/tissuetwin_recovery/recovery_evidence.json`
- Size: 34258 prey · Trace: (b) tissue behavior
- Covered by cell: no
- Numbers read from the file (literally, no rounding):
  - `gate_b.7.off` = 0.00930209402892701 fractional change
  - `gate_b.90.off` = 0.044920874571524605 fractional change
  - `gate_b.365.off` = 0.051584989856551244 fractional change
  - `ceiling.k_basal_per_day` = 0.0033333333333333335 1/day
  - `ceiling.ceiling_lo` = 0.00155 1/day
  - `ceiling.ceiling_hi` = 0.004766666666666667 1/day
  - `ceiling.basal_over_strict` = 2.1505376344086025 ratio
  - `grid[0].tau_cat_h` = 18.0 h
  - `grid[0].tau_ana_h` = 36.0 h
  - `grid[0].gross_to_net` = 38.26582426416812 ratio
- What a cell can calculate: Healing and recovery over time: net tissue build-up as a function of session frequency (every_days, session_min, sets_per_wk) and of catabolic/anabolic time constant. Caps how fast tissue can be built and when higher frequency becomes counterproductive — hi_pct is negative at 25 sets per week.
- Note: Directory also has cadence_evidence.json and cadence_noclip.json, not opened.

### 8. `source_repository/data/erythrocyte_membrane/erythrocyte_membrane_results.json`
- Size: 16467 prey · Track: (b) tissue behavior
- Covered by cell: no
- Numbers read from the file (verbatim, no rounding):
  - `1_shear_modulus.parameters.Lmax_nm` = 190.0 nm
  - `1_shear_modulus.parameters.ell0_native_nm` = 46.0 nm
  - `1_shear_modulus.parameters.p_modeling_nm` = 7.5 nm
  - `1_shear_modulus.parameters.p_measured_nm` = 10.0 nm
  - `1_shear_modulus.parameters.mu_anchor_band_uN_m` = [6.6, 8.3] uN/m
  - `mu0_native_x0_uN_m.p_modeling_7.5nm` = 2.795649939374008 uN/m
  - `mu0_native_x0_uN_m.p_measured_10nm` = 2.0967374545305058 uN/m
  - `void_floor_grid.frac_in_tight_anchor_6to9` = 0.05988455988455989 fraktion av 1386 punkter
  - `2_bilayer_adversary.mu_bilayer_uN_m` = 0.0 uN/m
- What a cell can calculate: Shear modulus of a spectrin network out of WLC-parameters (persistence length, contour length, prestretch) — a bottom-up elasticity cell generalizable to other network tissues, with built-in control that the appendix itself yields zero shear modulus.

### 9. `source_repository/data/calcium_pth_vitd/calcium_pth_vitd_results.json`
- Size: 40448 prey · Track: (b) tissue behavior
- Covered by cell: no
- Numbers read from the file (verbatim, no rounding):
  - `falsifier_1_setpoint.task_band_mmol_l` = [1.1, 1.2] mmol/L
  - `brown_setpoint_controls_mmol_l` = 1.13 mmol/L
  - `parfitt_setpoint_controls_mmol_l` = 1.25 mmol/L
  - `brown_vs_parfitt_correlation_controls_r` = 0.85 korrelation
  - `brown_vs_parfitt_correlation_patients_r` = 0.91 korrelation
  - `brown_setpoint_patients_1oHPT_mmol_l` = 1.32 mmol/L
  - `parfitt_setpoint_patients_1oHPT_mmol_l` = 1.42 mmol/L
  - `falsifier_2_clamp_perturbation.step_mmol_l` = 0.05 mmol/L
- What a cell can figure out: PTH–calcium setpoint as a sigmoid regulator: the response to a calcium-clamp step, and how a setpoint shift (healthy vs. primary HPT) alters bone resorption flux. The entrance to a bone remodeling cell with speech.

### 10. `source_repository/data/msk_smoketest/moment_arm_validation/moment_arm_validation_evidence.json`
- Size: 283742 prey · Track: (b) tissue behavior
- Covered by cell: no
- Numbers read from the file (verbatim, no rounding):
  - `_meta.sign_check_vasint_r_30deg_cm` = -5.054691308100564 cm
  - `sweeps.knee_angle_r.range_deg` = [0.1, 139.89999999999992] grader
  - `sweeps.knee_angle_r.step_deg` = 2.0 grader
  - `groups.vasti.mean_peak_cm` = -5.256708453956803 cm
  - `groups.vasti.min_peak_cm` = -5.403911030681704 cm
  - `groups.vasti.max_peak_cm` = -5.131499657704872 cm
  - `groups.vasti.mean_angle_of_peak_deg` = 12.100000000000001 grader
  - `groups.rectus_femoris.mean_peak_cm` = -5.6368169073348735 cm
- What a cell can calculate: Torque arm as a function of joint angle per muscle group — converts muscle force into joint torque over the entire range of motion. Base for a muscle-joint cell and for sanity checking external MSK models.
- Note: msk_smoketest has 469 files; only this one opened.

### 11. `source_repository/data/pulmonary_surfactant/pulmonary_surfactant_results.json`
- Size: 15749 prey · Track: (b) tissue behavior
- Covered by cell: no
- Numbers read from the file (verbatim, no rounding):
  - `schurch_1982_calibration.area_ratio_from_pct_TLC_via_A~V^2_3` = 1.4521964333909259 ratio
  - `schurch_1982_calibration.gamma_ratio` = 10.0 ratio
  - `schurch_1982_calibration.stability_threshold_n` = 0.5 dimensionless
  - `part1_laplace_table.r0_m` = 0.0001 m
  - `bare_saline_rows[0].gamma_mN_m` = 50.0 mN/m
  - `bare_saline_rows[0].P_Pa` = 1000.0 Pa
  - `bare_saline_rows[0].P_cmH2O` = 10.197162129779281 cmH2O
  - `bare_saline_rows[1].gamma_mN_m` = 60.0 mN/m
- What a cell can calculate: Laplace pressure in the alveolus as a function of surface tension and radius, plus the stability criterion dgamma/dA. Calculates when alveoli collapse without surfactant — mechanical interface between tissue and fluid.

### 12. `source_repository/data/bcell_affinity_maturation/bcell_affinity_maturation_results.json`
- Size: 13568 loot · Track: (c) immune system
- Covered by cell: no
- Numbers read from the file (verbatim, no rounding):
  - `part1_shm_rate.mu_shm_per_bp_per_division` = 0.001 mutations/bp/division
  - `part1_shm_rate.L_V_bp` = 500.0 bp
  - `part1_shm_rate.N_gc_divisions` = 10.0 delningar
  - `part1_shm_rate.E_mutations_per_lineage` = 5.0 mutationer
  - `part1_shm_rate.baseline_replication_rate_per_bp_per_division` = 1e-09 mutationer/bp/delning
  - `part1_shm_rate.fold_vs_baseline` = 1000000.0 ggr
  - `part2...kd0_germline_M` = 1e-06 M
  - `part2...kd_ceiling_anchor_M` = 1e-10 M
- What a cell can calculate: Affinity maturity: number of mutations per line in the germinal center and how far Kd can be driven (1e-6 to 1e-10 M) before the ceiling is reached. A cell can calculate how many GC-rounds are required for a desired antibody affinity and where selection saturates — directly relevant to underactivation, i.e. insufficient affinity.

### 13. `source_repository/data/complement_cascade/complement_cascade_results.json`
- Size: 20592 prey · Track: (c) immune system
- Covered by cell: Yes. IMMUNITY/complement_cascade.py and IMMUNITY/complement_ap_bistability_eculizumab.py are in the workspace cell list, so the domain is covered by filename (this deviates from the background task). The data file itself is not yet linked to any cell.
- Numbers read from the file (literally, no rounding):
  - `constants.T_HALF_UNSTABILIZED_MIN` = 3.0 min
  - `constants.T_HALF_ACTIVATOR_MIN` = 30.0 min
  - `constants.T_HALF_HOST_MIN` = 1.05 min
  - `constants.FOLD_H` = 4.761904761904762 ggr
  - `constants.KM_C5_NM` = 1400.0 nM
  - `constants.KCAT_C5_PER_MIN` = 0.288 1/min
  - `constants.C5_0_NM` = 400.0 nM
  - `constants.K_TICK_per_min` = 0.00016666666666666666 1/min
- What a cell can calculate: C3b amplification and flipping between quiescent and attack states: the half-life of C3b on host surface (1.05 min) versus activator surface (30 min) gives the overactivation threshold, and Km/kcat of C5 gives the MAC flux. Calculating the threshold for complement-driven tissue damage and the effect of C5 blockade.

### 14. `source_repository/data/renal_c3g_joint_lesion_probability/results.json`
- Size: 3985 prey · Track: (c) immune system
- Covered by cell: no
- Numbers read from the file (literally, no rounding):
  - `inputs_cited.c3nef_rate_pct.DDD` = 80.0 %
  - `inputs_cited.c3nef_rate_pct.C3GN` = 45.0 %
  - `inputs_cited.mutation_rate_blended_pct` = 17.9 %
  - `inputs_cited.ic_mpgn_comparator_pct` = 53.06 %
  - `sweep_by_substratum.DDD.frechet_lower_bound_pct` = 80.0 %
  - `sweep_by_substratum.DDD.frechet_upper_bound_pct` = 97.9 %
  - `sweep_by_substratum.DDD.rho_sweep[0].joint_pct` = 83.58 %
- What a cell can figure out: Coprobability of two independent complement defects (autoantibody plus mutation) given only the margins, with Fréchet bounds and correlation sweeps. A cell can calculate what proportion of a patient population carries dual overactivation cause — a pure overactivation question in track (c).

### 15. `source_repository/data/hematopoiesis/hematopoiesis_results.json`
- Size: 35618 prey · Track: (c) immune system
- Covered by cell: no
- Numbers read from the file (literally, no rounding):
  - `part1_quiescence_falsifier.task_band_days` = [145.0, 193.0] day
  - `wilson2008.interval_days` = 145.0 day
  - `wilson2008.divisions_per_lifetime` = 5.0 divisions
  - `foudi2009.rate_range_per_day` = [0.008, 0.018] 1/day
  - `foudi2009.implied_interval_days_range` = [55.6, 125.0] day
  - `forced_adversary_wilson_own_bulk_claim.adversary_days_range` = [21.0, 42.0] day
  - `measured_dormant_subset_days` = 145.0 day
- What a cell can figure out: Stem cell turnover: division interval and divisions per lifetime for HSC, which sets how quickly the immune system can be replenished after depletion. The production ceiling behind underactivation.

### 16. `source_repository/data/iron_hepcidin/iron_hepcidin_results.json`
- Size: 27924 prey · Track: (c) immune system
- Covered by cell: no
- Numbers read from the file (literally, no rounding):
  - `couples_to_siblings_readonly.rbc_lifespan_days` = 115.0 day
  - `couples_to_siblings_readonly.blood_volume_l` = 5.0 L
  - `step1_iron_pools.fe_mg_per_g_hb_derived` = 3.4658 mg Fe/g Hb
  - `step1_iron_pools.hb_tetramer_mw_g_mol_derived` = 64453.2 g/mol
  - `step1_iron_pools.lit_hb_fraction_of_total_band` = [0.667, 0.75] fraction
  - `step2_mass_balance_falsifier_F1.recycling_mg_day_derived` = 22.603 mg/day
  - `lit_recycling_mg_day_band` = [20.0, 25.0] mg/day
  - `absorption_mg_day_lit_mid` = 1.5 mg/day
- What a cell can calculate: Iron mass balance: feedback from macrophages' RBC-breakdown (22.603 mg/day) towards intestinal absorption (1.5 mg/day). A cell can figure out how long a given blood loss takes to compensate and where the hepcidin brake binds — linking immune cells to systemic physiology.

### 17. `source_repository/data/erythropoiesis/erythropoiesis_results.json`
- Size: 17664 prey · Track: (c) immune system
- Covered by cell: no
- Numbers read from the file (literally, no rounding):
  - `couples_to_siblings_readonly.hct_ref_pct` = 45.0 %
  - `cao2_ref_ml_dl` = 19.758710131993222 mL O2/dL
  - `rbf_l_min` = 1.0768 L/min
  - `renal_o2_delivery_ml_min` = 212.76 mL/min
  - `step1_rbc_lifespan.mpl_biotin_low_densities_days` = [115.0, 113.0] day
  - `t50_biotin_low_densities_days` = [58.0, 57.0] day
  - `mpl_biotin_single_density_1999_days` = 103.0 day
  - `mpl_51cr_elution_corrected_days` = 116.0 day
- What a cell can calculate: RBC lifespan from multiple independent measurements plus renal O2 delivery as EPO-signal. A cell can figure out how hematocrit responds to a hypoxic stimulus and which measurement method is binding.

### 18. `source_repository/data/socket_contact_winkler/evidence.json`
- Size: 47667 prey · Track: (d) surgery and medical robots
- Covered by cell: no
- Numbers read from the file (literally, no rounding):
  - `inputs.F_bayonet_governing_N` = 11.849592622820806 N
  - `inputs.M_bayonet_Nm` = 1.9302986382575094 Nm
  - `sweeps.E_tissue_kPa_list_ENGINEERING_ASSUMPTION` = [50.0, 600.0] kPa
  - `sweeps.T_TISSUE_MM_ENGINEERING_ASSUMPTION` = 8.0 mm
  - `sweeps.k_list_Pa_per_m` = [6250000.0, 75000000.0] Pa/m
  - `sweeps.worst_case.p_max_kPa` = 113.35360117736319 kPa
  - `sweeps.worst_case.peak_to_mean` = 3.3352531669326084 kvot
  - `sweeps.worst_case.contact_frac` = 0.49250000000000005 fraktion
  - `sweeps.nominal_case.p_max_kPa` = 5.600922791759508 kPa
  - `mesh_convergence_self_test.finest_vs_prev_relchange` = 0.00025526541691094534 relative
- What a cell can calculate: Tool to tissue: a Winkler bed gives contact pressure p_max and contact proportion from applied force and moment and the tissue stiffness E_tissue. A cell can calculate the grip force a robotic cuff must apply before the pressure ulcer limit is crossed, and how sensitive p_max is to tissue stiffness (k_invariance_relspread 2.7e-15).
- Note: Mesh convergence and comparison against closed form are in the same file as self-tests.

### 19. `source_repository/data/socket_contact_edge_singularity/evidence.json`
- Size: 14095 prey · Track: (d) surgery and medical robots
- Covered by cell: no
- Numbers read from the file (verbatim, no rounding):
  - `inputs.F_bayonet_governing_N` = 11.036249999999997 N
  - `inputs.M_bayonet_Nm` = 1.7978051249999996 Nm
  - `inputs.b_nominal_m` = 0.12645 m
  - `self_test_0_naive_port.true_answer_p_centre_Pa` = 1591549.4309189534 Pa
  - `self_test_0_naive_port.true_answer_p_near_edge_Pa` = 5097037.441251784 Pa
  - `self_test_0_naive_port.naive_port_p_centre_Pa` = -0.0 Pa
  - `self_test_2_winkler_limit_reduction.uniform_kernel_p_Pa` = 1268.2685123035417 Pa
  - `self_test_1_kernel_quadrature_nominal.F_relerr` = 1.6095656037152574e-16 relative
- What a cell can figure out: The edge singularity in contact between hard and soft: the pressure at the edge is roughly three times the center value (5.097 MPa vs. 1.592 MPa) — just the spike that determines whether a surgical tool or a sheath edge damages tissue. A cell can calculate the edge radius required to keep p below a tissue boundary.

### 20. `source_repository/data/prosth_real_task_grasp_sweep/evidence.json`
- Size: 24911 change · Track: (d) surgery and medical robots
- Covered by cell: no
- Numbers read from the file (verbatim, no rounding):
  - `delivered_force.F0_frictionless_N_at_T_high` = 11.849592622820806 N
  - `delivered_force.mu_cable_central` = 0.35 friktionskoefficient
  - `delivered_force.wrist_retained_fraction_central` = 0.9396468777237291 fraktion
  - `delivered_force.F_delivered_conservative_incl_pin_friction_N_at_T_high` = 9.822621826975045 N
  - `delivered_force.MU_FINGERTIP` = 0.4 friktionskoefficient
  - `geometric_envelope.hinge_wobble_slop_mm` = 5.407451092865333 mm
  - `geometric_envelope.approach_error_deg_chest_robust` = 17.1 grader
  - `tasks.door_key_turn_lock_cylinder.T_turn_Nm_bracket` = [0.05, 0.35] Nm
  - `tasks.door_key_turn_lock_cylinder.sweep.frac_force_closed` = 0.6791979949874687 fraction
- What a cell can calculate: Force transmission through a cable-driven manipulator with frictional losses, and what proportion of an approach error sweep produces force slope on a real object. A cell can calculate the force delivered at the grip tip and the angular tolerance of a tool–object grip — the same math as a surgical instrument in a holder.

### 21. `source_repository/data/socket_two_patch_free_body_resolution/evidence.json`
- Size: 11388 change · Track: (d) surgery and medical robots
- Covered by cell: no
- Numbers read from the file (verbatim, no rounding):
  - `inputs.F_bayonet_governing_N` = 11.849592622820806 N
  - `inputs.M_bayonet_Nm` = 1.9302986382575094 Nm
  - `inputs.cuff_length_mm` = 252.9 mm
  - `inputs.ARC_GAP_DEG` = 60.0 grader
  - `inputs.ceiling_kPa` = 83.5 kPa
  - `q1_free_body_from_cad.contact_arc_deg` = 300.0 grader
  - `gate_b_field_reconstruction.rows[0].reconstructed_p_max_Pa` = 5600.919304805415 Pa
- What a cell can figure out: Exposing two contact surfaces: how force and torque are distributed between two separate tissue patches, and whether the reconstructed pressure field is consistent with static equilibrium. A cell can calculate the load distribution when a tool has two support points against tissue, against a pressure ceiling of 83.5 kPa.

### 22. `source_repository/data/prosth_socket_wholebody_thermal_transient_v1/trajectories.json`
- Size: 462826 prey · Track: (d) surgery and medical robots
- Covered by cell: no
- Numbers read from the file (verbatim, no rounding):
  - `A_ordinary_indoor_22C.baseline.native.t_min` = [1.0, 836.0] min (168 points)
  - `B_hot_day_rest_32C.baseline.native.T_core_c` = [36.618250520833335, 36.782505208333305] degC (168 points)
- What a cell can calculate: Thermal transient in tissue under a non-breathable interface: core temperature and skin temperature above 836 min for three interfaces (baseline, gel_liner, ventilated) in four environments. A cell can calculate how long a tool or interface can be attached before the tissue's temperature limit is reached.
- Note: Only the t_min and T_core_c vectors for the first two scenarios opened; the file has more channels (native, grafted).

### 23. `source_repository/data/deglutition_swallowing/deglutition_swallowing_results.json`
- Size: 34835 prey · Track: (e) other physiology
- Covered by cell: no
- Numbers read from the file (verbatim, no rounding):
  - `margin_computation.apnea_duration_s` = 1.0 s
  - `margin_computation.ues_relax_dry_s` = 0.37 s
  - `margin_computation.ues_relax_20ml_s` = 0.65 s
  - `margin_computation.margin_dry_swallow_s` = 0.63 s
  - `margin_computation.margin_20ml_swallow_s` = 0.35 s
  - `tau_close_sweep.critical_tau_close_s` = 0.35 s
  - `random_onset_monte_carlo.tau_close_representative_s` = 0.15 s
- What a cell can calculate: Timing margin between airway closure and UES-relaxation. A cell can calculate aspiration risk as a function of bolus volume and closure delay — a pure margin cell with seconds as a unit.

### 24. `source_repository/data/mucociliary_clearance/mucociliary_clearance_results.json`
- Size: 12809 prey · Track: (e) other physiology
- Covered by cell: no
- Numbers read from the file (verbatim, no rounding):
  - `healthy.f0_hz` = 13.0 Hz
  - `healthy.x0_delta_over_L` = 0.9285714285714286 dimensionslos
  - `healthy.U_um_s` = 110.5 um/s
  - `healthy.U_mm_min` = 6.63 mm/min
  - `healthy.measured_band_mm_min` = [4.0, 10.0] mm/min
  - `cbf_sweep_measured_values.cbf_hz` = [8.0, 16.8] Hz (13 punkter)
  - `cbf_sweep_measured_values.U_mm_min` = [4.08, 8.568] mm/min (13 punkter)
  - `slope_check.model_slope_um_s_per_hz` = 8.5 (um/s)/Hz
- What a cell can calculate: Transport speed from cilia beat frequency and mucus layer thickness. A cell can calculate the clearance time for a particle and how much frequency drop is needed to knock out the transport.

### 25. `source_repository/data/fluid_compartments/fluid_compartments_results.json`
- Size: 25919 prey · Track: (e) other physiology
- Covered by cell: no
- Numbers read from the file (literally, no rounding):
  - `man_73kg.skeletal_muscle (mass kg, ECF-water L)` = [28.5, 3.47] kg and L respectively
  - `man_73kg.skin` = [4.3, 2.94] kg and L respectively
  - `man_73kg.brain_viscera` = [6.0, 1.51] kg and L respectively
  - `man_73kg.bone` = [5.75, 0.86] kg and L respectively
  - `man_73kg.active_marrow` = [1.17, 0.19] kg and L respectively
  - `man_73kg.connective_tissue` = [3.7, 2.96] kg and L respectively
  - `man_73kg.transcellular` = [1.1, 1.05] kg and L respectively
- What a cell can calculate: ICF/ECF-partition per tissue for a 73 kg reference. A cell can calculate volume of distribution for a substance and osmotic volume shifts per organ—the basis of all transport cells.

