import json, collections, hashlib, time, platform, resource
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lab_port import ROOT, sha

def dropout(rows, screen):
    reasons = collections.Counter((r['reason'] for r in screen['records'] if r['decision'] == 'REJECT'))
    reasons['no_keyword'] = screen['all_xml_files'] - screen['keyword_files']
    r1 = collections.Counter()
    r2 = collections.Counter()
    for r in rows:
        if not r['independent_first_state']:
            r1['repeated_recemented_state'] += 1
        elif r['height_mm'] is None:
            r1['unreported_or_unmapped_bonded_height'] += 1
        elif r['total_convergence_deg'] is None:
            r1['unknown_total_convergence_convention'] += 1
        elif r['thermal_cycles'] is None:
            r1['unknown_thermal_cycle_count'] += 1
        if not r['independent_first_state']:
            r2['repeated_recemented_state'] += 1
        elif r['substrate'] != 'titanium':
            r2['non_titanium_substrate'] += 1
        elif r['height_mm'] is None:
            r2['unreported_or_unmapped_bonded_height'] += 1
    studies = [{'study': 'PMC10054898', 'reason': 'Geometry only, no primary force endpoint', 'groups_rejected': 0}, {'study': 'PMC10064688', 'reason': 'Mean/SD only in force figure; pairwise mean differences are not absolute group means', 'groups_rejected': 8}, {'study': 'PMC10315082', 'reason': 'Relevant numeric Table2 is image-only, absent local image; discussion pooled means lack per-group SD; resin endpoints dominated by screw fracture', 'groups_rejected': 6}, {'study': 'PMC10413795', 'reason': 'Table1 explicitly labels MPa but methods describe force; no matched area to resolve unit conflict', 'groups_rejected': 4}, {'study': 'PMC10416713', 'reason': 'Laser debonding time/temperature instead of unassisted axial peak force', 'groups_rejected': 0}, {'study': 'PMC10873732', 'reason': 'Stress reported; tooth-specific tested area unavailable; overlapping jaw subsets not independent new groups', 'groups_rejected': 3}]
    return {'corpus_files': screen['all_xml_files'], 'keyword_files': screen['keyword_files'], 'title_abstract_primary_reads': 18, 'source_screen_rejected': screen['all_xml_files'] - 18, 'source_screen_rejected_fraction': (screen['all_xml_files'] - 18) / screen['all_xml_files'], 'screen_disjoint_reasons': dict(reasons), 'primary_read_studies': 18, 'force_studies_kept': 12, 'read_study_rejected_fraction': 6 / 18, 'primary_read_exclusions': studies, 'candidate_primary_force_or_stress_group_cells': 101, 'groups_kept': 80, 'rejected_measurement_group_cells': 21, 'rejected_group_cell_fraction': 21 / 101, 'R1_model_excluded_groups': sum(r1.values()), 'R1_disjoint_reasons': dict(r1), 'R2_model_excluded_groups': sum(r2.values()), 'R2_disjoint_reasons': dict(r2), 'coverage_note': 'Targeted finite local corpus screening, not a systematic review. Counts exclude repeat tables and dependent jaw subsets. Study rank/aging-identifiability failures are additional gates, not cell deletion.'}

def plot(r2, r3):
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
    (fig, ax) = plt.subplots(1, 3, figsize=(12, 3.9), layout='constrained')
    p = r2['sufficiency_test']
    ax[0].bar([0, 1], [p['force_A_N'], p['force_B_N']], color=['#667788', '#168a84'])
    ax[0].set_xticks([0, 1], ['0 dimples', '4 dimples'])
    ax[0].set_ylabel('Published group mean pull-off force (N)')
    ax[0].set_title('Same height / TOC / cement / aging\nChoi 2023, Table 1')
    ax[0].text(0.5, 680, f"+{p['relative_difference'] * 100:.1f}%\nsummary error = 0", ha='center')
    ax[0].set_ylim(0, 810)
    colors = ['#596c80', '#168a84', '#a45c58', '#bc954b']
    for (i, f) in enumerate(r2['loso']['folds']):
        q = f['predictions']
        x = [p['observed_N'] for p in q]
        y = [p['predicted_group_mean_N'] for p in q]
        ax[1].scatter(x, y, s=25, color=colors[i], label=f['held_out_study'])
    ax[1].plot([40, 1000], [40, 1000], color='#999999', lw=1)
    ax[1].set_xscale('log')
    ax[1].set_yscale('log')
    ax[1].set_xlabel('Published group mean (N)')
    ax[1].set_ylabel('Study-held-out prediction (N)')
    ax[1].set_title(f"External transfer fails\nlogRMSE {r2['loso']['log_rmse']:.3f}; coverage {r2['loso']['95pct_coverage'] * 100:.0f}%")
    ax[1].legend(fontsize=6, loc='lower right')
    vals = [48, 12, 12]
    labels = ['Fracture\nunaged', 'Retention\nsham 37°C', 'Retention\nTC10000']
    ax[2].bar(np.arange(3), vals, color=['#667788', '#168a84', '#a45c58'])
    ax[2].set_xticks(np.arange(3), labels)
    ax[2].set_ylabel('Allocated crowns (count)')
    ax[2].set_title('Same total 72 crowns\nProspective; no lab measurement')
    for (i, v) in enumerate(vals):
        ax[2].text(i, v + 1, str(v), ha='center')
    ax[2].set_ylim(0, 56)
    fig.savefig(ROOT / 'crown_retention.png', dpi=180)
    fig.savefig(ROOT / 'crown_retention.pdf')
    plt.close(fig)

def run():
    rows = json.loads((ROOT / 'raw/RETENTION_GROUPS.json').read_text())
    s = json.loads((ROOT / 'raw/SCREENING.json').read_text())
    r1 = json.loads((ROOT / 'rounds/RESULTS_R1.json').read_text())
    r2 = json.loads((ROOT / 'rounds/RESULTS_R2.json').read_text())
    r3 = json.loads((ROOT / 'rounds/RESULTS_R3.json').read_text())
    v = json.loads((ROOT / 'VERIFICATION.json').read_text())
    drop = dropout(rows, s)
    plot(r2, r3)
    edges = [{'producer': 'raw/RETENTION_GROUPS.json', 'consumer': 'source-conditioned retention model', 'shared_resolution': 'POPULATION', 'timescale': 'SIMULTANEOUS', 'quantity': 'Protocol-specific group mean force_N; no invented pointwise traction.'}, {'producer': 'X1b preparation.stl', 'consumer': 'nominal retention interface port', 'shared_resolution': 'PER_POINT', 'timescale': 'SIMULTANEOUS', 'quantity': 'Triangle normals/centroids/area retained; consumer performs regional aggregation.'}, {'producer': 'published and prospective aging protocol', 'consumer': 'post-aging retention assay', 'shared_resolution': 'PER_TOOTH', 'timescale': 'HANDOVER', 'quantity': 'Specimen history -> interface start state; chemistry UNKNOWN.'}, {'producer': 'prospective terminal force measurements', 'consumer': 'LAB_PROTOCOL ratio query', 'shared_resolution': 'PER_TOOTH', 'timescale': 'SIMULTANEOUS', 'quantity': 'One ultimate endpoint per crown; group means aggregate downstream.'}, {'producer': 'local Ti studies', 'consumer': 'X1b composite die', 'shared_resolution': 'POPULATION', 'timescale': 'SIMULTANEOUS', 'status': 'REFUSED_OUT_OF_SUPPORT', 'quantity': 'No direct force transfer edge admitted.'}]
    out = {'round_tag': 'X47-crown-retention', 'claim_type': 'information_link', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'scientific_admission': False, 'capability': 'Whole-crown retention is now an explicit source-conditioned endpoint with disjoint aging assay on same72 proposed crowns. Calibrated absolute force and aging ratio remain UNKNOWN.', 'external_referent': {'kind': 'independent_measurement', 'locator': 'https://doi.org/10.4047/jap.2023.15.2.63 Table1; all12 DOI/PMCID+table coordinates in raw/RETENTION_GROUPS.json', 'compared_quantity': 'Whole-crown axial pull-off mean force_N, summary-identical local-feature contrast; study-held-out absolute force', 'refutes_us': True}, 'results': {'published_groups': 80, 'published_studies': 12, 'R1_full_covariate_groups': 36, 'R1_full_covariate_studies': 3, 'R1_rank': 5, 'R1_columns': 6, 'R1_gate': r1['gate'], 'R2_titanium_groups': 25, 'R2_studies': 4, 'R2_gate': r2['gate'], 'R2_loso_log_rmse': r2['loso']['log_rmse'], 'R2_loso_95pct_coverage': r2['loso']['95pct_coverage'], 'summary_identity_error': 0.0, 'summary_identical_force_difference_N': r2['sufficiency_test']['difference_N'], 'summary_identical_force_relative_difference': r2['sufficiency_test']['relative_difference'], 'aged_unaged_ratio': 'UNKNOWN', 'matched_numeric_aging_pairs': 0, 'absolute_X1b_retention_force_N': 'UNKNOWN', 'allocation': r3['allocation_counts'], 'all_required_controls_pass': v['pass']}, 'resolution': {'published_force_group_statistics': 'POPULATION', 'individual_assay_force': 'PER_TOOTH', 'nominal_face_geometry': 'PER_POINT', 'nominal_regions': 'PER_SURFACE_REGION', 'model_coefficients': 'PHENOMENOLOGICAL', 'screen_and_allocation_counts': 'POPULATION operational counts; not patient prevalence'}, 'dropout': drop, 'same_information_control': {'kind': 'Conventional marginal GLS with same REML study covariance', 'max_beta_difference': r2['model']['gls_control_max_abs'], 'scope': 'Numerical parity only; information-link value is the new measured endpoint and refuted sufficiency, not algorithm superiority'}, 'sufficiency_test': r2['sufficiency_test'], 'edges': edges, 'phenomenological_debts': [{'quantity': 'log-force coefficients and Gaussian study effects', 'replacement_measurement': 'Whole-crown force traces, exact bonded prep geometry, substrate/cement/surface labels and held-out batches', 'empirical_transfer': 'FAILED'}, {'quantity': 'local traction-separation and post-aging damage field', 'replacement_measurement': 'Interface coupon/region measurements, water/thermal history and failure map; whole-crown peaks alone insufficient'}, {'quantity': 'new assay aged/unaged ratio', 'replacement_measurement': '24 disjoint retention measurements in same material batch and locked protocol, then independent held-out batch'}], 'rigorous_error_enclosure': 'MISSING for empirical response/linear-log model and small-n delta intervals; no affine sensitivity certificate claimed', 'largest_obstacle': 'No supported composite-die absolute force or numeric matched aging contrast; unknown cohesive law and unvalidated tensile grip.', 'next_operation': 'Verify drag adapter; manufacture according to frozen72 allocation; measure12sham+12TC retention, retain competing modes; freeze calibrated contrasts before held-out second batch.', 'frozen_predictions': {'path': 'FROZEN_PREDICTIONS.json', 'sha256': sha(ROOT / 'FROZEN_PREDICTIONS.json'), 'quantitative_candidate_force': 'UNKNOWN', 'quantitative_candidate_aging_ratio': 'UNKNOWN', 'null_ratio': 1, 'is_physical_prediction_validation': False}, 'data_artifacts': [{'path': g['local_face_artifact'], 'bytes': (ROOT / g['local_face_artifact']).stat().st_size, 'sha256': g['local_face_sha256'], 'resolution': 'PER_POINT'} for g in r3['geometry']], 'full_cost': {'preparation': 'Local XML screen, manual extraction and source replay; orientation time and reasoning tokens not metered', 'fit_and_LOSO_wall_s_R1': r1['wall_s'], 'fit_and_LOSO_wall_s_R2': r2['wall_s'], 'geometry_and_allocation_wall_s': r3['timing_wall_s'], 'discovery': 'Three frozen constructions plus preserved screening/metadata/row-index/runtime failures', 'validation': 'External LOSO,40+fault/consistency checks, physical lab not performed', 'query': 'Source support mask and prospective scoring; fresh run timed in RUN_RECEIPT.json', 'fallback': 'No clinical-year conversion or invented retention force; disjoint protocol supplied', 'peak_RSS_MiB_report_process': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'physical_cost': 'UNKNOWN; 72crowns already in proposed batch, reusable tensile adapter and ~194.44h nominal TC duration plus handling'}, 'claims': {'algorithm_superiority': False, 'calibrated_absolute_retention': False, 'calibrated_age_response': False, 'clinical_lifetime': False, 'prospective_assay_delivered': True}, 'figures': ['crown_retention.png', 'crown_retention.pdf']}
    (ROOT / 'results.json').write_text(json.dumps(out, indent=2, ensure_ascii=False) + '\n')
    return out
if __name__ == '__main__':
    r = run()
    print(r['results'])
    print(r['dropout']['R1_disjoint_reasons'], r['dropout']['R2_disjoint_reasons'])
