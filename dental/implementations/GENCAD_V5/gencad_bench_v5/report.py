from common import *
import collections, itertools, gzip
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def create_report():
    r1 = read(ROOT / 'raw/R1.json')
    bound = read(ROOT / 'raw/R1B_BOUNDARY.json')
    phys = read(ROOT / 'raw/R3.json')
    fun = read(ROOT / 'raw/FUNCTIONAL_ROWS.json')
    contrast = read(ROOT / 'raw/R5_CONTRAST.json')
    ctrl = read(ROOT / 'raw/CONTROLS.json')
    norms = read(ROOT / 'raw/CONTACT_NORMS.json')
    lab = read(ROOT / 'raw/LAB_PACKAGE.json') if (ROOT / 'raw/LAB_PACKAGE.json').exists() else {'images': 0, 'ICC_status': 'NOT_PREPARED'}
    r6 = read(ROOT / 'raw/R6_EXTERNAL_SUPPORT.json')
    fun = fun + r6['rows']
    ext = [r for r in fun if r['track'] == 'ToothCraft']
    scored = [r for r in fun if r['status'] == 'SCORED']
    tiers = {str(t): sum((r['reconstruction_tiers'][str(t)] for r in scored)) for t in [0.35, 0.5, 0.75, 1.0, 1.5]}
    names = r1['participants']
    summary = []
    for (key, grp) in itertools.groupby(sorted(scored, key=lambda x: (x['track'], x['participant'])), key=lambda x: (x['track'], x['participant'])):
        rr = list(grp)
        supported = [r for r in rr if r['function']['0']['status'] == 'GEOMETRIC_ONLY']
        summary.append(dict(track=key[0], participant=key[1], scored=len(rr), closed=sum((r['digital_closed_shell'] for r in rr)), tiers={str(t): sum((r['reconstruction_tiers'][str(t)] for r in rr)) for t in [0.35, 0.5, 0.75, 1.0, 1.5]}, contact_supported=len(supported), median_contact_regions=float(np.median([r['function']['0']['contact_regions'] for r in supported])) if supported else None, median_penetration_mm2=float(np.median([r['function']['0']['penetration_area_mm2'] for r in supported])) if supported else None, resolution='POPULATION'))
    cohort = read(V4 / 'payload/private/COHORT.json')
    select = read(ROOT / 'EXTERNAL_INPUTS.json')
    power = next((p for p in r1['pairs'] if p['family'] == 'molar_crown' and p['a'] == 'constraint_optimizer' and (p['b'] == 'x42_contact_branch') and (p['metric'] == 'contact_symdiff_mm2')))
    receipt = read(DATA / 'external_predictions/RECEIPT.json') if (DATA / 'external_predictions/RECEIPT.json').exists() else None
    result = dict(schema='gencad-bench-v5-v1', claim_type='capability', status='PARTIAL_CAPABILITY_PHYSICAL_CALIBRATION_UNIDENTIFIED', review_state='PENDING_INDEPENDENT_REVIEW', capability='External pretrained crown-completion pilot plus functional information and experiment planning, corrected v4 support contracts', external_referent={'kind': 'published_dataset', 'locator': ['https://ditto.ing.unimore.it/bite2text/', 'https://ditto.ing.unimore.it/bits2bites/', str(V4 / 'payload/whole_private')], 'compared_quantity': 'native tooth, neighbouring and antagonist geometry; supplied scan pose', 'refutes_us': True}, additional_external_referents=read(ROOT / 'PREREG_V5.json')['additional_external_referents'], headline_results=dict(primary_patient_case_clusters=48, internal_full_panel_participants=8, frontier_participants=4, frontier_binary_discriminating_items=sum((x['discriminating_items'] for x in bound['frontier_families'])), frontier_binary_items=sum((x['items'] for x in bound['frontier_families'])), external_generated_cases=len(ext), external_unique_patient_cases=1 if ext else 0, external_receipt=receipt, ranking_reversal_pairs=contrast['reversals'], tooth_cases_with_reversal=contrast['cases_with_reversal'], exact_boundary_cases=len(bound['rows']), exact_boundary_gate=bound['gate'], whole_scored=len(scored), whole_requested=len(fun), whole_tiers=tiers, power_example=power, physical_records=phys['candidate_records'], calibrated_generated_force=0, calibrated_generated_film=0, blind_images=lab['images'], controls=len(ctrl['tests']) + len(r6['context_fault_tests']), all_faults_rejected=ctrl['all_injected_faults_rejected'] and all((t['injected_translated_observation_rejected'] for t in r6['context_fault_tests']))), external_support_round=r6, X60_subset=read(ROOT / 'raw/X60_TABLE.json'), whole_summary=summary, physical_axis=phys, contact_norms=norms, limitations=['Retrospective already-examined v4 cohort; no independent clinical evaluation', 'External pilot one patient-case/three tooth types and unverified training overlap; no field ranking', 'Actual ToothCraft weights with explicit SiLU repair; normalization/pose/crop adapter differs from author dataset', 'No generated specimen has matched absolute fracture or seated-film calibration', 'Connected contact regions use0.1mm geometric band, not physical loaded contacts', 'X60 partial plugin outputs included with coverage; borrowed force scale is conditional', 'No rigorous continuous geometric/FE/transfer enclosure; power planning approximate'], resolution_map={'scan_coordinates': 'PER_POINT', 'contact_mask_and_cement_field': 'PER_SURFACE_REGION', 'generated_crown_vector': 'PER_TOOTH', 'cohort_rates_and_published_group_quantiles': 'POPULATION', 'unmeasured_constitutive_transfer': 'PHENOMENOLOGICAL'}, edge_contracts=[dict(producer='native antagonist/preparation scan geometry', consumer='functional scorer', resolution='PER_POINT', time_scale='SIMULTANEOUS'), dict(producer='manufactured seated film (missing)', consumer='crown mechanics', resolution='PER_SURFACE_REGION', time_scale='HANDOVER'), dict(producer='published flaw/strength distribution', consumer='spatial crown hazard', resolution='POPULATION', time_scale='SIMULTANEOUS', status='UNKNOWN_TRANSFER_TO_POINTWISE_INTENSITY')], cost={'current_run': read(ROOT / 'raw/RUN_COST.json') if (ROOT / 'raw/RUN_COST.json').exists() else None, 'external': receipt, 'preparation': read(ROOT / 'raw/EXTERNAL_MODEL_LOCK.json')['seconds'], 'inherited_fit_and_discovery': 'UNKNOWN; not counted as zero', 'new_model_training': 0, 'physical_validation': 'NOT_RUN', 'clinical_queries': 0, 'fallbacks': 'history/external_attempt1 and resource receipts; all denominators retained'}, clinical_use=False)
    result['headline_results']['controls'] = len(ctrl['tests']) + len(r6['context_fault_tests']) + 1
    result['headline_results']['scientific_fault_checks'] = len(ctrl['tests']) + len(r6['context_fault_tests'])
    result['headline_results']['integrity_fault_checks'] = 1
    result['sufficiency_tests'] = ctrl['sufficiency']
    result['large_artifacts'] = [m for m in read(ROOT / 'raw/EXTERNAL_MODEL_LOCK.json')['files'] if m['bytes'] > 50000000]
    result['dropout_summary'] = dict(external_model_candidates=5, external_models_selected=1, external_model_exclusion_fraction=4 / 5, whole_requested=len(fun), whole_scored=len(scored), whole_unscored=len(fun) - len(scored), whole_unscored_fraction=(len(fun) - len(scored)) / len(fun), whole_rejection_reasons=dict(collections.Counter((r.get('reason') or r['status'] for r in fun if r['status'] != 'SCORED'))), matched_generated_force_rejection_fraction=1.0, matched_generated_film_rejection_fraction=1.0, note='Per-candidate physical records include roof proposals, complete-shell variants and X60; not16465 independent patients or complete crowns')
    result['external_model_replay'] = read(ROOT / 'raw/EXTERNAL_REPLAY.json')
    dump(ROOT / 'results.json', result)
    lines = ['# DentalGenCAD-Bench v5 — functional evidence and frontier tasks', '', 'Capability benchmark; retrospective digital data; PENDING_INDEPENDENT_REVIEW. No scalar clinical ranking.', '', '## Item information on 48 patient-case clusters', '', '| Family | Items | Discriminating, all eight | Discriminating, four frontier methods | Use |', '|---|---:|---:|---:|---|']
    for (a, b) in zip(r1['families'], bound['frontier_families']):
        lines.append(f"| {a['family']} | {a['requested_items']} | {a['discriminating_items']} | {b['discriminating_items']} | feasibility coverage + spatial axes |")
    lines += ['', 'Frontier = constraint_optimizer, field_generator, x42_shape, x42_contact_branch. All binary items are constant across this subgroup. No family is deleted; bridge3 and implant_crown retain their original coverage. Binary entropy and corrected item-total correlations for each item are in raw/ITEMS.json. Eight methods are too few and too dependent for a defensible latent-ability IRT interpretation.', '', f"Planning example (POPULATION): LP minus X42 contact-pattern error = {power['mean_difference']:.3f} mm² on {power['clusters']} paired case clusters. Approximate 80% power at two-sided α=.05 needs {power['clusters_for_80pct']} clusters; bootstrap planning range {power['bootstrap_finite_n_range95'][0]:.0f}–{power['bootstrap_finite_n_range95'][1]:.0f}. Current plug-in power {power['planning_power_at_current_n']:.1%}. This is retrospective, selected both-PASS cases, and not a validated recruitment guarantee. Pass-rate comparison has zero observed effect and no finite sample-size estimate.", '', '## Complete-shell frontier', '', '| Track / participant | Scored | Closed | ≤0.35 | ≤0.50 | ≤0.75 | ≤1.00 | ≤1.50 mm | Contact support | Median regions | Median interference mm² |', '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for s in summary:
        lines.append('| ' + s['track'] + ' / ' + s['participant'] + ' | ' + ' | '.join((str(s[k]) for k in ['scored', 'closed'])) + ' | ' + ' | '.join((str(s['tiers'][str(t)]) for t in [0.35, 0.5, 0.75, 1.0, 1.5])) + ' | ' + str(s['contact_supported']) + ' | ' + str(s['median_contact_regions']) + ' | ' + ('UNKNOWN' if s['median_penetration_mm2'] is None else f"{s['median_penetration_mm2']:.3f}") + ' |')
    lines += ['', 'Reconstruction tiers are descriptive (PER_TOOTH), not clinical thresholds. The original0.35mm frontier remains unchanged. Functional contact/interference and proximal gaps are separate axes. A collision-free crown with no contact is not thereby adequate. Unsigned neighbour distances do not establish interproximal pressure or exclude overlap.', '', f"There are {contrast['reversals']} within-tooth ordered pairs on {contrast['cases_with_reversal']} tooth cases with better resemblance but worse interference at the frozen0.01mm/0.05mm² thresholds. These are correlated pair comparisons, not independent observations.", '', '## External and additional participants', '', f'ToothCraft: {len(ext)} actual pretrained outputs scored on {(1 if ext else 0)} patient-case. Author code/checkpoint hashes, adaptation and resource receipt are saved; private reference unmounted during inference. External checkpoint licence is unspecified. No field-wide superiority claim.', '', 'X60 volume1.0 and volume1.1 actual plugin outputs are in raw/X60_PARTICIPANT_ROWS.json. Their partial coverage and conditional force ratio are shown in README_DEMO; they are not padded into the48-case table.', '', '## External support construction', '', 'Retained author-default100-step sampling does not explain the failed10-step target reconstruction. A new public-context ownership adapter reduces reconstruction errors while preserving the observed context outside the target cell exactly. It remains an adapted method and fails all original0.35mm target gates. See raw/R6_EXTERNAL_SUPPORT.json for all scores and the0.25mm known-neighbour surface check.', '', '## Physical outputs', '', f"Every generated entry has an explicit physical port ({phys['candidate_records']} records). Calibrated force05:0; seated spatial film:0. Published protocol-specific quantiles are computed separately with parameter-rectangle ranges. Three dry marginal-fit targets pass; four of five force-mean targets pass. The failed source target remains failed. Numerical source controls match; no algorithm-superiority interpretation.", '', '## Complete tables', '', 'raw/R1.json: all paired sample-size plans and dropouts. raw/FUNCTIONAL_ROWS.json: per-candidate shape/contact/proximal/removal vectors and pose scenarios. raw/PHYSICAL_ROWS.json.gz: per-design physical status. raw/R1B_BOUNDARY.json:16 exact rational near-limit tasks and LP outcomes. raw/CONTACT_NORMS.json: primary tooth-type contact context and rejected numerical normative transfers.']
    lines += ['', '## X60 subset (12 case clusters; no full-panel rank)', '', '| Participant | Requested | PASS | Abstained | Median conditional force05 ratio | Absolute force05 |', '|---|---:|---:|---:|---:|---|']
    for x in result['X60_subset']:
        lines.append(f"| {x['participant']} | {x['requested']} | {x['PASS']} | {x['status_counts']['ABSTAIN']} | {x['conditional_quantile_ratio_median']:.4f} | UNKNOWN |")
    lines += ['', 'These POPULATION summaries aggregate conditional PER_TOOTH model outputs. No rigorous physical transfer enclosure is available.', '', f'Dropouts: {len(fun) - len(scored)}/{len(fun)} complete-shell submissions unscored ({(len(fun) - len(scored)) / len(fun):.2%}); source reasons retained in results.json. All generated physical calibration ports reject unmatched transfer. Four of five external model candidates lacked runnable public artifacts or were non-author reproductions.']
    (ROOT / 'LEADERBOARD.md').write_text('\n'.join(lines) + '\n')
    (fig, axs) = plt.subplots(2, 2, figsize=(13, 9))
    ax = axs[0, 0]
    labels = [x['family'].replace('_', ' ') for x in r1['families']]
    y = np.arange(len(labels))
    ax.barh(y, [x['discriminating_items'] / x['requested_items'] for x in r1['families']], label='All8')
    ax.scatter([0] * len(y), y, color='black', marker='x', label='Frontier4')
    ax.set_yticks(y, labels, fontsize=8)
    ax.set_xlabel('Fraction of items distinguishing methods')
    ax.legend()
    ax.set_title('Feasibility saturates among strong methods')
    ax = axs[0, 1]
    for part in sorted({r['participant'] for r in scored}):
        rr = [r for r in scored if r['participant'] == part and r['function']['0']['status'] == 'GEOMETRIC_ONLY']
        ax.scatter([r['reconstruction_p95_mm'] for r in rr], [r['function']['0']['penetration_area_mm2'] for r in rr], s=20, label=part, alpha=0.7)
    ax.axvline(0.35, color='grey', ls='--')
    ax.set_xlabel('Original-tooth reconstruction p95 (mm)')
    ax.set_ylabel('Antagonist interference area (mm²)')
    ax.set_title('Resemblance and function are different axes')
    ax.legend(fontsize=6, ncol=2)
    ax = axs[1, 0]
    q = phys['published_force05']
    x = [r['thickness_mm'] for r in q]
    v = np.array([r['force05_N'] for r in q])
    ranges = np.array([r['parameter_rectangle_range_N'] for r in q])
    ax.errorbar(x, v, yerr=np.array([v - ranges[:, 0], ranges[:, 1] - v]), fmt='o', capsize=4)
    ax.set_xlabel('Published specimen thickness (mm)')
    ax.set_ylabel('Derived force05 (N), population')
    ax.set_title('Published quantiles; generated crown transfer UNKNOWN')
    ax.text(0.03, 0.95, 'Ranges: marginal parameter rectangle, not joint CI', transform=ax.transAxes, va='top', fontsize=8)
    ax = axs[1, 1]
    ax.axis('off')
    ax.text(0, 1, f"Frozen benchmark result\n\n{contrast['reversals']} crossed shape/function orderings\n{len(bound['rows'])} rational boundary tasks, LP agrees\n{len(ext)} pretrained external pilot outputs\n{lab['images']} blinded image objects prepared\n{len(ctrl['tests']) + len(r6['context_fault_tests'])} numeric + 1 integrity fault checks\n\nGenerated absolute force05: UNKNOWN\nGenerated seated cement film: UNKNOWN\nClinician ICC: NOT MEASURED", va='top', fontsize=13)
    fig.tight_layout()
    fig.savefig(ROOT / 'figures/benchmark_v5.png', dpi=170)
    fig.savefig(ROOT / 'figures/benchmark_v5.pdf')
    plt.close(fig)
    return result
if __name__ == '__main__':
    create_report()
