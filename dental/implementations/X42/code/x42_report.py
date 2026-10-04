"""Report frozen measurements, failure gates and all family/level comparisons."""
from common import *
import math, csv

def ff(x, d=3):
    return 'UNKNOWN' if x is None else f'{x:.{d}f}'

def ci(s, d=3):
    if s is None:
        return 'UNKNOWN'
    (a, b) = s['case_bootstrap_95']
    return f"{ff(s['mean_difference'], d)} [{ff(a, d)}, {ff(b, d)}]"

def report(agg):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    candidate = agg['selected_candidate']
    control = 'control_full_kernel'
    summary = agg['primary_test_summary']
    contrasts = agg['paired_contrasts']
    family = sorted({d['family'] for d in contrasts if d['family'] != 'all'})
    official = agg['official_evaluation']
    g = agg['gates']
    rs = {n: json.load(open(ROOT / 'rounds' / f'{n}.json')) for n in ['R1', 'R2', 'R3']}
    gen = json.load(open(DATA / 'PREDICTIONS/COST.json'))
    fit = json.load(open(ROOT / 'raw/FINAL_FIT.json'))
    isolation = {'first_failed': json.load(open(ROOT / 'raw/ISOLATION_FAILED_01.json')), 'resumed': json.load(open(ROOT / 'raw/ISOLATION_RESUME_RECEIPT.json')), 'extra_control': json.load(open(ROOT / 'raw/FULL_KERNEL_ISOLATION_RECEIPT.json'))}
    full = dict(preparation={'orientation_discovery_wall_seconds': 'UNKNOWN; research reading/code design not instrumented', 'inherited_X11_label_fit': 'UNKNOWN', 'inherited_primary_train_template_fit': 'UNKNOWN; never counted free', 'own_data_byte_count': sum((p.stat().st_size for p in DATA.rglob('*') if p.is_file() and 'replay_evaluation' not in p.relative_to(DATA).parts))}, fit={'R2_all_models_wall_seconds': rs['R2']['fit_seconds'], 'R3_contact_models_wall_seconds': rs['R3']['fit_seconds'], 'final_refit_wall_seconds': sum((r['seconds'] for r in fit['cost'])), 'full_kernel_control_fit_wall_seconds': sum((r['seconds'] for r in json.load(open(ROOT / 'raw/FULL_KERNEL_FIT.json'))['fits']))}, discovery={'R1_wall_seconds': rs['R1']['wall_seconds'], 'R2_wall_seconds': rs['R2']['wall_seconds'], 'R3_wall_seconds': rs['R3']['wall_seconds'], 'R1_failed_first_invocation': 'preserved; precise wall cost UNKNOWN', 'unprojected_model_candidates': {'kernel': 9, 'forest': 3, 'ridge': 4, 'mean': 1}, 'projected_R3_beta_candidates': 8, 'no_claim_of_10x_cost_gain': True}, generation={'resumed_remaining': gen, 'first_failed_attempt': json.load(open(ROOT / 'raw/ISOLATION_FAILED_01.json')), 'first279_prefix_cost': json.load(open(ROOT / 'raw/GENERATION_FIRST_PART_COST.json')), 'full_kernel_control': json.load(open(ROOT / 'raw/full_kernel_cost/COST.json'))}, validation={'official_scorer_wall_seconds': official['seconds'], 'official_scorer_peak_rss_MiB': official['peak_rss_MiB'], 'same_freeze_replay_wall_seconds': json.load(open(ROOT / 'DELIVERY_CHECK_COST.json'))['replay_scorer_wall_seconds'] if (ROOT / 'DELIVERY_CHECK_COST.json').exists() else None, 'public_dev_generator_replay_wall_seconds': json.load(open(ROOT / 'DELIVERY_CHECK_COST.json'))['dev_generator_wall_seconds'] if (ROOT / 'DELIVERY_CHECK_COST.json').exists() else None, 'controls': '12 geometry/schema/native source fault probes with positive restored plus8 contact-band probes', 'export_receipt': 'raw/EXPORT_RECEIPT.json', 'CPU_field_port': json.load(open(ROOT / 'raw/CPU_FIELD_PORT.json'))}, questions={'hidden_joint_queries': 1, 'dev_queries': '24 R1 plus104 R2 plus104 R3, reused cases and full grids explicitly counted', 'no_post_test_model_selection': True}, fallback={'operation': 'Exact constrained LP or ABSTAIN; all calls included in method generation seconds'}, hardware={'gpu_used': False, 'gpu_seconds': 0, 'max_numeric_threads_per_worker': 1, 'max_workers': 3, 'max_total_cpu_capacity': 4, 'actual_cpu_user_sys_seconds': 'UNKNOWN; recorded method seconds are elapsed, not CPU work'}, license={'Bite2Text': 'precise redistribution terms UNKNOWN in local metadata/public primary page; local research only', 'Bits2Bites': 'CC BY-NC-SA per user brief; exact license file not present in V3 metadata', 'ToothFairy2': 'NOT_USED', 'Teeth3DS': 'No meshes used; inherited X11 training lineage/terms remain upstream', 'mandible_defects': 'NOT_USED'})
    results = dict(schema='X42-generative-roof-study-v1', claim_type='algorithm', review_state='PENDING_INDEPENDENT_REVIEW', outcome='PARENT_GATE_FAILED' if not all((d['combined_gate'] == 'PASS' for d in g.values())) else 'PARENT_GATE_PASS_IN_ROOF_SCOPE_ONLY', capability='Executable source-conditioned, contact-state roof generator and exact constrained controls evaluated against published native IOS geometry', selected_candidate=candidate, external_referent=agg['external_referent'], benchmark_sha256=agg['benchmark_sha256'], prediction_freeze_sha256=agg['prediction_freeze_sha256'], primary_test_case_clusters=588, primary_test_tasks=588 * 36, strict_test_summary=summary, gates=g, per_family_contrasts=contrasts, per_family_by_level_contrasts=agg['paired_by_level'], family_level_pass=agg['family_level'], official_evaluation=official, operational_pass_ceiling=json.load(open(ROOT / 'raw/PASS_CEILING_CERTIFICATE.json')) if (ROOT / 'raw/PASS_CEILING_CERTIFICATE.json').exists() else None, development_rounds={k: {j: v for (j, v) in r.items() if j not in ['failed_executions', 'contrasts', 'paired_comparison']} for (k, r) in rs.items()}, dropout={'source_records': official['dropout']['source_records'], 'task_sites': official['dropout'], 'final_fit': {k: fit[k] for k in ['attempted', 'rejected', 'rejection_fraction', 'rejections']}}, full_cost=full, isolation=isolation, scope=agg['scope'], edges=[{'source': 'registered native antagonist triangles', 'consumer': 'constrained generator', 'quantity': 'height inequality by projected overlap vertex', 'resolution': 'PER_POINT', 'timescale': 'SIMULTANEOUS', 'status': 'EXECUTED'}, {'source': 'native source height reference', 'consumer': 'source error and proximity evaluator', 'quantity': 'height mm; proximity state per point', 'resolution': 'PER_POINT', 'timescale': 'SIMULTANEOUS', 'status': 'EXECUTED'}, {'source': 'frozen generated mesh', 'consumer': 'field/manufacture/metrology', 'quantity': 'vertices/faces in mm and local site frame', 'resolution': 'PER_POINT', 'timescale': 'HANDOVER', 'status': 'CPU_EXPORT_EXECUTED_NATIVE_FIELD_NOT_RUN'}], phenomenological_debt=[{'quantity': 'film/preparation/clearance scenario', 'resolution': 'PHENOMENOLOGICAL', 'replace_with': 'registered marked prepared/unprepared specimen and material IFU protocol'}, {'quantity': '0.1mm proximity as contact proxy', 'resolution': 'PHENOMENOLOGICAL', 'replace_with': 'loaded bite/contact map at known registration uncertainty'}, {'quantity': 'whole-shape kernel and probability threshold', 'resolution': 'PHENOMENOLOGICAL', 'replace_with': 'external-domain validation and probability calibration; no force law inferred'}], field_generator_port='FIELD_GENERATOR_PORT.md', data_manifest={'path': str(ROOT / 'DATA_MANIFEST.json'), 'sha256': sha(ROOT / 'DATA_MANIFEST.json'), 'total_bytes': json.load(open(ROOT / 'DATA_MANIFEST.json'))['total_bytes']} if (ROOT / 'DATA_MANIFEST.json').exists() else None, data_artifacts=[{'path': str(p), 'bytes': p.stat().st_size, 'sha256': sha(p)} for p in [DATA / 'PREDICTIONS/FROZEN_PREDICTIONS.json', DATA / 'evaluation/SCORED_ROWS.csv.gz', DATA / 'evaluation/SCORE_DIGEST.json']], negative_results=['R1 contact bias improves contact but worsens height', 'R2 height advantage uncertain versus same-data RF', 'R3 contact-state branch helps controls too', 'Fixed preparation/roof constraints cannot add observed feasible space', 'Full milling, approximate contacts, substance removal and measured contralateral mirror absent from V3; full X42 goal not established'])
    dump(ROOT / 'results.json', results)
    p = ROOT / 'raw/FAMILY_LEVEL.csv'
    p.parent.mkdir(exist_ok=True)
    with p.open('w') as f:
        w = csv.writer(f)
        w.writerow(['dataset', 'family', 'level', 'method', 'case_clusters', 'PASS', 'ABSTAIN', 'FAIL', 'UNKNOWN', 'UNKNOWN_SITE', 'strict_pass_fraction', 'case_bootstrap_low', 'case_bootstrap_high', 'resolution'])
        for r in agg['family_level']:
            w.writerow([r['dataset'], r['family'], r['level'], r['method'], r['case_clusters'], *[r['counts'].get(k, 0) for k in ['PASS', 'ABSTAIN', 'FAIL', 'UNKNOWN', 'UNKNOWN_SITE']], r['strict_pass_fraction'], *r['case_bootstrap_95'], r['resolution']])

    def con(ctrl, fam):
        return next((d for d in contrasts if d['control'] == ctrl and d['family'] == fam))
    lp = con('constraint_optimizer', 'all')
    rf = con(control, 'all')
    counts = summary[candidate]['counts']
    s_pass = counts.get('PASS', 0)
    lp_pass = summary['constraint_optimizer']['counts'].get('PASS', 0)
    text = f"""# Source-based crown generator tested against registered IOS - geometry\n\nThe generator combines a conditioned model of the entire tooth surface , a separate contact field and exact limitations of the locked height field. It has been tested once on 588 held out Bite2Text case, 21168 data, against the original published IOS coordinates. Models and all predictions was frozen prior to comparison.\n\nRun from this directory: `./run_all.sh` . run command the frozen generator on 36 public development tasks and exports nine research shell (STL ), controls frozen generator/model replacement and benchmark release, recalculates the same frozen predictions and recreates table and `figures/x42.png` . This is a replay of the same measurement; new test variants are rejected by the query-leader. The local, locked `../PROOF_LANE_GENCAD_V3` delivery with its Python runtime and lane data must remain. No data is downloaded.\n\n| Primary result, POPULATION with patient cases as cluster | Candidate minus control [95% cluster range] |\n|---|---|\n| Strict passed data | {s_pass}/21168; LP {lp_pass}/21168 |\n| Native height error against original LP , mm | {ci(lp['paired']['rmse_mm'])} |\n| Kontaktareafel against ursprunglig LP, mm² | {ci(lp['paired']['contact_error_mm2'])} |\n| Native height errors to direct standard kernel with the entire form field and the same contact field, mm | {ci(rf['paired']['rmse_mm'])} |\n| Contact area error against the same full-kernel, mm² | {ci(rf['paired']['contact_error_mm2'])} |\n\nNegative difference improves the error. Errors are compared on the same data where both designs get strict PASS ; loss and cluster number follow each difference in results.json . All data remains in the PASS denominator. Matched controls receives the same 208 development donors, all public fields and the same contact operation.\n\n| Familj, POPULATION; 588 patientfall vardera | Strikt PASS kandidat/LP | Δ height error against full-kernel, mm [ 95% ] | Δkontaktfel against full-kernel, mm² [95%] |\n|---|---:|---|---|\n"""
    official_ix = {(r['family'], r['participant']): r for r in official['summaries'] if r['split'] == 'test'}
    for fam in family:
        c = official_ix[fam, candidate]
        l = official_ix[fam, 'constraint_optimizer']
        d = con(control, fam)
        text += f"| {fam} | {c['counts'].get('PASS', 0)}/{l['counts'].get('PASS', 0)} | {ci(d['paired']['rmse_mm'], 4)} | {ci(d['paired']['contact_error_mm2'], 6)} |\n"
    text += "\nThe full gate requires more passed designs and better height/contact-Pareto than any strong control. It is counted as cases where no requirement is available. Result for all twelve participants, all nine families and easy/normal/hard/boundary are available in raw/FAMILY_LEVEL.csv and results.json . LP achieves the operating PASS ceiling in this release: 5234 data is impossible according to the facitore, 300 lacks location data and 16 lacks measurable antagonist overlap. This is an explanation within the scoring contract, no physical impossibility; see raw/PASS_CEILING_CERTIFICATE.json. The candidate's height error decreases against full-kernel in all nine families, but the contact error lacks certain improvement: four families have exactly identical contact area errors, implant has worse point estimates, and veneer has no contact details. The family intervals are descriptive and lack simultaneous multiplicity adjustment. Constant bootstrap is no guarantee of population safety; benchmark case-Hoeffding bands remain.\n\nWhat does not holds : benchmark level 1 applies only to the ceiling surface wall , nominal flat film, static non-penetration and the fixed deposit model. Tool radius, shaft, cervical margin and axial retention are not assessed. Level 3 produces native height errors and a static occlusal proximity area at fixed sparse assessment points, not approximate contacts, tooth removal or contact force. The measured scanner/bit registration uncertainty is UNKNOWN, so a decrease of approximately 0.006 mm in height errors is not evidence of the corresponding physical precision. The proximity area error assesses the total area and may miss the wrong contact position; a correct area number does not validate the local distribution of the contact field. A measured counterlateral mirror is missing in the public interface. Therefore, we cannot claim that all named baselines or a complete clinical crown design was hit.\n\nX1b runs with the explicit frozen V2 radii 0.25/0.5/0.8/1.0mm ; V3 does not have these requirements. R1 : s original X1b run error is preserved. X18 is its original QP and the same outer field is reused between the four levels where inputs is identical. The checks can reject incorrect values : twelve geometry -/schedule/source failure and eight contact tape injections are run with restored positive controls .\n\nData/ reference : [ Bite2Text ](https://ditto.ing.unimore.it/bite2text/) and [ Bits2Bites ](https://ditto.ing.unimore.it/bits2bites/), local files only. Bits2Bites specified CC BY - NC - SA in the assignment; exact license file is missing in V3 metadata. Bite2Text : s exact distribution conditions are UNKNOWN ; this is local research. No images, journals, personal names or new data sets are used. ToothFairy2 , mandible defects and Teeth3DS -mesher are not used; inherited X11-training specimen and unknown FDI-accuracy persists. Bits2Bites is reported separately and is never pooled with the primary test.\n\nThe field gate is located in FIELD_GENERATOR_PORT.md . CPU -export of nine research shell is executed (8 / 9 waterproof; the open grid shell is refused by the field port); the actual CPU field call has a separate fixture control; native/ GPU field call is NOT_RUN. A frozen design is a HANDOVER - state to manufacturing/measurement, no measured final product. Each individual geometric edge is preserved PER_POINT before the consumer aggregates.\n\nThe next design needs free preparation and a source supported cervical boundary in a separate benchmark version. Minimum external reference is a registered unprepared/prepared specimen with marked finish line, then a manufactured and seat-measured crown. Full cost and unknown legacy fit/discovery costs are in results.json . Status PENDING_INDEPENDENT_REVIEW.\n"
    (ROOT / 'README_DEMO.md').write_text(text)
    (ROOT / 'RESULTS.md').write_text(text.split('Run from this directory:')[0] + f"\nPrimary test: {s_pass}/21168 strikt PASS, LP {lp_pass}/21168. Original full X42 - gate REJECTED .\n\nAgainst LP : Δ height error {ci(lp['paired']['rmse_mm'])}mm; Δkontaktfel {ci(lp['paired']['contact_error_mm2'])}mm². Towards direct standard kernel with the entire shape field and the same contact field: {ci(rf['paired']['rmse_mm'])}mm and {ci(rf['paired']['contact_error_mm2'])}mm².\n\nSee README_DEMO.md , results.json and raw/FAMILY_LEVEL.csv for all families/levels, loss, cost and limitations. Only published native scan-geometry is reference; full manufacturing/approximate contact/removal is missing. PENDING_INDEPENDENT_REVIEW.\n")
    (fig, ax) = plt.subplots(1, 3, figsize=(15, 6))
    methods = [candidate, control, 'constraint_optimizer']
    color = ['#1b9e77', '#7570b3', '#666666']
    y = np.arange(len(family))
    for (j, (m, col)) in enumerate(zip(methods, color)):
        mu = [official_ix[f, m]['mean_case_pass_fraction'] for f in family]
        ax[0].plot(mu, y + j * 0.08, 'o', color=col, label=m)
    ax[0].set_yticks(y, family)
    ax[0].set_xlim(0, 1.03)
    ax[0].set_xlabel('Strict PASS fraction / patient-case')
    ax[0].legend(fontsize=7)
    ax[0].set_title('Fixed roof digital gates')
    for (axis, k, label) in [(ax[1], 'rmse_mm', 'Δ native-height RMSE (mm)'), (ax[2], 'contact_error_mm2', 'Δ proximity-area error (mm²)')]:
        for (j, (c, col)) in enumerate([(control, color[1]), ('constraint_optimizer', color[2])]):
            dd = [con(c, f)['paired'][k] for f in family]
            v = np.array([d['mean_difference'] if d else np.nan for d in dd])
            lo = np.array([d['case_bootstrap_95'][0] if d else np.nan for d in dd])
            hi = np.array([d['case_bootstrap_95'][1] if d else np.nan for d in dd])
            axis.errorbar(v, y + j * 0.12, xerr=np.array([v - lo, hi - v]), fmt='o', color=col, label='minus ' + c, markersize=4)
        axis.axvline(0, color='k', linewidth=0.6)
        axis.set_yticks(y, [])
        axis.set_xlabel(label)
        axis.set_title('Common PASS; negative improves')
        axis.legend(fontsize=7)
    fig.suptitle('X42 : 588 heldout cases; source reference , no complete restoration claim')
    fig.tight_layout()
    fig.savefig(ROOT / 'figures/x42.png', dpi=180)
    fig.savefig(ROOT / 'figures/x42.svg')
    plt.close(fig)
    (fig, axes) = plt.subplots(1, 2, figsize=(11, 5))
    for (axis, k, label) in [(axes[0], 'rmse_mm', 'Δ native-height RMSE (mm)'), (axes[1], 'contact_error_mm2', 'Δ proximity-area error (mm²)')]:
        dd = [con(control, f)['paired'][k] for f in family]
        v = np.array([d['mean_difference'] if d else np.nan for d in dd])
        lo = np.array([d['case_bootstrap_95'][0] if d else np.nan for d in dd])
        hi = np.array([d['case_bootstrap_95'][1] if d else np.nan for d in dd])
        axis.errorbar(v, y, xerr=np.array([v - lo, hi - v]), fmt='o', color=color[0], markersize=4)
        axis.axvline(0, color='k', linewidth=0.6)
        axis.set_yticks(y, family if k == 'rmse_mm' else [])
        axis.set_xlabel(label)
    fig.suptitle('Candidate minus same-data full-output kernel; 95% patient-case intervals')
    fig.tight_layout()
    fig.savefig(ROOT / 'figures/strongest_control.png', dpi=180)
    fig.savefig(ROOT / 'figures/strongest_control.svg')
    plt.close(fig)
    state('HIDDEN_EVALUATION_COMPLETE', {'parent_gate': 'FAIL', 'candidate_test_PASS': s_pass, 'LP_test_PASS': lp_pass}, 'Independent review; next version needs controllable preparation/native boundary', hidden_queries=1, result_sha256=sha(ROOT / 'results.json'))
