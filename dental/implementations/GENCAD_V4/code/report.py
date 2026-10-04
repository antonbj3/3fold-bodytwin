from common import *
from score_panel import boot
import collections

def fmt(x, d=3):
    return 'UNKNOWN' if x is None else f'{x:.{d}f}'

def stat(x):
    return 'UNKNOWN' if x['mean'] is None else f"{x['mean']:.3f} [{x['CI95'][0]:.3f}, {x['CI95'][1]:.3f}]"

def whole_summary(r):
    groups = collections.defaultdict(list)
    for row in r['rows']:
        groups[row['dataset'], row['split'], row['family'], row['participant']].append(row)
    out = []
    for ((ds, sp, f, p), rows) in sorted(groups.items()):
        good = [x for x in rows if x.get('full_surface_p95_mm') is not None]
        metrics = {'full_surface_p95_mm': boot([x['full_surface_p95_mm'] for x in good]), 'contact_symdiff_mm2': boot([x['contact']['contact_symdiff_mm2'] for x in good if x['contact']['contact_symdiff_mm2'] is not None])}
        for side in ['mesial', 'distal']:
            metrics[side + '_gap_error_mm'] = boot([x['proximal'][side]['mean_abs_gap_error_mm'] for x in good if side in x['proximal']])
        out.append(dict(dataset=ds, split=sp, family=f, participant=p, requested=len(rows), scored=len(good), closed=sum((x['digital_complete_shell'] for x in good)), anatomy_pass=sum((x['anatomy_gate'] for x in good)), quality=metrics, scope='One tooth per family per patient-case; bootstrap is patient-case clustered. Only4 primary test clusters; descriptive, broad and nonconfirmatory.'))
    return out

def semantic(x):
    if isinstance(x, dict):
        return {k: semantic(v) for (k, v) in x.items() if k not in ['seconds', 'cost_seconds', 'cost', 'costs', 'updated_utc', 'frozen_utc', 'wall_seconds', 'peak_rss_MiB']}
    if isinstance(x, list):
        return [semantic(v) for v in x]
    return x

def run():
    r = {k: read(ROOT / 'rounds' / f'{k}.json') for k in ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7']}
    panel = r['R2']
    full = whole_summary(r['R3'])
    prep = whole_summary(r['R5'])
    test = [x for x in panel['summaries'] if x['split'] == 'test']
    families = sorted({x['family'] for x in test})
    md = ['# DentalGenCAD v4 — spatial quality on registered scans', '', 'Primary test:48 patient-case clusters,9 families,4 difficulty tasks each. Eight actual participants;64 total cases including8 dev and8 Bits2Bites auxiliary. Dev is in sample for X42. Quality below is conditional on each participant passing L1, so use the paired contrasts for comparisons. These are roofs; complete-crown results follow separately.', '', '| Family | Participant | PASS / requested | PASS fraction [95% CI] | Shape RMSE mm [95% CI] | Contact mismatch mm² [95% CI] |', '|---|---|---:|---:|---:|---:|']
    for s in test:
        md.append(f"| {s['family']} | {s['participant']} | {s['counts'].get('PASS', 0)}/{s['requested']} | {stat(s['pass_fraction'])} | {stat(s['quality_on_pass']['anatomy_rmse_mm'])} | {stat(s['quality_on_pass']['contact_symdiff_mm2'])} |")
    md += ['', 'Paired test differences: positive values mean lower error for the second participant. Bootstrap resamples patient-case clusters after averaging paired difficulty rows within each case. No multiplicity correction; these intervals describe this retrospective panel. No single weighted quality score.', '', '| Family | Contrast | Metric | Paired patients / tasks | Mean difference [95% CI] |', '|---|---|---|---:|---:|']
    for p in panel['paired_contrasts']:
        if p['split'] == 'test' and p['a'] == 'constraint_optimizer' and (p['b'] in ['field_generator', 'x42_contact_branch']) and (p['metric'] != 'negative_gap_area_mm2'):
            md.append(f"| {p['family']} | {p['a']} − {p['b']} | {p['metric']} | {p['patients']}/{p['paired_tasks']} | {stat(p)} |")
    md += ['', 'Complete virtual crowns: actual X1B preparation, frozen roof lift, all generated exterior faces scored. Each primary family has4 patient-case clusters; these are not a clinical validation cohort. No complete clinical L1 is established.', '', '| Family | Participant | Scored / requested | Closed | Anatomy ≤0.35 mm | Full exterior p95 mm [95% CI] | Mesial gap error mm [95% CI] | Distal gap error mm [95% CI] |', '|---|---|---:|---:|---:|---:|---:|---:|']
    for s in full + prep:
        if s['split'] == 'test':
            md.append(f"| {s['family']} | {s['participant']} | {s['scored']}/{s['requested']} | {s['closed']} | {s['anatomy_pass']} | {stat(s['quality']['full_surface_p95_mm'])} | {stat(s['quality']['mesial_gap_error_mm'])} | {stat(s['quality']['distal_gap_error_mm'])} |")
    md += ['', '`prep_uniform` and `prep_regional` receive the complete prepared field and form a separate information track. All physical cement fields and generated-crown fracture margins remain UNKNOWN. Whole-preparation removal is a virtual volume, shared by participants; it does not discriminate their roof proposals. Every retained and missing row is in results.json and raw/QUALITY_ROWS.csv.gz.']
    (ROOT / 'LEADERBOARD.md').write_text('\n'.join(md) + '\n')
    original = read(ROOT / 'revisions/EXECUTION_ATTEMPT_1/rounds/R2.json')
    cost = dict(preparation=read(ROOT / 'FROZEN_WHOLE_INPUTS.json')['seconds'], participant_generation=read(ROOT / 'raw/GENERATION_PROCESS.json')['seconds'], participant_components=read(DATA / 'predictions/COST.json'), field_process=read(DATA / 'predictions/FIELD_COST.json'), whole_generation=read(ROOT / 'raw/WHOLE_GENERATION_PROCESS.json')['seconds'], prepared_field_generation=read(ROOT / 'raw/PREP_GENERATION_PROCESS.json')['seconds'], roof_scoring=panel['cost_seconds'], whole_scoring=r['R3']['seconds'], prepared_field_scoring=r['R5']['seconds'], preparation_replay=read(ROOT / 'raw/PREPARATION_REPLAY.json')['seconds'], failed_execution_roof_scoring=original['cost_seconds'], inherited_X11_fit='UNKNOWN', inherited_v3_fit='Reused195 training groups; historical cost not remeasured', inherited_X42_fit='Reused208 dev groups; historical fit recorded in source lane, not free', author_search_and_discovery='UNKNOWN', physical_measurement='NOT_RUN', scope='Metered components include repeated failed and corrected computation. They are not total research cost or a speedup claim.')
    res = dict(claim_type='capability', status='PARTIAL_CAPABILITY_WITH_NEGATIVE_FULL_CROWN_RESULT', review_state='PENDING_INDEPENDENT_REVIEW', headline='Reference-owned spatial quality distinguishes feasible generators; complete virtual shells expose a failed anatomy requirement.', external_referent=panel['external_referent'], additional_external_referents=[r['R4']['external_referent'], r['R6']['external_referent']], rounds=r, whole_family_summaries=full, prepared_field_family_summaries=prep, cost=cost, source_model_split=read(DATA / 'private/MODEL_SPLIT_PROOF.json'), input_integrity='Final frozen release checks all consumed inputs and runtime; first exploratory scorer code was not separately frozen before initial scores, though prereg/predictions were. Final replay is post hoc verification, not fresh blind evidence.', field_execution_correction='Attempt1 shared HiGHS global scheduler caused false abstentions. All attempts retained. Current run uses a separate Field process with byte-identical source.', resolution_map={'geometric_fields': 'PER_POINT', 'contacts_and_gap_regions': 'PER_SURFACE_REGION', 'whole_shell_and_volume': 'PER_TOOTH', 'patient_cluster_estimates': 'POPULATION', 'virtual_preparation_and_pose_scenario': 'PHENOMENOLOGICAL'}, debts=[{'quantity': 'virtual finish line / source insertion envelope', 'resolution': 'PHENOMENOLOGICAL', 'replacement_measurement': 'independent cervical annotation and registered original/prepared specimen'}, {'quantity': 'actual cement film', 'resolution': 'PER_POINT', 'replacement_measurement': 'spatial seated gap field with calibrated scanner or microCT and process metadata'}, {'quantity': 'fracture margin', 'resolution': 'PER_TOOTH', 'replacement_measurement': 'same crown geometry, material, support, cement, load location, angle/rate/aging with two calibration and held validation thicknesses'}, {'quantity': 'pose uncertainty', 'resolution': 'PHENOMENOLOGICAL', 'replacement_measurement': 'repeated registered loaded bite plus independent contact film'}], edge_contracts=read(ROOT / 'DECOMPOSITION.json')['edges'], data_manifest={'whole_inputs_sha256': sha(ROOT / 'FROZEN_WHOLE_INPUTS.json'), 'predictions_sha256': sha(ROOT / 'FROZEN_PREDICTIONS.json'), 'whole_predictions_sha256': sha(ROOT / 'FROZEN_WHOLE_PREDICTIONS.json'), 'prepared_field_predictions_sha256': sha(ROOT / 'FROZEN_PREP_PREDICTIONS.json')}, clinical_use=False)
    dump(ROOT / 'results.json', res)
    digest = hashlib.sha256(json.dumps(semantic(res), sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    dump(ROOT / 'raw/SEMANTIC_DIGEST.json', dict(sha256=digest, scope='Scientific values and decisions excluding timing/cost fields'))
    figure(r, test)
    return res

def figure(r, test):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 9})
    (fig, ax) = plt.subplots(1, 3, figsize=(16, 5), layout='constrained')
    ps = [p for p in r['R2']['paired_contrasts'] if p['split'] == 'test' and p['a'] == 'constraint_optimizer' and (p['b'] == 'x42_contact_branch') and (p['metric'] == 'contact_symdiff_mm2')]
    for (i, p) in enumerate(ps):
        if p['mean'] is not None:
            ax[0].errorbar(p['mean'], i, xerr=[[p['mean'] - p['CI95'][0]], [p['CI95'][1] - p['mean']]], fmt='o', color='#197a7d')
    ax[0].set_yticks(range(len(ps)), [p['family'] for p in ps])
    ax[0].axvline(0, color='grey', lw=1)
    ax[0].set_xlabel('LP − X42 contact mismatch (mm²)\npatient-cluster 95% bootstrap CI')
    ax[0].set_title('Same feasibility, spatial quality differs')
    names = ['constraint_optimizer', 'field_generator', 'x42_shape', 'x42_contact_branch']
    colors = ['#555555', '#7d5495', '#bb7229', '#197a7d']
    for (n, c) in zip(names, colors):
        s = next((s for s in test if s['family'] == 'molar_crown' and s['participant'] == n))
        a = s['quality_on_pass']['anatomy_rmse_mm']['mean']
        b = s['quality_on_pass']['contact_symdiff_mm2']['mean']
        ax[1].scatter(a, b, label=n, color=c, s=60)
    ax[1].set_xlabel('Molar shape RMSE (mm)')
    ax[1].set_ylabel('Contact symmetric difference (mm²)')
    ax[1].set_title('Quality vector; lower on each axis')
    ax[1].legend(fontsize=7)
    values = [[x['full_surface_p95_mm'] for x in r[k]['rows'] if x.get('full_surface_p95_mm') is not None and (p is None or x['participant'] == p)] for (k, p) in [('R3', None), ('R5', 'prep_uniform'), ('R5', 'prep_regional')]]
    ax[2].boxplot(values, tick_labels=['roof lift', 'prepared\nuniform', 'prepared\nregional'], showfliers=True)
    ax[2].axhline(0.35, color='crimson', ls='--', label='frozen 0.35 mm limit')
    ax[2].set_ylabel('Whole-exterior sampled p95 error (mm)')
    ax[2].set_title('All whole-crown anatomy gates fail')
    ax[2].legend(fontsize=7)
    fig.suptitle('DentalGenCAD v4 — registered IOS geometry; physical fit and strength UNKNOWN', fontsize=13)
    (ROOT / 'figures').mkdir(exist_ok=True)
    fig.savefig(ROOT / 'figures/quality_v4.png', dpi=180)
    fig.savefig(ROOT / 'figures/quality_v4.pdf')
    plt.close(fig)
if __name__ == '__main__':
    run()
