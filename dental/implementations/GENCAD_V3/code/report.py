"""Researcher-facing outputs; only aggregated heldout results are public."""
from pathlib import Path
import numpy as np
from util import *

def generate(scores, controls, source, r2, r3, r4, r5, r6, r7, bundle):
    from integrity import safe_bytes, h, strict_json, IntegrityError
    from release_anchor import PREDICTIONS_SHA256
    pred_blob = safe_bytes(PAYLOAD / 'predictions', 'FROZEN_PREDICTIONS.json', 2000000)
    if h(pred_blob) != PREDICTIONS_SHA256:
        raise IntegrityError('prediction receipt changed before report')
    pred = strict_json(pred_blob)
    cost_blob = safe_bytes(ROOT, 'raw/GENERATION_PROCESS.json', 2000000)
    if h(cost_blob) != pred['controller_cost_sha256']:
        raise IntegrityError('generator cost receipt changed')
    generator_cost = strict_json(cost_blob)
    primary = [r for r in scores['summaries'] if r['dataset'] == 'Bite2Text' and r['split'] == 'test']
    cases = primary[0]['case_clusters']
    count = scores['task_count']
    cohort = bundle.json('raw/COHORT_SUMMARY.json')
    lines = ['# DentalGenCAD-Bench v3', '', f"{count:,} surface tasks, {scores['case_count']} evaluation case records; {cases} heldout Bite2Text patient-case clusters. These are scoped digital checks, not clinical success rates. Each case has four levels per family. UNKNOWN_SITE counts zero in strict PASS fraction; abstentions are separate.", '', 'All intervals below are 95% case-cluster bootstrap intervals, POPULATION resolution. They exclude scanner, FDI and domain uncertainty. The bounded-case Hoeffding intervals and task ICC sensitivity are in `results.json`.', '', '| Family | Participant | PASS fraction [95% CI] | Cases | Tasks | Estimated task n_eff | RMSE median (mm) | Correct abstentions |', '|---|---|---:|---:|---:|---:|---:|---:|']
    for r in primary:
        (lo, hi) = r['case_cluster_bootstrap_95']
        ne = r['n_eff_tasks_estimated']
        rm = r['height_rmse_median_mm']
        ne_text = f'{ne:.1f}' if ne is not None else 'UNKNOWN constant outcomes'
        rm_text = f'{rm:.3f}' if rm is not None else 'UNKNOWN'
        lines.append(f"| {r['family']} | {r['participant']} | {r['mean_case_pass_fraction']:.3f} [{lo:.3f}, {hi:.3f}] | {r['case_clusters']} | {r['task_rows']} | {ne_text} | {rm_text} | {r['correct_abstentions']} |")
    lines += ['', 'Task n_eff is a method-of-moments ICC diagnostic. The number of independent patient-case units is always the Cases column, conditional on upstream case/patient documentation. A degenerate bootstrap is flagged and accompanied by a nondegenerate bounded-case interval in the JSON.', '', f"Anatomical volume track: original small-domain preparation {r2['geometric_pass']}/{r2['requested']}; extended-domain preparation {r3['geometric_pass']}/{r3['requested']}; height-envelope research crowns {r4['geometric_pass']}/{r4['requested']}; source-topology crowns {r6['geometric_pass']}/{r6['requested']}. Failed requests remain in denominators. R7 source-support eligibility: {r7['exported']}/{r7['requested']}. These retrospective geometry tracks are separate from the hidden surface test.", '', 'Control: the inherited equal-information constrained LP is executed; trivial optimality certificates were compared against actual LP solves. No algorithm superiority is claimed.']
    for (ds, split) in [('Bite2Text', 'dev'), ('Bits2Bites', 'auxiliary')]:
        lines += ['', f'## {ds}: {split}', '', '| Family | Participant | Strict PASS fraction [95% CI] | Case records |', '|---|---|---:|---:|']
        for r in scores['summaries']:
            if r['dataset'] == ds and r['split'] == split:
                (lo, hi) = r['case_cluster_bootstrap_95']
                lines.append(f"| {r['family']} | {r['participant']} | {r['mean_case_pass_fraction']:.3f} [{lo:.3f}, {hi:.3f}] | {r['case_clusters']} |")
    (ROOT / 'LEADERBOARD.md').write_text('\n'.join(lines) + '\n')
    out = dict(claim_type='capability', review_state='PENDING_INDEPENDENT_REVIEW', capability='Integrity-checked portable evaluator on a larger within-dataset patient-case split, with a full-volume source/preparation/crown track.', external_referent=scores['external_referent'], additional_external_referent=dict(kind='external_review', locator='review_evidence/XREVIEW.md', compared_quantity='C01 copied-run failure and C02 clearance edit accepting FAIL as PASS; reproduced as regression controls', refutes_us=True), resolutions=dict(source_geometry='PER_POINT', restoration_metrics='PER_TOOTH', exterior_percentiles_and_source_coverage='PER_SURFACE_REGION', summary_statistics='POPULATION', virtual_preparation_parameters='PHENOMENOLOGICAL'), benchmark_sha256=scores['benchmark_sha256'], prediction_freeze_sha256=scores['prediction_freeze_sha256'], cohort=cohort, surface_benchmark=scores, integrity=controls, source_reference_audit={k: v for (k, v) in source.items() if k != 'rows'}, anatomical_track=dict(R2={k: v for (k, v) in r2.items() if k != 'rows'}, R3={k: v for (k, v) in r3.items() if k != 'rows'}, R4={k: v for (k, v) in r4.items() if k != 'rows'}, R5={k: v for (k, v) in r5.items() if k != 'rows'}, R6={k: v for (k, v) in r6.items() if k != 'rows'}, R7={k: v for (k, v) in r7.items() if k != 'rows'}), artifact_manifests=dict(benchmark='BENCHMARK_LOCK.json', predictions='FROZEN_PREDICTIONS.json', anatomical=['FROZEN_PREPARATIONS.json', 'FROZEN_PREPARATIONS_R3.json', 'FROZEN_CROWNS.json', 'FROZEN_CROWNS_R5.json', 'FROZEN_CROWNS_R6.json', 'FROZEN_CROWNS_R7.json']), protocol_limits=['R3–R7 are retrospective refinements', 'Expanded R3–R6 first-principles correspondence written retrospectively as DECOMPOSITION_MAP.json; original mechanisms/equations/closures remain in the pre-construction PREREG files and initial decomposition'], hidden_test=dict(within_Bite2Text_patient_case_disjoint=True, patient_case_clusters=cases, participant_reference_isolation=True, author_blind=False, global_cross_dataset_patient_disjointness='UNKNOWN', source_archives_available_to_participant=False, preparation_track_retrospective=True, query_policy='One frozen baseline release. New submissions require referee intake receipt; no test reference or per-case feedback in public leaderboard.'), limitations=['No observed preparation or clinical finish line', 'X11 target-domain FDI accuracy UNKNOWN', 'No setup-matched manufactured/seated-gap/force measurements', '0.25mm preparation field does not validate an0.08mm physical cement film', 'Bits2Bites not pooled with Bite2Text as independent patients', 'Linux x86_64 ABI and unprivileged namespaces required; other OS/architectures not tested', 'Trusted evaluator and externally retained release digest required; owner/root replacing both is out of scope', 'V2 published fracture/cement transfer failures are not overturned; those physical tracks are not rerun in v3'], debts=bundle.json('PREREG_R7.json')['debts'], edges=[dict(producer='published STL source triangle', consumer='reference-height scorer', resolution='PER_POINT', time_scale='SIMULTANEOUS', binding='private/audit/PROVENANCE.json and immutable scene/reference hashes'), dict(producer='X18b source surface/envelope', consumer='X1B preparation', resolution='PER_POINT', time_scale='HANDOVER', status='virtual preparation closure'), dict(producer='R3 preparation field', consumer='R4 crown cavity/mesh generator', resolution='PER_POINT', time_scale='HANDOVER', status='executed; physical fit UNKNOWN')], cost=dict(build=bundle.json('raw/BUILD_COST.json'), generator=generator_cost, reused_X11_training_and_prior_research='UNKNOWN, not charged as free', developer_discovery_time='UNKNOWN', offline_replay='raw/LAST_RUN.json', query_count='One frozen baseline execution; controls and replays logged separately', fallback='None: invalid inputs, failed geometry and missing patient links remain rejected/UNKNOWN'))
    out['cost']['primary_template_fit_seconds'] = bundle.json('payload/private/TEMPLATE_FIT.json')['seconds']
    out['cost']['anatomical_construction_record_times'] = {}
    for fn in out['artifact_manifests']['anatomical']:
        records = bundle.json(fn)['payload']['records']
        observed = [r['seconds'] for r in records if 'seconds' in r]
        out['cost']['anatomical_construction_record_times'][fn] = dict(sum_observed_seconds=sum(observed) if observed else None, records_with_timing=len(observed), records_total=len(records), missing_times='UNKNOWN, not zero')
    out['cost']['accounting_limits'] = 'BUILD_COST is the final resumed cohort segment, not total discovery/build cost. Earlier interrupted segments, Python/package preparation, development and prior X11 fitting are not fully metered; retained command/log records are in COMMANDS.md and raw. R6 interrupted resource probe is separate from its completed construction. No missing cost is treated as free.'
    relocation = ROOT / 'raw/RELOCATION_RECEIPT.json'
    if relocation.exists():
        out['relocation'] = read(relocation)
    relocation_status = out.get('relocation', {}).get('status', 'PENDING_UNTIL_RELOCATION_RECEIPT')
    dump(ROOT / 'results.json', out)
    (ROOT / 'RESULTS.md').write_text(f"# A portable evaluator with case-level uncertainty\n\nDentalGenCAD-Bench v3 evaluates {count:,} locked surface tasks from {scores['case_count']} evaluation records against native published IOS geometry. The primary hidden process test contains {cases} Bite2Text patient-case groups; its family-level estimates and case bootstrap intervals are in LEADERBOARD.md. Bits2Bites is reported separately. Source rows excluded: {len(cohort['excluded'])}/{cohort['attempted']}; site-level exclusions and exact reasons are retained in results.json.\n\nIntegrity: {controls['integrity']['count']} injected post-freeze/schema attacks rejected, with restored positive controls. The original C02 witness changes FAIL to PASS in the unguarded geometric check and is refused by the supported v3 score and submission entry points. C01 relocation status: {relocation_status}. The relocation test moves the entire data/runtime payload into a standalone delivery and runs with original workspaces, source archives and network inaccessible; see raw/RELOCATION_RECEIPT.json.\n\nSource-linked full-volume preparation: R2 failed {r2['failed']}/{r2['requested']}; R3 passes {r3['geometric_pass']}/{r3['requested']}, retaining its two bounded-grid exclusions. R4 consumes those fields and constructs research crown meshes: {r4['geometric_pass']}/{r4['requested']} pass the frozen combined digital gates. Native exterior identity, mesh topology and source-distance checks are separate from physical cement fit or strength. Every claimed check has a concrete rejecting fault. R5 exports {r5['exported']}/{r5['requested']} because the source boundary branches. R6 splits disconnected vertex fans without moving observed triangles: {r6['exported']}/{r6['requested']} exported and {r6['geometric_pass']}/{r6['requested']} pass all frozen geometric gates. These eight retrospective cases are not an independent generalization test. R7 requires a clipping margin above every virtual closure and50% retained native area: {r7['exported']}/{r7['requested']} are eligible. Refusal is a scoped answer, not successful crown generation.\n\nThese are capability claims. The constrained LP is an equally informed executed control, with no claim of algorithm superiority. Every test-family estimate has a case-cluster interval; task ICC estimates are diagnostics and never counts of independent patients. Global patient identity across datasets/pretraining remains UNKNOWN.\n\nAll thresholds remain in the original PREREG files. R2's negative result and its exact predictions are retained. Preparation/crown refinements reuse the same24 requests and therefore carry no fresh heldout generalization claim. The primary surface algorithms use only the frozen public inputs and primary training templates.\n\nExternal facit: published Bite2Text/Bits2Bites native STL coordinates, with source member hashes and independently enumerated triangle witnesses. Physical prepared anatomy, finish lines, seated gap, fabrication error and force remain unmeasured. PENDING_INDEPENDENT_REVIEW.\n")
    return out

def figure(scores, r3, r4, r6, bundle):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'svg.hashsalt': 'PROOF_LANE-gencad-v3'})
    fig = plt.figure(figsize=(13.5, 8.5))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.2, 1.0])
    ax = fig.add_subplot(gs[0, :])
    names = ['population', 'parametric', 'constraint_optimizer']
    colors = ['#8a8a8a', '#4375b9', '#c85b30']
    for (k, (name, color)) in enumerate(zip(names, colors)):
        rs = [next((r for r in scores['summaries'] if r['dataset'] == 'Bite2Text' and r['split'] == 'test' and (r['family'] == f) and (r['participant'] == name))) for f in FAMILIES]
        mean = np.array([r['mean_case_pass_fraction'] for r in rs])
        ci = np.array([r['case_cluster_bootstrap_95'] for r in rs])
        x = np.arange(9) + (k - 1) * 0.19
        ax.errorbar(x, mean, yerr=[mean - ci[:, 0], ci[:, 1] - mean], fmt='o', capsize=3, color=color, label=name, ms=4)
    ax.set(xticks=np.arange(9), xticklabels=[f.replace('_', '\n') for f in FAMILIES], ylim=(-0.03, 1.04), ylabel='Strict scoped PASS fraction', title='Patient-case heldout test: 95% case bootstrap intervals (population resolution)')
    ax.grid(axis='y', alpha=0.25)
    ax.legend(loc='upper left', bbox_to_anchor=(0, 1.14), ncol=3, frameon=False)
    choices = [r['case_key'] for r in r6['rows'] if r['split'] == 'dev' and r['status'] == 'EXPORTED']
    prep_folder = 'preparations_r6'
    crown_folder = 'crowns_r6'
    if not choices:
        choices = [r['case_key'] for r in r4['rows'] if r['split'] == 'dev' and r['status'] == 'EXPORTED']
        prep_folder = 'preparations_r3'
        crown_folder = 'crowns'
    if choices:
        key = choices[0]
        a = bundle.npz('payload/' + prep_folder + '/' + key + '/fields.npz')
        b = bundle.npz('payload/' + crown_folder + '/' + key + '/prediction.npz')
        meshes = [(a['source_vertices'], a['source_faces'], 'Measured source exterior', '#aaaaaa'), (a['preparation_vertices'] @ a['R'].T + a['center'], a['preparation_faces'], 'Virtual X1B preparation', '#548ab5'), (b['vertices'] @ a['R'].T + a['center'], b['faces'], 'R6 crown: source-distance FAIL', '#cf9e52')]
        for (j, (v, f, title, color)) in enumerate(meshes):
            ax = fig.add_subplot(gs[1, j], projection='3d')
            idx = np.linspace(0, len(f) - 1, min(8000, len(f)), dtype=int)
            vv = v - a['center']
            pc = Poly3DCollection(vv[f[idx]], facecolor=color, edgecolor='none', alpha=0.97)
            ax.add_collection3d(pc)
            lo = vv.min(0)
            hi = vv.max(0)
            mid = (lo + hi) / 2
            rad = max(hi - lo) / 2
            ax.set(xlim=(mid[0] - rad, mid[0] + rad), ylim=(mid[1] - rad, mid[1] + rad), zlim=(mid[2] - rad, mid[2] + rad), title=title, xlabel='mm', ylabel='mm')
            ax.view_init(26, -62)
            ax.set_box_aspect((1, 1, 1))
    fig.text(0.02, 0.015, f"Volume track: {r3['geometric_pass']}/{r3['requested']} closed virtual preparations; {r6['geometric_pass']}/{r6['requested']} source-topology crowns pass. Physical fit and strength UNKNOWN.", fontsize=9)
    fig.subplots_adjust(left=0.07, right=0.99, top=0.9, bottom=0.08, hspace=0.4, wspace=0.1)
    out = ROOT / 'figures'
    out.mkdir(exist_ok=True)
    fig.savefig(out / 'benchmark_v3.png', dpi=160)
    fig.savefig(out / 'benchmark_v3.svg', metadata={'Date': None})
    plt.close(fig)
