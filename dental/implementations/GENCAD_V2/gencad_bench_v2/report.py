import collections, subprocess, sys
from pathlib import Path
from .common import *

def run():
    from .lab_export import run as export_lab
    lab_pairs = export_lab()
    r = read(ROOT / 'rounds/R1.json')
    m = read(ROOT / 'raw/MEASUREMENT_SCORES.json')
    ctrl = read(ROOT / 'raw/CONTROL_REPORT.json')
    b = read(ROOT / 'data/BENCHMARK_LOCK.json')['payload']
    rows = read(ROOT / 'raw/SCORED_ROWS_R1.json')
    cost = {}
    for name in r['totals']:
        cost[name] = read(DATA / 'predictions' / name / '_cost.json')
    latest = read(ROOT / 'rounds/R2.json') if (ROOT / 'rounds/R2.json').exists() else None
    r3 = read(ROOT / 'rounds/R3.json') if (ROOT / 'rounds/R3.json').exists() else None
    result = dict(claim_type='capability', status='PARTIAL_BENCHMARK_EXECUTED_PENDING_REVIEW', capability='Generative surface completion and explicit abstention measured against withheld native IOS in the paired antagonist frame', requested_tasks=b['task_count'], source_ready_tasks=b['ready'], unknown_source_tasks=b['unknown'], source_case_records=b['source_records'], split_case_counts=b['case_counts'], external_referent=r['external_referent'], R1=r, R2=latest, external_measurement_tracks=m, controls=ctrl, physical_generated_part_validation='UNKNOWN_NO_MATCHED_MANUFACTURE_MEASUREMENT_PAIR', whole_restoration_geometry='INCOMPLETE_AXIAL_WALLS_MARGIN_RETENTION_EMERGENCE_IMPLANT_CONNECTION', patient_disjointness='Case records and exact source hashes disjoint; latent cross-dataset identity UNKNOWN', hidden_test='Held out from process-isolated plugins; author-visible retrospective references; not a prospectively secret public competition', cost=dict(preparation_seconds=b['preparation_seconds'], fit_seconds_included_in_preparation=b['fit_seconds'], generators={k: {f: v[f] for f in ['seconds', 'max_rss_kib']} for (k, v) in cost.items()}, scoring_seconds=r['seconds'], measurement_fit=read(ROOT / 'raw/measurement_fit_cost.json'), discovery='Implementation and design effort: UNKNOWN wall-clock attribution; no free discovery claim', inherited_data_acquisition_and_training='UNKNOWN; X11 trained elsewhere', fallback='Per-task exceptions, unknown sites and abstentions retained', lane_data_bytes=budget(), thread_limit=4, gpu=False), artifacts=dict(data_root=str(DATA), benchmark_lock_sha256=sha(ROOT / 'data/BENCHMARK_LOCK.json'), frozen_predictions_sha256=sha(ROOT / 'FROZEN_PREDICTIONS_R1.json'), figure='figures/benchmark_v2.png', rows='raw/SCORED_ROWS_R1.json'), review_state='PENDING_INDEPENDENT_REVIEW')
    result['R3'] = r3
    result['frozen_lab_geometry_pairs'] = dict(count=len(lab_pairs), manifest='FROZEN_PREDICTIONS_FOR_LAB.json', sha256=sha(ROOT / 'FROZEN_PREDICTIONS_FOR_LAB.json'), physical_prediction='UNKNOWN; digital geometry predicates only')
    if (ROOT / 'raw/RECOMPUTATION_RECEIPT.json').exists():
        result['fresh_generator_recomputation'] = read(ROOT / 'raw/RECOMPUTATION_RECEIPT.json')
    if (ROOT / 'raw/EXTERNAL_PLUGIN_SMOKE.json').exists():
        result['external_plugin_smoke'] = read(ROOT / 'raw/EXTERNAL_PLUGIN_SMOKE.json')
    if (ROOT / 'raw/INTEGRATION_CONTROLS.json').exists():
        result['additional_integration_controls'] = read(ROOT / 'raw/INTEGRATION_CONTROLS.json')
    if (ROOT / 'raw/MATERIAL_REPLAY_RECEIPT.json').exists():
        result['inherited_material_witnesses'] = read(ROOT / 'raw/MATERIAL_REPLAY_RECEIPT.json')
    result['whole_restoration_geometry'] = 'R1 surface tasks incomplete. R3 has 31 complete digital crown shells with conical cavities; clinical margin/retention/implant connection and other full restoration classes remain unvalidated.'
    dump(ROOT / 'results.json', result)
    line = ['# DentalGenCAD-Bench v2 — per-family and difficulty outcomes', '', 'These are partial surface-restoration tasks. PASS refers only to the stated digital model. Complete clinical restoration and material/process calibration remain UNKNOWN. Counts retain unavailable sites and abstentions. Natural teeth are morphology references, not optimal restorations.', '', '| Participant | Digital PASS | FAIL | UNKNOWN | ABSTAIN | INVALID |', '|---|---:|---:|---:|---:|---:|']
    for (name, c) in r['totals'].items():
        line.append(f"| {name} | {c.get('PASS', 0)} | {c.get('FAIL', 0)} | {c.get('UNKNOWN', 0) + c.get('UNKNOWN_SITE', 0)} | {c.get('ABSTAIN', 0)} | {c.get('INVALID', 0)} |")
    line += ['', 'No scalar leaderboard score. Pareto fronts use height RMSE, contact-area error, nominal film error and volume only for scoped PASS designs with measured overlap.', '', '| Split | Family | Level | Participant | N | PASS/FAIL/UNKNOWN/ABSTAIN/INVALID | Median height RMSE mm | Median contact-area error mm² | Pareto count |', '|---|---|---|---|---:|---|---:|---:|---:|']
    for a in r['summaries']:
        c = a['counts']
        ct = '/'.join((str(c.get(k, 0) + (c.get('UNKNOWN_SITE', 0) if k == 'UNKNOWN' else 0)) for k in ['PASS', 'FAIL', 'UNKNOWN', 'ABSTAIN', 'INVALID']))
        fmt = lambda v: '—' if v is None else f'{v:.4f}'
        line.append(f"| {a['split']} | {a['family']} | {a['level']} | {a['participant']} | {a['n']} | {ct} | {fmt(a['height_rmse_median_mm'])} | {fmt(a['contact_area_error_median_mm2'])} | {a['pareto_count']} |")
    if r3:
        line += ['', '## Full axial-crown extension (R3)', '', f"32 requested normal molar/premolar tasks; {r3['designed_meshes']} meshes = 31 paired designs plus one retained unavailable roof. Complete-cap wall: unflared 31 FAIL; axial flare 31 PASS. All 31 gains retain minimum/maximum film and measured-antagonist nonpenetration. Cervical rim excluded from wall. No clinical/physical admission.", '', 'R2 is preserved: both arms passed 31/31 because the cavity was overly conservative; the intended new-pass gate failed. R3 corrects the finite-boundary-segment operation without changing any material or film threshold.']
    (ROOT / 'LEADERBOARD.md').write_text('\n'.join(line) + '\n')
    text = ['# What can now be measured', '', f"The locked benchmark executes {b['task_count']} requested tasks over nine restoration families and four difficulty levels; {b['ready']} have usable source geometry and {b['unknown']} remain explicit UNKNOWN. Five participants run without access to held-out reference meshes. The independent referent is native geometry from Bits2Bites/Bite2Text, with the antagonist pose preserved.", '', f"Scoped surface-model feasibility: `{r['feasibility_counts']}`. Complete crown/bridge construction is not certified: the present representation lacks clinical axial walls, cervical margins and implant connections. Published-rule cards do not make the horizontal preparation patches clinically faithful preparations.", '', '| Participant | PASS | FAIL | ABSTAIN | Other |', '|---|---:|---:|---:|---:|']
    for (name, c) in r['totals'].items():
        text.append(f"| {name} | {c.get('PASS', 0)} | {c.get('FAIL', 0)} | {c.get('ABSTAIN', 0)} | {sum(c.values()) - sum((c.get(k, 0) for k in ['PASS', 'FAIL', 'ABSTAIN']))} |")
    text += ['', 'External measurement tracks retain their negative results: fracture transport fails the 25% all-group error gate; the cement regional transfer gate is ' + str(m['cement']['gate_all_primary_MAE_le_20um']) + '; the sintering angle RMSE is ' + str(m['sinter']['angles']['candidate']['rmse']) + ' degrees. These are literature prediction instruments, not generated-part strength/fit scores.', '', f"{m['cement']['cells']} cement cells, {m['fracture']['groups']} fracture groups, {m['sinter']['angle_cells']} angle and {m['sinter']['shrinkage_cells']} shrinkage cells are scored with source locators. Pulpy/TF2 tissue intervals use {m['tissue']['untruncated_rows']} untruncated rows in their own CT cohort. STS has no independently verified pulp mask in the reused audit; IOS-to-pulp safety is UNKNOWN.", '', 'Controls: the patched exact predicates reject all four previously published dimension exploits; positive fixtures and injected failures pass the executable test suite. The plugin sandbox denies a real private-file read. Equal-information LP/QP results are numerical controls, not a method-superiority claim.', '', f'Data/intermediates: {budget() / 1000000.0:.1f} MB. See results.json for all runtime and memory costs, LEADERBOARD.md for every family/level, and raw/SCORED_ROWS_R1.json for retained individual failures.', '', 'Status: PENDING_INDEPENDENT_REVIEW. No clinical recommendation, manufacture, measured force or physical prospective validation.']
    if r3:
        text += ['', '## A roof pass is insufficient for a full crown', '', f"R3 constructs actual axial walls, conical intaglio and annular margins for 31 of 32 requested crown sites. All 31 unflared complete caps fail the exact material wall requirement, despite their R1 roofs passing. The changed axial flare gives {len(r3['new_wall_pass_without_film_or_occlusal_failure'])} full-cap wall passes with no minimum/maximum-film or measured-antagonist regression. These are rational checks on serialized geometry and independent native-antagonist comparisons. No clinical finish-line accuracy or physical fit is established.", '', 'R2 reached no new passes because an infinite-line inradius made the cavity too small; its negative gate and meshes are preserved. R3 uses finite boundary segments and an independently checked scalar minimizer. Same cohort, retrospective construction refinement; no fresh test-patient generalization.']
    (ROOT / 'RESULTS.md').write_text('\n'.join(text) + '\n')
    fig = dict(totals=r['totals'], summaries=r['summaries'], feasibility=r['feasibility_counts'], measurements=m['cement']['summary'], example=None)
    for row in rows:
        if row['family'] == 'molar_crown' and row['level'] == 'normal' and (row['participant'] == 'constraint_optimizer') and (row['checks']['validity'] == 'PASS'):
            from .tasks import load_task
            t = load_task(DATA / 'public' / (row['task_id'] + '.json'))
            d = read(DATA / 'predictions/constraint_optimizer' / (row['task_id'] + '.json'))
            fig['example'] = dict(xy=t['xy'], faces=t['faces'], outer=d['outer_vertices'], inner=d['inner_vertices'], antagonist=t['antagonist'], task_id=row['task_id'])
            break
    dump(ROOT / 'raw/FIGURE_DATA.json', fig)
    if r3:
        tid = r3['new_wall_pass_without_film_or_occlusal_failure'][0]
        examples = {}
        for arm in ['unflared', 'axial_flare']:
            a = np.load(DATA / 'axial_R3' / (tid + '_' + arm + '.npz'))
            examples[arm] = {k: a[k] for k in ['vertices', 'faces', 'inner_vertices', 'inner_faces']}
        dump(ROOT / 'raw/AXIAL_FIGURE_DATA.json', dict(task_id=tid, arms=examples))
    p = subprocess.run(['/usr/bin/python3', '-s', str(ROOT / 'gencad_bench_v2/plot.py')], cwd=ROOT)
    if p.returncode:
        raise RuntimeError('figure rendering failed')
    state('THREE_ROUNDS_REPORTED', 'R1 576 tasks; R2 negative; R3 31 complete-cap wall gains', 'Independent review and native clinical preparation/finish-line measurement')
    return result
