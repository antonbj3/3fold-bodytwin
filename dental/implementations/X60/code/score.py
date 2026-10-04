from common import *
import sys, importlib.util, time, collections
sys.path.insert(0, str(V4 / 'payload/participant_code'))
from task_io import attach
from legacy.checks import all_checks
from legacy.generators import submit
from legacy.geometry import shell, export_stl
from mechanics import roof_model, gradient_check
spec = importlib.util.spec_from_file_location('reference_quality', V4 / 'code/quality.py')
q = importlib.util.module_from_spec(spec)
spec.loader.exec_module(q)

def verified_npz(path, expected):
    if sha(path) != expected:
        raise ValueError('NPZ hash ' + str(path))
    with np.load(path, allow_pickle=False) as a:
        return dict(a)

def run():
    start = time.perf_counter()
    panel = read(ROOT / 'PANEL_SELECTION.json')
    lock = read(ROOT / 'INPUT_LOCK.json')['files']
    v4freeze = read(V4 / 'FROZEN_PREDICTIONS.json')
    rows = []
    checks = []
    fem = []
    preds = {}
    diags = {}
    names = ['population', 'parametric', 'X1B_morphology', 'X18_antagonist', 'constraint_optimizer', 'field_generator', 'x42_shape', 'x42_contact_branch']
    qref = next((r for r in read(ROOT / 'raw/R1.json')['published_force05'] if r['thickness_mm'] == 0.8))['force05_N']
    for tag in ['volume_1.0', 'volume_1.1']:
        pf = read(ROOT / f'FROZEN_PREDICTIONS_{tag}.json')
        preds[tag] = pf['files']
        diags[tag] = {r['task_id']: r for r in read(DATA / 'predictions' / tag / 'DIAGNOSTICS.json')}
    for case in panel['selected']:
        key = case['case_key']
        tasks = read(V4 / 'payload/public/tasks' / f'{key}.json')
        target = verified_npz(V4 / 'payload/private/references' / f'{key}.npz', lock[str(V4 / 'payload/private/references' / f'{key}.npz')]['sha256'])
        scenes = {}
        newidx = 0
        packs = {}
        for name in names:
            path = V4 / 'payload/predictions' / name / f'{key}.npz'
            packs[name] = verified_npz(path, v4freeze['files'][name + '/' + path.name]['sha256'])
        for tag in preds:
            packs[tag] = verified_npz(DATA / 'predictions' / tag / f'{key}.npz', preds[tag][key + '.npz']['sha256'])
        for (oldidx, t0) in enumerate(tasks):
            if t0['family'] not in panel['families']:
                continue
            t = t0
            if t['status'] == 'READY':
                rel = t['geometry_file']
                path = V4 / 'payload/public' / rel
                if rel not in scenes:
                    scenes[rel] = verified_npz(path, lock[str(path)]['sha256'])
                t = attach(t, scenes[rel])
            taskrows = []
            for name in names + list(preds):
                a = packs[name]
                i = newidx if name in preds else oldidx
                if str(a['task_ids'][i]) != t['task_id']:
                    raise ValueError('wrong prediction task')
                st = str(a['status'][i])
                row = dict(task_id=t['task_id'], case_key=key, split=case['split'], dataset=case['dataset'], family=t['family'], level=t['level'], participant=name, status=st, resolution='PER_TOOTH', calibrated_fracture_force05_N=None, physical_substance_removal_mm3=None, physical_Pareto='UNKNOWN_MISSING_CALIBRATED_FORCE_AND_PHYSICAL_REMOVAL', time_scale='SIMULTANEOUS')
                if t['status'] != 'READY':
                    row.update(L1='UNKNOWN_SITE', reason=t['site_error'])
                elif st != 'DESIGN':
                    row.update(L1=st, reason=diags.get(name, {}).get(t['task_id'], {}).get('reason'))
                else:
                    z = a['outer_' + str(i)]
                    inner = a['inner_' + str(i)]
                    d = submit(t, z, inner)
                    ck = all_checks(t, d)
                    row.update(L1=ck['validity'], checks=ck['checks'], complete_restoration_eligibility=ck['complete_restoration_eligibility'])
                    row.update(q.roof_metrics(t, z, inner, target[t['family']]))
                    row['digital_roof_volume_mm3'] = float(np.asarray(t['weights']) @ (z - inner))
                    good = np.isfinite(target[t['family']])
                    row['virtual_roof_removal_quadrature_mm3'] = float(np.asarray(t['weights'])[good] @ np.maximum(target[t['family']][good] - t['preparation_z'][good], 0))
                    row['removal_scope'] = 'Area-weighted virtual original roof minus fixed horizontal preparation; not retained dentin or full tooth volume. Fixed across participants.'
                    if name in preds:
                        dd = diags[name][t['task_id']]
                        row['optimization'] = dd
                        row['conditional_model_force05_N'] = qref * dd['model_force_quantile_ratio']
                        row['conditional_model_force05_scope'] = 'Borrowed Prott G0.8 scale normalized to this task morphology baseline; unvalidated force-law/volume-flaw transfer, NOT calibrated absolute force'
                        row['meaningful_conditional_gain'] = bool(ck['validity'] == 'PASS' and dd['model_force_quantile_ratio'] > 1.02 and (dd['equilibrium_relative_residual'] <= 1e-08))
                        if len(fem) < 4:
                            basez = packs['constraint_optimizer']['outer_' + str(oldidx)]
                            bi = packs['constraint_optimizer']['inner_' + str(oldidx)]
                            (model, M) = roof_model(t, basez, bi)
                            gc = gradient_check(model, np.asarray(dd['coefficients_mm']))
                            gc.update(task_id=t['task_id'])
                            fem.append(gc)
                    if len(checks) < 16 and name == 'constraint_optimizer' and (ck['validity'] == 'PASS'):
                        faults = []
                        faults.append(('wall', submit(t, inner + float(t['requirements']['wall_mm']) / 2, inner)))
                        faults.append(('film_min', submit(t, z, t['preparation_z'] - 1)))
                        faults.append(('film_max', submit(t, z, t['preparation_z'] + 1)))
                        faults.append(('nesting', submit(t, inner - 1, inner)))
                        if t['A'].shape[0]:
                            faults.append(('antagonist', submit(t, z + 100, inner)))
                        if t['family'] == 'bridge3':
                            faults.append(('connector_area', submit(t, inner + 0.001, inner)))
                        for (field, bad) in faults:
                            r = all_checks(t, bad)
                            checks.append(dict(task_id=t['task_id'], field=field, rejected=r['validity'] == 'FAIL', check=r['checks'].get(field)))
                taskrows.append(row)
                rows.append(row)
            valid = [r for r in taskrows if r.get('L1') == 'PASS' and all((r.get(k) is not None for k in ['anatomy_rmse_mm', 'contact_symdiff_mm2', 'negative_gap_area_mm2']))]
            if valid:
                mask = q.nondominated([[r[k] for k in ['anatomy_rmse_mm', 'contact_symdiff_mm2', 'negative_gap_area_mm2']] for r in valid], [0.01, 0.05, 0.05])
                for (r, v) in zip(valid, mask):
                    r['quality_Pareto'] = bool(v)
            newidx += 1
        print('scored', key, len(rows), flush=True)
    summary = []
    for family in panel['families']:
        for level in panel['levels']:
            for name in names + list(preds):
                rr = [r for r in rows if r['family'] == family and r['level'] == level and (r['participant'] == name)]
                pp = [r for r in rr if r['L1'] == 'PASS']
                nr = [r for r in pp if r.get('meaningful_conditional_gain')]
                summary.append(dict(family=family, level=level, participant=name, tasks=len(rr), L1_PASS=len(pp), statuses=dict(collections.Counter((r['L1'] for r in rr))), quality_Pareto=sum((r.get('quality_Pareto', False) for r in rr)), conditional_gain_gt_2pct=len(nr), median_force_ratio=float(np.median([r['optimization']['model_force_quantile_ratio'] for r in pp])) if pp and name in preds else None, mean_anatomy_rmse_mm=float(np.mean([r['anatomy_rmse_mm'] for r in pp if r.get('anatomy_rmse_mm') is not None])) if pp else None))
    dump(ROOT / 'raw/SCORED_ROWS.json', rows)
    dump(ROOT / 'raw/R2.json', dict(summary=summary, task_count=len(rows) // 10, participant_count=10, controls=checks, adjoint_checks=fem, all_faults_rejected=all((c['rejected'] for c in checks)), all_adjoint_checks=all((c['pass_gate'] for c in fem)), seconds=time.perf_counter() - start, source_unavailable_tasks=sum((r['L1'] == 'UNKNOWN_SITE' for r in rows)) // 10))
    state('R2_SCORED', 'See raw/R2.json; all failures retained', 'R3 whole-crown geometry and batch-cancelling contrast')
    (ROOT / 'HANDOFF_R2.md').write_text('R2: real isolated callable plugin ran at all frozen levels with every inherited projected obstacle constraint. Prediction hashes precede withheld-reference scoring. All per-task geometry/quality and optimizer failures are in raw/SCORED_ROWS.json. Only conditional rigid-support volume-hazard gain is supported; calibrated force and physical Pareto are UNKNOWN. Next construction changes from affine roof to whole STS crown, with matched batch-cancelling prospective measurement.\n')
    return rows
if __name__ == '__main__':
    run()
