import sys, time, json, datetime, resource
import numpy as np
from geometry import *
from experiment import local_tri, contact, score, summarize
from obstacle import constraints, hessian, qp, lp_control

def predict():
    pr = json.loads((H / 'PREREG_R2.json').read_text())
    start = time.perf_counter()
    template = np.load(D / 'donor_template.npz')
    rows = json.loads((H / 'raw/PREDICTIONS_R2_checkpoint.json').read_text()) if (H / 'raw/PREDICTIONS_R2_checkpoint.json').exists() else []
    done = {(r['case'], r['fdi']) for r in rows}
    for case in pr['test_cases']:
        state('R2_DESIGNING', 'R1 FAIL preserved; new cases and local continuous constraints', 'Cusp-lift obstacle ' + str(case), completed=len(rows))
        (data, man) = pair(case)
        for fdi in pr['tooth_fdi']:
            if (case, fdi) in done:
                continue
            t = time.perf_counter()
            try:
                s = site(data['lower'], fdi, pr['grid_n'])
                xy = s['xy']
                U = local_tri(data['upper']['tri'], xy)
                (ceiling, uf) = query_height(U, xy, True)
                z0 = s['z_top'] + template[str(fdi % 10)]
                (A, b) = constraints(U, xy, s['faces'], pr['clearance_mm'])
                HH = hessian(s['faces'], len(xy), pr['morphology_regularizer'])
                attempts = []
                selected = None
                first_contact = None
                for lift in pr['lift_candidates_mm']:
                    tc = time.perf_counter()
                    (z, info) = qp(A, b, HH, z0 + lift, pr['dual_maxiter'], pr['dual_gtol'])
                    (z, repairinfo) = repair(U, s, z, pr['clearance_mm'])
                    co = contact(ceiling - z, z, xy, s['index'], s['area_weight'], 0.1, 0.1)
                    sel = pr['selection']
                    eligible = sel['patch_min'] <= co['patch_count'] <= sel['patch_max'] and sel['area_min_mm2'] <= co['area_mm2'] <= sel['area_max_mm2']
                    attempts.append(dict(lift_mm=lift, qp=info, facet_repair=repairinfo, area_mm2=co['area_mm2'], patch_count=co['patch_count'], selected_gate=eligible, seconds=time.perf_counter() - tc))
                    if first_contact is None and co['patch_count'] > 0:
                        first_contact = (z, lift)
                    if selected is None:
                        selected = (z, lift)
                    if eligible:
                        selected = (z, lift)
                        break
                else:
                    if first_contact is not None:
                        selected = first_contact
                (zz, lift) = selected
                tlp = time.perf_counter()
                (ctrl, ci) = lp_control(A, b, z0 + lift)
                if ctrl is None:
                    ctrl = np.full_like(zz, np.nan)
                    controlrepair = None
                else:
                    (ctrl, controlrepair) = repair(U, s, ctrl, pr['clearance_mm'])
                stem = str(case) + '_' + str(fdi) + '_R2'
                path = D / (stem + '_pred.npz')
                np.savez_compressed(path, xy=xy, uv=s['uv'], faces=s['faces'], border=s['border'], index=s['index'], weights=s['area_weight'], ceiling=ceiling, z_practice=z0, z_informed=zz, z_control=ctrl, center=s['center'], half=s['half'])
                for (arm, z) in [('practice', z0), ('informed', zz), ('control', ctrl)]:
                    if np.isfinite(z).all():
                        (v, f) = shell(xy, z, s['faces'], 1.2)
                        write_stl(D / (stem + '_' + arm + '.stl'), v, f)
                rows.append(dict(case=case, fdi=fdi, status='DESIGNED', file=str(path), sha256=sha(path), source=man, lift_selected_mm=lift, attempts=attempts, control=ci, control_facet_repair=controlrepair, control_seconds=time.perf_counter() - tlp, constraints=A.shape[0], design_seconds=time.perf_counter() - t))
                dump(H / 'raw/PREDICTIONS_R2_checkpoint.json', rows)
            except (ValueError, IndexError) as e:
                rows.append(dict(case=case, fdi=fdi, status='UNKNOWN_SITE', reason=str(e)))
        print('R2 designed', case, len(rows), flush=True)
    dump(H / 'raw/PREDICTIONS_R2.json', rows)
    dump(H / 'FROZEN_PREDICTIONS_R2.json', dict(round='R2', frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), prereg_sha256=sha(H / 'PREREG_R2.json'), predictions_sha256=sha(H / 'raw/PREDICTIONS_R2.json'), files=[dict(path=r['file'], sha256=r['sha256']) for r in rows if r['status'] == 'DESIGNED'], code={p.name: sha(p) for p in (H / 'code').glob('*.py')}, target_contacts_read=False))
    dump(H / 'raw/PREDICT_R2_COST.json', dict(wall_seconds=time.perf_counter() - start, peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    state('R2_FROZEN', 'New target contacts unread before prediction freeze', 'Evaluate second-round withheld original roofs')

def evaluate():
    pr = json.loads((H / 'PREREG_R2.json').read_text())
    fr = json.loads((H / 'FROZEN_PREDICTIONS_R2.json').read_text())
    assert sha(H / 'PREREG_R2.json') == fr['prereg_sha256']
    assert sha(H / 'raw/PREDICTIONS_R2.json') == fr['predictions_sha256']
    preds = json.loads((H / 'raw/PREDICTIONS_R2.json').read_text())
    rows = []
    current = None
    start = time.perf_counter()
    for p in preds:
        if p['status'] != 'DESIGNED':
            rows.append(p)
            continue
        if current != p['case']:
            (data, _) = pair(p['case'])
            current = p['case']
            state('R2_EVALUATING', 'Predictions hash matched', 'Second-round original contacts ' + str(current))
        assert sha(p['file']) == p['sha256']
        z = np.load(p['file'])
        xy = z['xy']
        a = data['lower']
        tri = a['tri'][a['owner'] == p['fdi']]
        (original, of) = query_height(tri, xy, False)
        U = local_tri(data['upper']['tri'], xy)
        ceiling = z['ceiling']
        raw = dict(xy=xy, original_z=original, original_faces=of, original_gap=ceiling - original, weights=z['weights'], faces=z['faces'], index=z['index'])
        row = dict(case=p['case'], fdi=p['fdi'], type={4: 'P1', 5: 'P2', 6: 'M1', 7: 'M2'}[p['fdi'] % 10], status='SCORED', arms={}, continuous={}, design_seconds=p['design_seconds'], lift_selected_mm=p['lift_selected_mm'], constraints=p['constraints'])
        for band in [0.1, 0.05, 0.2]:
            truth = contact(ceiling - original, original, xy, z['index'], z['weights'], band, 0.1)
            key = str(band)
            row['arms'][key] = {'original': truth}
            for arm in ['practice', 'informed', 'control']:
                zz = z['z_' + arm]
                co = contact(ceiling - zz, zz, xy, z['index'], z['weights'], band, 0.1)
                if arm == 'control' and (not np.isfinite(zz).all()):
                    co['errors'] = None
                    co['status'] = 'UNKNOWN_LP_TIMEOUT_OR_FAILURE'
                else:
                    co['errors'] = score(co, truth, xy)
                    co['status'] = 'SCORED'
                row['arms'][key][arm] = co
            for (arm, co) in row['arms'][key].items():
                raw[key + '_' + arm + '_contact'] = co.pop('mask')
        for arm in ['practice', 'informed', 'control']:
            row['continuous'][arm] = signed_gap(U, xy, z['z_' + arm], z['faces']) if np.isfinite(z['z_' + arm]).all() else dict(maximum_penetration_mm=None)
        row['pose_sensitivity'] = {str(delta): contact(ceiling + delta - z['z_informed'], z['z_informed'], xy, z['index'], z['weights'], 0.1, 0.1)['area_mm2'] for delta in [-0.05, 0.05]}
        np.savez_compressed(D / (str(p['case']) + '_' + str(p['fdi']) + '_R2_reference.npz'), **raw)
        rows.append(row)
    summary = summarize(rows, pr)
    dump(H / 'raw/RESULTS_R2_ROWS.json', rows)
    dump(H / 'rounds/R2.json', dict(claim_type='information_link', summary=summary, external_referent=pr['external_referent'], cost=dict(wall_seconds=time.perf_counter() - start, peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)))
    state('R2_DECIDED', summary['decision'], 'Independent controls, export and geometry-to-load uncertainty')
    print(json.dumps(summary, indent=2))
if __name__ == '__main__':
    {'predict': predict, 'evaluate': evaluate}[sys.argv[1]]()
