"""Design first, hash-freeze; evaluate withheld contacts in a separate phase."""
import sys, time, json, datetime, resource
import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree
from geometry import *

def local_tri(tri, xy):
    lo = xy.min(0)
    hi = xy.max(0)
    return tri[np.all(tri[:, :, :2].max(1) >= lo, axis=1) & np.all(tri[:, :, :2].min(1) <= hi, axis=1)]

def predict():
    pr = json.loads((H / 'PREREG_R1.json').read_text())
    start = time.perf_counter()
    D.mkdir(exist_ok=True)
    state('R1_FITTING', 'PREREG frozen; target contacts unread', 'Fit independent donor roof')
    template = donor_template(pr)
    rows = []
    for case in pr['test_cases']:
        state('R1_DESIGNING', 'No target contact scores read', 'Design crowns ' + str(case), completed=len(rows))
        (data, manifest) = pair(case)
        a = data['lower']
        for fdi in pr['tooth_fdi']:
            t = time.perf_counter()
            try:
                s = site(a, fdi, pr['grid_n'])
                U = local_tri(data['upper']['tri'], s['xy'])
                (ceiling, uf) = query_height(U, s['xy'], True)
                z0 = s['z_top'] + template[str(fdi % 10)]
                z1 = np.minimum(z0, ceiling - pr['clearance_mm'])
                (z1, repair_info) = repair(U, s, z1, pr['clearance_mm'])
                zctrl = z0 - np.maximum(z0 - ceiling + pr['clearance_mm'], 0)
                (zctrl, ctrl_repair) = repair(U, s, zctrl, pr['clearance_mm'])
                stem = str(case) + '_' + str(fdi)
                p = D / (stem + '_pred.npz')
                np.savez_compressed(p, xy=s['xy'], uv=s['uv'], faces=s['faces'], border=s['border'], index=s['index'], weights=s['area_weight'], ceiling=ceiling, upper_faces=uf, z_practice=z0, z_informed=z1, z_control=zctrl, center=s['center'], half=s['half'], z_cervical=s['z_cervical'])
                for (arm, z) in [('practice', z0), ('informed', z1)]:
                    (v, f) = shell(s['xy'], z, s['faces'], pr['fe_contract']['thickness_mm'])
                    write_stl(D / (stem + '_' + arm + '.stl'), v, f)
                rows.append(dict(case=case, fdi=fdi, status='DESIGNED', file=str(p), sha256=sha(p), source=manifest, repair=repair_info, control_repair=ctrl_repair, control_max_height_difference_mm=float(np.max(np.abs(z1 - zctrl))), site=dict(center_mm=s['center'], halfspan_mm=s['half'], z_top_from_neighbors_mm=s['z_top'], z_cervical_mm=s['z_cervical']), design_seconds=time.perf_counter() - t))
            except (ValueError, IndexError) as e:
                rows.append(dict(case=case, fdi=fdi, status='UNKNOWN_SITE', reason=str(e)))
        print('Designed', case, 'rows', len(rows), flush=True)
    dump(H / 'raw/PREDICTIONS_R1.json', rows)
    dump(H / 'FROZEN_PREDICTIONS.json', dict(round='R1', frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), prereg_sha256=sha(H / 'PREREG_R1.json'), code={p.name: sha(p) for p in (H / 'code').glob('*.py')}, predictions_sha256=sha(H / 'raw/PREDICTIONS_R1.json'), files=[dict(path=r['file'], sha256=r['sha256']) for r in rows if r['status'] == 'DESIGNED'], target_contacts_read=False, measurement_status='Independent scan-derived original reference evaluated only after this freeze; future physical measurement NOT_RUN'))
    dump(H / 'raw/PREDICT_COST.json', dict(wall_seconds=time.perf_counter() - start, peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    state('R1_PREDICTIONS_FROZEN', 'All design files hashed before original contact scoring', 'Evaluate original measured surfaces')

def contact(gap, z, xy, index, weights, band, min_area):
    mask = np.isfinite(gap) & (np.abs(gap) <= band)
    im = np.zeros(index.shape, bool)
    good = index >= 0
    im[good] = mask[index[good]]
    (labels, n) = ndimage.label(im)
    patch = []
    filtered = np.zeros(len(mask), bool)
    for k in range(1, n + 1):
        ids = index[labels == k]
        area = float(weights[ids].sum())
        if area < min_area:
            continue
        filtered[ids] = True
        patch.append(dict(area_mm2=area, centroid_xyz_mm=np.average(np.c_[xy[ids], z[ids]], axis=0, weights=weights[ids]), point_indices=ids.tolist()))
    return dict(mask=filtered, area_mm2=float(weights[filtered].sum()), patch_count=len(patch), patches=patch, minimum_vertex_gap_mm=float(np.min(gap[np.isfinite(gap)])) if np.isfinite(gap).any() else None, penetrating_vertex_area_mm2=float(weights[np.isfinite(gap) & (gap < 0)].sum()))

def score(pred, truth, xy):
    (p, t) = (pred['mask'], truth['mask'])
    union = np.sum(p | t)
    iou = float(np.sum(p & t) / union) if union else 1.0
    pc = np.array([x['centroid_xyz_mm'][:2] for x in pred['patches']])
    tc = np.array([x['centroid_xyz_mm'][:2] for x in truth['patches']])
    if len(pc) and len(tc):
        loc = float((cKDTree(tc).query(pc)[0].mean() + cKDTree(pc).query(tc)[0].mean()) / 2)
    elif not len(pc) and (not len(tc)):
        loc = 0.0
    else:
        loc = float(np.linalg.norm(xy.max(0) - xy.min(0)))
    return dict(IoU=iou, absolute_area_error_mm2=abs(pred['area_mm2'] - truth['area_mm2']), absolute_patch_count_error=abs(pred['patch_count'] - truth['patch_count']), centroid_error_mm=loc)

def evaluate():
    pr = json.loads((H / 'PREREG_R1.json').read_text())
    fr = json.loads((H / 'FROZEN_PREDICTIONS.json').read_text())
    assert sha(H / 'PREREG_R1.json') == fr['prereg_sha256']
    assert sha(H / 'raw/PREDICTIONS_R1.json') == fr['predictions_sha256']
    preds = json.loads((H / 'raw/PREDICTIONS_R1.json').read_text())
    rows = []
    start = time.perf_counter()
    data = None
    current = None
    for r in preds:
        if r['status'] != 'DESIGNED':
            rows.append(r)
            continue
        if current != r['case']:
            (data, _) = pair(r['case'])
            current = r['case']
            state('R1_EVALUATING', 'Predictions frozen and hash matched', 'Read withheld original surface ' + str(current), completed=len(rows))
        assert sha(r['file']) == r['sha256']
        z = np.load(r['file'])
        xy = z['xy']
        f = z['faces']
        weights = z['weights']
        index = z['index']
        U = local_tri(data['upper']['tri'], xy)
        a = data['lower']
        C = a['tri'][a['owner'] == r['fdi']]
        (original, of) = query_height(C, xy, False)
        ceiling = z['ceiling']
        gap = ceiling - original
        raw = dict(xy=xy, original_z=original, original_faces=of, original_gap=gap, ceiling=ceiling, faces=f, weights=weights, index=index)
        out = dict(case=r['case'], fdi=r['fdi'], type={4: 'P1', 5: 'P2', 6: 'M1', 7: 'M2'}[r['fdi'] % 10], status='SCORED', design_seconds=r['design_seconds'], control_max_height_difference_mm=r['control_max_height_difference_mm'], arms={})
        for band in [pr['contact_band_mm']] + pr['sensitivity_bands_mm']:
            truth = contact(gap, original, xy, index, weights, band, pr['minimum_patch_area_mm2'])
            key = str(band)
            out['arms'][key] = {'original': truth}
            for arm in ['practice', 'informed']:
                zz = z['z_' + arm]
                co = contact(ceiling - zz, zz, xy, index, weights, band, pr['minimum_patch_area_mm2'])
                co['errors'] = score(co, truth, xy)
                out['arms'][key][arm] = co
            for (arm, co) in out['arms'][key].items():
                raw[key + '_' + arm + '_contact'] = co.pop('mask')
        out['continuous'] = {arm: signed_gap(U, xy, z['z_' + arm], f) for arm in ['practice', 'informed']}
        out['original_scan'] = {'finite_roof_fraction': float(np.isfinite(original).mean()), 'negative_gap_fraction': float(np.mean(gap[np.isfinite(gap)] < 0)) if np.isfinite(gap).any() else None}
        s = dict(xy=xy, faces=f)
        out['pose_sensitivity'] = {str(delta): contact(ceiling + delta - z['z_informed'], z['z_informed'], xy, index, weights, 0.1, 0.1)['area_mm2'] for delta in pr['metrics']['rigid_pose_shift_mm']}
        np.savez_compressed(D / (str(r['case']) + '_' + str(r['fdi']) + '_reference.npz'), **raw)
        rows.append(out)
    dump(H / 'raw/RESULTS_R1_ROWS.json', rows)
    summary = summarize(rows, pr)
    dump(H / 'rounds/R1.json', dict(claim_type='information_link', summary=summary, rows_file='raw/RESULTS_R1_ROWS.json', external_referent=pr['external_referent'], cost=dict(wall_seconds=time.perf_counter() - start, peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)))
    state('R1_DECIDED', summary['decision'], 'Write handoff; change clipping-only construction to sparse contact-constrained cusp editing')
    print(json.dumps(summary, indent=2))

def summarize(rows, pr):
    rows = [r for r in rows if r['status'] == 'SCORED']
    key = str(pr['contact_band_mm'])
    out = {}
    for kind in ['all', 'P1', 'P2', 'M1', 'M2']:
        rr = [r for r in rows if kind == 'all' or r['type'] == kind]
        group = {}
        for arm in ['practice', 'informed']:
            a = [r['arms'][key][arm] for r in rr]
            group[arm] = {q: float(np.median([x['errors'][q] for x in a])) for q in ['IoU', 'absolute_area_error_mm2', 'absolute_patch_count_error', 'centroid_error_mm']}
            group[arm]['median_patch_count'] = float(np.median([x['patch_count'] for x in a]))
            group[arm]['max_penetration_mm'] = float(max((r['continuous'][arm]['maximum_penetration_mm'] for r in rr)))
        group['n'] = len(rr)
        group['median_IoU_gain'] = float(np.median([r['arms'][key]['informed']['errors']['IoU'] - r['arms'][key]['practice']['errors']['IoU'] for r in rr]))
        group['area_error_ratio'] = group['informed']['absolute_area_error_mm2'] / max(group['practice']['absolute_area_error_mm2'], 0.1)
        out[kind] = group
    a = out['all']
    passes = dict(nonpenetration=a['informed']['max_penetration_mm'] <= pr['metrics']['penetration_tolerance_mm'], IoU=a['median_IoU_gain'] >= pr['metrics']['contact_IoU_gain_min'], area=a['area_error_ratio'] <= pr['metrics']['contact_area_abs_error_ratio_max'])
    return dict(by_type=out, gates=passes, decision='PASS' if all(passes.values()) else 'FAIL', scope='Scan-derived projected contacts conditional on predicted FDI; no physical contact-force accuracy')
if __name__ == '__main__':
    {'predict': predict, 'evaluate': evaluate}[sys.argv[1]]()
