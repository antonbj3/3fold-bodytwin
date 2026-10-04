import json, time, zipfile
import numpy as np
from occlusion_operator import H, DATA, ZIP, read_stl, make_map, analyze_map, dump, sha, state
from arc_regions import predict

def main():
    models = json.loads((H / 'raw/segmentation_r2_model.json').read_text())
    z = zipfile.ZipFile(ZIP)
    rows = []
    refine = []
    selected = list(range(1, 13))
    for p in sorted((H / 'raw/cases').glob('*.json')):
        r = json.loads(p.read_text())
        case = r['case']
        out = H / 'raw/cases_r2' / p.name
        if out.exists():
            rows.append(json.loads(out.read_text()))
            continue
        st = time.perf_counter()
        (U, uh) = read_stl(z, r['upper_member'])
        (L, lh) = read_stl(z, r['upper_member'].replace('upper.stl', 'lower.stl'))
        (ur, uf) = predict(U.mean(1), models['upper'], True)
        (lr, lf) = predict(L.mean(1), models['lower'], True)
        with np.load(r['map_path']) as arr:
            m = {k: arr[k] for k in ['gap', 'xy', 'upper_z', 'lower_z', 'upper_face', 'lower_face', 'origin', 'shape']}
        m['grid_mm'] = 0.2
        r2 = dict(r, region_version='R2_arc', metrics=analyze_map(m, ur, lr), upper_frame=uf, lower_frame=lf, relabel_s=time.perf_counter() - st)
        np.savez_compressed(DATA / f'{case:03d}_regions_r2.npz', upper_region=ur[m['upper_face']], lower_region=lr[m['lower_face']])
        r2['regions_path'] = str(DATA / f'{case:03d}_regions_r2.npz')
        r2['regions_sha256'] = sha(r2['regions_path'])
        if case in selected:
            st = time.perf_counter()
            fm = make_map(U, L, 0.1)
            fm['grid_mm'] = 0.1
            fine = analyze_map(fm, ur, lr)
            comparison = []
            for (coarse, f) in zip(r2['metrics'], fine):
                if coarse['mechanics'] is None or f['mechanics'] is None:
                    comparison.append(dict(band_mm=coarse['band_mm'], status='NO_CONTACT', area_L1_pp=None, mechanics_L1_pp=None, gate=False))
                    continue
                ae = float(np.abs(np.asarray(coarse['area_shares_pp']) - f['area_shares_pp']).sum() / 2)
                me = float(np.abs(np.asarray(coarse['mechanics']['shares_pp']) - f['mechanics']['shares_pp']).sum() / 2)
                comparison.append(dict(band_mm=coarse['band_mm'], area_L1_pp=ae, mechanics_L1_pp=me, gate=bool(ae <= 5 and me <= 5)))
            r2['refinement'] = dict(fine_metrics=fine, comparison=comparison, wall_s=time.perf_counter() - st, coarse_mm=0.2, fine_mm=0.1)
            np.savez_compressed(DATA / f'{case:03d}_fine_map.npz', **{k: v for (k, v) in fm.items() if k != 'grid_mm'}, upper_region=ur[fm['upper_face']], lower_region=lr[fm['lower_face']])
            r2['fine_map_path'] = str(DATA / f'{case:03d}_fine_map.npz')
            r2['fine_map_sha256'] = sha(r2['fine_map_path'])
        dump(out, r2)
        rows.append(r2)
        if case % 25 == 0:
            print('relabel/refine', case, flush=True)
    dump(H / 'raw/refinement.json', dict(cases=[dict(case=x['case'], **x['refinement']) for x in rows if 'refinement' in x], all_fine_caseband_L1_gate=all((y['gate'] for x in rows if 'refinement' in x for y in x['refinement']['comparison']))))
    state('R2_RELABEL_FINISHED', 'Fine-grid gates saved', 'Score externally with unchanged R1 thresholds')
if __name__ == '__main__':
    main()
