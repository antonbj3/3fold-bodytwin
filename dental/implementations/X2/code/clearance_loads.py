"""Area control and conditional NNLS tangent loads on the fixed 3D proximity graph."""
import json, zipfile, time
import numpy as np
from occlusion_operator import H, DATA, ZIP, read_stl, predict_regions, dump, sha
from feasible_forces import nnls_tangent, bounds_and_basis

def main():
    z = zipfile.ZipFile(ZIP)
    models = json.loads((H / 'raw/segmentation_model.json').read_text())
    k = np.tile([2, 1, 2, 2, 2, 1, 2, 2], 2) * np.sqrt(500 * 1130)
    feasible = []
    cache = {}
    start = time.perf_counter()
    for case in range(1, 201):
        outfile = H / 'raw/cases_3d' / f'{case:03d}.json'
        if outfile.exists():
            feasible.append(json.loads((H / 'raw/feasible_3d' / outfile.name).read_text()))
            continue
        r = json.loads((H / 'raw/cases' / f'{case:03d}.json').read_text())
        c = json.loads((H / 'raw/clearance' / f'{case:03d}.json').read_text())
        (L, _) = read_stl(z, r['upper_member'].replace('upper.stl', 'lower.stl'))
        (regions, _, _) = predict_regions(L.mean(1), models['lower'], True)
        with np.load(c['data_path']) as p:
            d = p['distance_mm']
            ur = p['upper_region']
            face = p['nearest_lower_face']
        metrics = []
        fr = []
        for m in c['metrics']:
            band = m['band_mm']
            ok = d <= band
            if not ok.any():
                metrics.append(dict(band_mm=band, mechanics=None, area_shares_pp=None, n_pixels=0))
                fr.append(dict(band_mm=band, status='NO_CONTACT'))
                continue
            pairs = np.unique(np.c_[ur[ok], regions[face[ok]]], axis=0)
            ms = time.perf_counter()
            n = nnls_tangent(pairs, k)
            elapsed = time.perf_counter() - ms
            n['pairs'] = pairs.tolist()
            key = str(pairs)
            if key not in cache:
                cache[key] = bounds_and_basis(pairs)
            fr.append(dict(band_mm=band, **cache[key]))
            metrics.append(dict(band_mm=band, mechanics=n, area_shares_pp=m['threeD_area_shares_pp'], n_pixels=int(ok.sum()), projected_area_mm2=float(ok.sum() * 0.2 ** 2), mechanics_query_s=elapsed, area_query_s=0, clearance='3D unsigned fixed-pose source-point distance<=band'))
        r.update(metrics=metrics, contact_version='R4_fixed_pose_3D_sampled', clearance_path=c['data_path'], clearance_sha256=c['data_sha256'])
        dump(outfile, r)
        frow = dict(case=case, metrics=fr)
        dump(H / 'raw/feasible_3d' / outfile.name, frow)
        feasible.append(frow)
    central = [m for x in feasible for m in x['metrics'] if m['band_mm'] == 0.1 and 'affine_dimension' in m]
    dump(H / 'round4/force_feasibility_3d.json', dict(n_cases=len(feasible), identified_joint_vectors=sum((m['affine_dimension'] == 0 for m in central)), no_contact_cases=[x['case'] for x in feasible if x['metrics'][1]['status'] == 'NO_CONTACT'], rank_distribution={str(i): sum((m['affine_dimension'] == i for m in central)) for i in sorted(set((m['affine_dimension'] for m in central)))}, minimum_measurement_channels_median=float(np.median([m['affine_dimension'] for m in central])), max_LP_parity_error=max((m['LP_parity_error'] for m in central)), wall_s=time.perf_counter() - start, scope='Necessary outer set on3Dproximity graph; anatomical transfer remains failed', rows=feasible))
if __name__ == '__main__':
    main()
