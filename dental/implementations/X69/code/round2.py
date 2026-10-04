import time, resource
import numpy as np
from registration import *

def run_r2():
    start = time.monotonic()
    base = DATA / 'hao_demo/Demo/Demo_1'
    cbct = stl(base / 'CBCT/reconstruction/out_smoothed.stl')
    cent = cbct.mean(1)
    rng = np.random.default_rng(69)
    rows = []
    witnesses = {}
    allresults = []
    for (jaw, folder) in [('lower', 'L'), ('upper', 'U')]:
        tris = cbct[cent[:, 2] < 39] if jaw == 'lower' else cbct[cent[:, 2] >= 39]
        points = []
        fdis = []
        localregions = []
        labels = []
        for file in sorted((base / f'IOS/teeth_seg/{folder}').glob('*.stl')):
            fdi = int(file.stem)
            tri = stl(file)
            p = sample(tri, 600, rng)
            mid = np.median(p[:, 2])
            points.append(p)
            fdis.extend([fdi] * len(p))
            localregions.extend((p[:, 2] >= mid).astype(int).tolist())
            labels.append({'fdi': fdi, 'source': str(file), 'sha256': sha(file), 'triangle_count': len(tri), 'unit': 'mm assumed, STL has no units'})
        p = np.concatenate(points)
        fdi = np.array(fdis)
        region = np.array(localregions)
        fit = ~np.isin(fdi, [11, 21, 31, 41])
        v = tris.reshape(-1, 3)
        band = (v[:, 2] >= 31) & (v[:, 2] < 39) if jaw == 'lower' else (v[:, 2] >= 39) & (v[:, 2] <= 47)
        qpool = v[band]
        q = qpool[rng.choice(len(qpool), min(len(qpool), 18000), replace=False)]
        starts = []
        best = None
        for (i, T) in enumerate(initializations(p[fit], q)):
            (T, loss, n) = icp(p[fit], q, T)
            starts.append({'id': i, 'loss_mm2': loss, 'iterations': n, 'transform': T.tolist()})
            if best is None or loss < best[0]:
                best = (loss, T)
            print('R2', jaw, 'start', i, 'loss', loss, flush=True)
        T = best[1]
        verts = tris.reshape(-1, 3)
        faces = np.arange(len(verts)).reshape(-1, 3)
        (d, w, cell) = nearest_surface(transform(p, T), verts, faces)
        (control, c_loss, c_iters) = icp(p, q, T)
        (dc, wc, cc) = nearest_surface(transform(p, control), verts, faces)
        toothrows = []
        for tooth in sorted(set(fdis)):
            mask = fdi == tooth
            r = {'jaw': jaw, 'fdi': int(tooth), 'heldout_from_fit': tooth in [11, 21, 31, 41], **stats(d[mask]), 'resolution': 'PER_TOOTH', 'gate_pass': bool(np.percentile(d[mask], 95) <= 0.5), 'control_p95_mm': float(np.percentile(dc[mask], 95))}
            toothrows.append(r)
            for regionid in [0, 1]:
                sel = mask & (region == regionid)
                rows.append({'jaw': jaw, 'fdi': int(tooth), 'spatial_z_half': regionid, 'heldout_from_fit': r['heldout_from_fit'], **stats(d[sel]), 'gate_pass': bool(np.percentile(d[sel], 95) <= 0.5)})
        witnesses.update({jaw + '_source_points': p, jaw + '_transformed_points': transform(p, T), jaw + '_closest_cbct_surface': w, jaw + '_surface_cell': cell, jaw + '_residual_mm': d, jaw + '_fdi': fdi, jaw + '_region': region, jaw + '_T': T, jaw + '_control_T': control})
        allresults.append({'jaw': jaw, 'source_to_target': T.tolist(), 'fit_summary': stats(d[fit]), 'heldout_summary': stats(d[~fit]), 'tooth_results': toothrows, 'heldout_gate_pass': all((r['gate_pass'] for r in toothrows if r['heldout_from_fit'])), 'starts': starts, 'source_labels': labels, 'control': {'kind': 'Conventional full-support ICP, same source information', 'fit_loss_mm2': c_loss, 'iterations': c_iters, 'heldout_summary': stats(dc[~fit]), 'source_to_target': control.tolist()}, 'target_triangles': len(tris), 'target_band_point_count': len(qpool), 'crown_fit_samples': len(q), 'rejected_target_band_fraction': 1 - len(qpool) / len(v)})
    output = DATA / 'r2_witnesses.npz'
    np.savez_compressed(output, **witnesses)
    r = {'construction': 'R2', 'claim_type': ['capability', 'information_link'], 'case': 'Hao2023_Demo_1', 'actual_source_pair': True, 'jaws': allresults, 'regions': rows, 'heldout_gate_pass': all((x['heldout_gate_pass'] for x in allresults)), 'independent_registration_TRE': 'UNKNOWN_NO_REFERENCE_TRANSFORM_OR_FIDUCIALS', 'uniform_error_bound_mm': None, 'metric_scale': 'mm assumption from publisher pipeline; no STL embedded units or independent calibration', 'seconds': time.monotonic() - start, 'maxrss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'witness_path': str(output), 'witness_sha256': sha(output), 'external_referent': {'kind': 'published_dataset', 'locator': 'https://doi.org/10.5281/zenodo.8027553', 'compared_quantity': 'separately acquired source-paired IOS crown and CBCT-only tooth reconstruction consistency', 'refutes_us': True}, 'publisher_fused_output': 'REJECTED_AS_INDEPENDENT_ACCURACY_FACIT: IOS crown replacement is circular'}
    dump(ROOT / 'raw/R2_RESULTS.json', r)
    print('R2 summary', [(j['jaw'], j['heldout_summary']['p95_mm']) for j in allresults], flush=True)
    small = cbct[::max(1, len(cbct) // 6000)].mean(1)
    np.savez_compressed(DATA / 'figure_r2.npz', target=small, **witnesses)
if __name__ == '__main__':
    run_r2()
