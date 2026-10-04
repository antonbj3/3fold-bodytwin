from dental_release.paths import expand as _release_expand
from pathlib import Path
import sys, json, gzip, struct, hashlib, time, zipfile, io
import numpy as np
P = Path(__file__).resolve().parent
PRE = json.loads((P / 'PREREG_R1.json').read_text())
SOURCE = Path(PRE['data_root'])
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X24_cbct_hu'))
sys.path.insert(0, str(Path(PRE['source_cell']).parent))
import cbct_multidevice_sensitometry_overdet as old

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda : f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def prepare_rois(device):
    rois = {}
    rows = []
    for (mat, kws) in old.KW.items():
        fs = sorted([p for p in device.glob('**/ROI/*.nii.gz') if any((k in p.name.lower() for k in kws)) and 'air2' not in p.name.lower()])
        if not fs:
            raise ValueError('missing ' + str(device) + ' ' + mat)
        roi = fs[0]
        rois[mat] = roi
        rows.append({'path': str(roi), 'sha256': sha(roi), 'bytes': roi.stat().st_size, 'material': mat})
    return (rois, rows)

def sample(img, A, rois):
    inv = np.linalg.inv(A)
    rows = []
    samples = {}
    for (mat, roi) in rois.items():
        (mask, Ar) = old.read_nii(roi)
        count = np.zeros(img.size, dtype=np.uint32)
        nroi = 0
        maxdisp = 0.0
        for z0 in range(0, mask.shape[0], 4):
            idx = np.argwhere(mask[z0:z0 + 4] > 0.5)
            idx[:, 0] += z0
            nroi += len(idx)
            xyz = np.column_stack([idx[:, 2], idx[:, 1], idx[:, 0], np.ones(len(idx))])
            vi = (inv @ (Ar @ xyz.T)).T[:, :3]
            near = np.round(vi).astype(int)
            (ii, jj, kk) = near.T
            ok = (ii >= 0) & (ii < img.shape[2]) & (jj >= 0) & (jj < img.shape[1]) & (kk >= 0) & (kk < img.shape[0])
            lin = (kk[ok] * img.shape[1] + jj[ok]) * img.shape[2] + ii[ok]
            (ids, weights) = np.unique(lin, return_counts=True)
            count[ids] += weights.astype(np.uint32)
            if ok.any():
                maxdisp = max(maxdisp, float(np.max(np.abs(vi[ok] - near[ok]))))
        ids = np.flatnonzero(count)
        weights = count[ids].astype(np.int64)
        values = img.ravel()[ids].astype(float)
        valid = np.isfinite(values)
        (ids, values, weights) = (ids[valid], values[valid], weights[valid])
        order = np.argsort(values)
        sv = values[order]
        cw = np.cumsum(weights[order])
        n = int(weights.sum())

        def kth(k):
            return float(sv[np.searchsorted(cw, k + 1)])

        def pct(q):
            pos = (n - 1) * q
            lo = int(np.floor(pos))
            return kth(lo) + (kth(int(np.ceil(pos))) - kth(lo)) * (pos - lo)
        mean = float(values @ weights / n)
        sd = float(np.sqrt((values - mean) ** 2 @ weights / n))
        rows.append({'material': mat, 'gray_median': pct(0.5), 'gray_mean': mean, 'gray_SD_voxels': sd, 'gray_p05_p95': [pct(0.05), pct(0.95)], 'n_roi': nroi, 'n_in_grid': n, 'n_finite': n, 'n_unique_image_points': len(ids), 'unique_point_median': float(np.median(values)), 'max_nearest_grid_displacement_vox': maxdisp, 'resolution': 'PER_SURFACE_REGION', 'uncertainty_note': 'legacy nearest mapping repeats image points; SD is not independent SE'})
        for (k, v) in [('image_flat_index', ids), ('gray', values), ('legacy_multiplicity', weights)]:
            samples[mat + '_' + k] = v
        del mask, count
    return (rows, samples)

def baseline():
    t = time.perf_counter()
    manifest = []
    out = []
    for device in sorted((SOURCE / 'extracted').iterdir()):
        (rois, rmanifest) = prepare_rois(device)
        manifest += rmanifest
        fs = sorted(device.glob('**/*test-retest1_*.nii.gz'))
        assert len(fs) == 1, fs
        f = fs[0]
        (img, A) = old.read_nii(f)
        (rows, samples) = sample(img, A, rois)
        arr = DATA / (device.name + '_baseline_ROI.npz')
        np.savez_compressed(arr, **samples)
        out.append({'device': device.name, 'scan': 'baseline_repetition1', 'source': str(f), 'sha256': sha(f), 'bytes': f.stat().st_size, 'image_shape_zyx': list(img.shape), 'image_affine': A.tolist(), 'spacing_mm': np.linalg.norm(A[:3, :3], axis=0).tolist(), 'inserts': rows, 'raw_samples': {'path': str(arr), 'sha256': sha(arr), 'bytes': arr.stat().st_size}})
        print(device.name, {r['material']: r['gray_median'] for r in rows}, flush=True)
        del img, samples
    (P / 'PHANTOM_BASELINE.json').write_text(json.dumps(out, indent=2) + '\n')
    (P / 'SOURCE_MANIFEST_ROIS.json').write_text(json.dumps(manifest, indent=2) + '\n')
    (P / 'MEASUREMENT_COST_R1.json').write_text(json.dumps({'wall_s': time.perf_counter() - t, 'scans_read': len(out), 'new_physical_acquisitions': 0, 'raw_bytes': sum((r['raw_samples']['bytes'] for r in out))}, indent=2) + '\n')
if __name__ == '__main__':
    baseline()
