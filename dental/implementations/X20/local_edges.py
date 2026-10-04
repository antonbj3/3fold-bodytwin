"""R4: stream raw CBCT observations at the decision's active bone constraint."""
import hashlib
import json
from pathlib import Path
import time
import zipfile
import numpy as np
from scipy import ndimage as ndi
from measure import HERE, ROOT, DATA, tf2_io, sha, dump, check_frozen, unit
from capacity_map import load_local, body_points

def image_crop(case, lo, hi):
    member = f'{tf2_io.ROOT}/imagesTr/{case}_0000.mha'
    with zipfile.ZipFile(tf2_io.ZIP) as z, z.open(member) as f:
        header = {}
        header_bytes = b''
        while True:
            line = f.readline()
            header_bytes += line
            if not line:
                raise ValueError('missing MHA local header')
            s = line.decode().strip()
            if ' = ' in s:
                (k, v) = s.split(' = ', 1)
                header[k] = v
                if k == 'ElementDataFile':
                    break
        assert header['ElementDataFile'] == 'LOCAL' and header.get('CompressedData', 'False') == 'False'
        shape = np.asarray([int(v) for v in header['DimSize'].split()][::-1])
        spacing = np.asarray([float(v) for v in header['ElementSpacing'].split()][::-1])
        dtype = np.dtype(tf2_io._DT[header['ElementType']]).newbyteorder('>' if header.get('BinaryDataByteOrderMSB') == 'True' else '<')
        assert np.all(hi <= shape) and np.all(lo >= 0)
        out = np.empty(tuple(hi - lo), dtype=np.float32)
        size = int(shape[1] * shape[2]) * dtype.itemsize
        for k in range(int(hi[0])):
            raw = f.read(size)
            assert len(raw) == size, 'truncated image member'
            if k >= lo[0]:
                plane = np.frombuffer(raw, dtype=dtype).reshape(tuple(shape[1:]))
                out[k - lo[0]] = plane[lo[1]:hi[1], lo[2]:hi[2]]
        info = {'member': member, 'member_crc32': z.getinfo(member).CRC, 'member_bytes': z.getinfo(member).file_size, 'header_sha256': hashlib.sha256(header_bytes).hexdigest(), 'shape': shape.tolist(), 'spacing': spacing.tolist(), 'crop_lo': lo.tolist(), 'crop_hi': hi.tolist(), 'image_dtype': str(dtype)}
    return (out, spacing, info)

def crossings(x, y, level, positive=False):
    ids = np.where((y[:-1] - level) * (y[1:] - level) <= 0)[0]
    out = []
    for i in ids:
        den = y[i + 1] - y[i]
        if abs(den) > 1e-12 and (not positive or den > 0):
            out.append(float(x[i] + (level - y[i]) * (x[i + 1] - x[i]) / den))
    return out

def run():
    start = time.monotonic()
    check_frozen('PREREG_R4')
    check_frozen('FROZEN_PREDICTIONS_R4')
    fr = json.loads((HERE / 'FROZEN_PREDICTIONS_R4.json').read_text())
    assert sha(HERE / 'RAW_R3_CAPACITY_MAPS.jsonl') == fr['source_map_sha256']
    rows = [json.loads(l) for l in (HERE / 'COMBINED_SITES.jsonl').read_text().splitlines()]
    maps = [json.loads(l) for l in (HERE / 'RAW_R3_CAPACITY_MAPS.jsonl').read_text().splitlines()]
    measurements = []
    manifest = []
    for row in [r for r in rows if r['valid']]:
        (local, anchor, axis, normal, tangent) = load_local(row)
        p = np.asarray(row['crest_mm'])
        (points, depths) = body_points(p, axis, normal, tangent, 4.0, 6.0)
        selected = points[depths >= 2 - 1e-09]
        values = local.sample('bone_sdf', selected, order=0, fill=-999)
        witness = selected[int(np.argmin(values))]
        gradient = np.array([(local.sample('bone_sdf', (witness + 0.15 * np.eye(3)[k])[None])[0] - local.sample('bone_sdf', (witness - 0.15 * np.eye(3)[k])[None])[0]) / 0.3 for k in range(3)])
        if np.linalg.norm(gradient) < 1e-06:
            measurements.append({'case': row['case'], 'fdi': row['fdi'], 'valid': False, 'reason': 'zero_sdf_normal'})
            continue
        n = unit(gradient)
        guess = witness - float(local.sample('bone_sdf', witness[None])[0]) * n
        t = np.arange(-2, 2.00001, 0.05)
        pp = guess + t[:, None] * n
        occ = local.sample('bone', pp)
        label_crossings = crossings(t, occ, 0.5, positive=True)
        if not label_crossings:
            measurements.append({'case': row['case'], 'fdi': row['fdi'], 'valid': False, 'reason': 'no_label_surface_crossing'})
            continue
        label_t = min(label_crossings, key=abs)
        surface = guess + label_t * n
        (image, sp, meta) = image_crop(row['case'], local.lo, local.hi)
        assert np.allclose(sp, local.sp)
        coords = ((surface + t[:, None] * n) / sp - local.lo).T
        intensity = ndi.map_coordinates(image, coords, order=1, mode='constant', cval=float('nan'))
        outside = float(np.median(intensity[(t >= -1.5) & (t <= -0.9)]))
        inside = float(np.median(intensity[(t >= 0.9) & (t <= 1.5)]))
        contrast = inside - outside
        half = (outside + inside) / 2
        edge_crossings = [a for a in crossings(t, intensity, half, positive=True) if abs(a) <= 1.5]
        rec = {'case': row['case'], 'fdi': row['fdi'], 'witness_mm': witness.tolist(), 'label_surface_mm': surface.tolist(), 'normal_into_bone': n.tolist(), 'profile_t_mm': t.tolist(), 'intensity': intensity.tolist(), 'outside_gray': outside, 'inside_gray': inside, 'contrast_gray': contrast, 'valid': bool(np.isfinite(intensity).all() and contrast >= 100 and edge_crossings), 'physical_truth': 'UNKNOWN', 'image_meta': meta}
        if rec['valid']:
            offset = min(edge_crossings, key=abs)
            smooth = ndi.gaussian_filter1d(intensity, 1)
            grad = np.gradient(smooth, t)
            eligible = np.where(np.abs(t) <= 1.5)[0]
            gpos = float(t[eligible[np.argmax(grad[eligible])]])
            rec.update(image_minus_label_offset_mm=offset, proxy_box_pass=abs(offset) <= 0.3, gradient_max_offset_mm=gpos, estimator_difference_mm=abs(gpos - offset), estimator_control_pass=abs(gpos - offset) <= 0.3, reason='ok')
        else:
            rec['reason'] = 'low_contrast_or_no_correct_crossing'
        path = DATA / f"{row['case']}_{row['fdi']}_image_crop.npz"
        np.savez_compressed(path, image=image, lo=local.lo, spacing=sp)
        manifest.append({'path': str(path), 'bytes': path.stat().st_size, 'sha256': sha(path), 'source': meta})
        measurements.append(rec)
        print(row['case'], row['fdi'], rec['reason'], rec.get('image_minus_label_offset_mm'), flush=True)
    with open(HERE / 'RAW_R4_LOCAL_EDGES.jsonl', 'w') as f:
        for r in measurements:
            f.write(json.dumps(r, allow_nan=False) + '\n')
    valid = [r for r in measurements if r['valid']]
    fraction = sum((r['proxy_box_pass'] for r in valid)) / len(valid) if valid else None
    summary = {'n_selected': len(measurements), 'n_qc_valid': len(valid), 'fraction_inside_proxy_box': fraction, 'proxy_consistency_gate': 'UNKNOWN_MINIMUM_N' if len(valid) < 8 else 'PASS' if fraction >= 0.9 else 'FAIL', 'estimator_control_failed': sum((not r['estimator_control_pass'] for r in valid)), 'median_image_minus_label_offset_mm': float(np.median([r['image_minus_label_offset_mm'] for r in valid])) if valid else None, 'physical_calibration_gate': 'UNKNOWN_NO_PAIRED_PHYSICAL_TRUTH', 'elapsed_s': time.monotonic() - start, 'external_referent': json.loads((HERE / 'PREREG_R4.json').read_text())['external_referent']}
    dump(HERE / 'DATA_MANIFEST_R4.json', manifest)
    dump(HERE / 'SUMMARY_R4.json', summary)
    dump(HERE / 'CURRENT_WORK_STATE.json', {'lane': 'X20-short-implant-sinus', 'milestone': 'R4_completed', 'latest_gate': summary, 'next_operation': 'Corruption controls, figure, one-command reproducibility, local paired physical measurement port'})
    print(json.dumps(summary, indent=2))
if __name__ == '__main__':
    run()
