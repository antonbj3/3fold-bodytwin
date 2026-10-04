from dental_release.paths import expand as _release_expand
import argparse
import csv
import hashlib
import itertools
import json
import os
import resource
import sys
import time
import zipfile
from pathlib import Path
import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from archive_io import index, nifti, read_member, ARCHIVE
ROOT = Path(os.environ.get('X12_RUN_ROOT', str(Path(__file__).resolve().parents[1])))
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X12'))
from tf2_io import ZIP, ROOT as TFROOT, read_mha_bytes
FDI = {v: v for v in list(range(31, 39)) + list(range(41, 49))}

def transformations(source_shape, target_shape):
    for perm in itertools.permutations(range(3)):
        if tuple((source_shape[i] for i in perm)) == tuple(target_shape):
            for flip in itertools.product((False, True), repeat=3):
                yield (perm, flip)

def transform(a, perm, flip):
    return a.transpose(perm)[tuple((slice(None, None, -1) if b else slice(None) for b in flip))]

def transform_coords(coords, source_shape, perm, flip):
    c = coords[:, perm].copy()
    for (axis, b) in enumerate(flip):
        if b:
            c[:, axis] = source_shape[perm[axis]] - 1 - c[:, axis]
    return c

def exact_match(source, target):
    samples = tuple((np.linspace(0, n - 1, min(n, 13), dtype=int) for n in target.shape))
    ix = np.ix_(*samples)
    found = []
    for (p, f) in transformations(source.shape, target.shape):
        candidate = transform(source, p, f)
        if np.array_equal(candidate[ix], target[ix]) and np.array_equal(candidate, target):
            found.append((p, f))
    return found

def align_pulp(pulp, tooth):
    coords = np.argwhere((pulp >= 31) & (pulp <= 38) | (pulp >= 41) & (pulp <= 48))
    labels = pulp[tuple(coords.T)]
    scores = []
    for (p, f) in transformations(pulp.shape, tooth.shape):
        c = transform_coords(coords, pulp.shape, p, f)
        vals = tooth[tuple(c.T)]
        score = float(np.mean(vals == labels)) if len(coords) else 0.0
        scores.append((score, p, f))
    scores.sort(reverse=True)
    if not scores or scores[0][0] < 0.95:
        return (None, scores)
    top = [s for s in scores if s[0] >= 0.95]
    mapped = transform(pulp, top[0][1], top[0][2])
    if any((not np.array_equal(mapped, transform(pulp, s[1], s[2])) for s in top[1:])):
        return (None, scores)
    return (mapped, scores)

def tooth_measure(case, fdi, tooth, pulp, spacing, header, validate=False, box=None):
    if box is not None:
        lo = np.array([max(0, s.start - 2) for s in box])
        hi = np.array([min(n, s.stop + 2) for (s, n) in zip(box, tooth.shape)])
        sl = tuple((slice(int(a), int(b)) for (a, b) in zip(lo, hi)))
        t = tooth[sl] == fdi
        p = pulp[sl] == fdi
        coords = np.argwhere(t) + lo
    else:
        whole = tooth == fdi
        coords = np.argwhere(whole)
        lo = np.maximum(coords.min(0) - 2, 0)
        hi = np.minimum(coords.max(0) + 3, tooth.shape)
        sl = tuple((slice(int(a), int(b)) for (a, b) in zip(lo, hi)))
        t = whole[sl]
        p = pulp[sl] == fdi
    if len(coords) < 200:
        return (None, {'reason': 'TOOTH_TOO_SMALL', 'case': case, 'fdi': fdi})
    npulp = int(np.sum(pulp == fdi))
    fraction = float(np.sum(p & t) / npulp) if npulp else 0
    if npulp < 30 or fraction < 0.95:
        return (None, {'reason': 'PULP_COUNT_OR_CONTAINMENT', 'case': case, 'fdi': fdi, 'pulp_voxels': npulp, 'containment': fraction})
    pc = np.argwhere(p)
    profile = t.sum(axis=(1, 2))
    occupied = np.flatnonzero(profile)
    k = max(1, int(np.ceil(0.4 * len(occupied))))
    lowarea = float(np.mean(profile[occupied[:k]]))
    higharea = float(np.mean(profile[occupied[-k:]]))
    area_contrast = abs(lowarea - higharea) / max(lowarea, higharea, 1)
    if area_contrast < 0.05:
        return (None, {'reason': 'UNKNOWN_CROWN_DIRECTION', 'case': case, 'fdi': fdi, 'area_contrast': area_contrast})
    crown_sign = -1 if lowarea > higharea else 1
    if crown_sign < 0:
        t = t[::-1].copy()
        p = p[::-1].copy()
        pc = np.argwhere(p)
    trunc = bool(np.any(coords.min(0) == 0) or np.any(coords.max(0) == np.array(tooth.shape) - 1))
    boundary = t & ~ndi.binary_erosion(t, structure=ndi.generate_binary_structure(3, 1), border_value=0)
    pb = p & ~ndi.binary_erosion(p, structure=ndi.generate_binary_structure(3, 1), border_value=0)
    dist = ndi.distance_transform_edt(~boundary, sampling=spacing)
    ztop = int(pc[:, 0].max())
    horns = pc[pc[:, 0] == ztop]
    occl = []
    ax0 = []
    ax15 = []
    for (z, y, x) in horns:
        col = np.flatnonzero(t[:, y, x])
        above = col[col >= z]
        if len(above):
            occl.append((int(above.max()) - int(z)) * spacing[0])
        for (zz, out) in [(z, ax0), (max(0, int(z) - round(1.5 / spacing[0])), ax15)]:
            if not t[zz, y, x]:
                continue
            for axis in [1, 2]:
                line = t[zz, :, x] if axis == 1 else t[zz, y, :]
                pos = int(y if axis == 1 else x)
                for sgn in [-1, 1]:
                    j = pos
                    while 0 <= j + sgn < len(line) and line[j + sgn]:
                        j += sgn
                    out.append(abs(j - pos) * spacing[axis])
    cp = pb.copy()
    cp[:max(0, ztop - round(2 / spacing[0]))] = False
    d = dist[cp]
    if not len(d):
        return (None, {'reason': 'EMPTY_HORN_REGION', 'case': case, 'fdi': fdi})
    physical = (coords - coords.mean(0)) * np.array(spacing)
    (eig, vec) = np.linalg.eigh(physical.T @ physical / max(1, len(coords)))
    tilt = float(np.degrees(np.arccos(np.clip(abs(vec[0, -1]), 0, 1))))
    check = None
    if validate:
        tree = cKDTree(np.argwhere(boundary) * np.array(spacing))
        vals = tree.query(np.argwhere(cp) * np.array(spacing), workers=1)[0]
        check = {'case': case, 'fdi': fdi, 'max_abs_error_mm': float(np.max(np.abs(vals - d))), 'n_query': len(vals), 'mutated_distance_rejected': bool(np.max(np.abs(vals - (d + 1))) > 1e-06)}
    r = dict(case=case, fdi=fdi, source_pulpy_fdi=fdi + 10 if fdi < 40 else fdi - 10, semantic_map='left_right_swap', crown_axis_sign=crown_sign, crown_area_contrast=area_contrast, tooth_type={1: 'central_incisor', 2: 'lateral_incisor', 3: 'canine', 4: 'premolar1', 5: 'premolar2', 6: 'molar1', 7: 'molar2', 8: 'molar3'}[fdi % 10], pulp_voxels=npulp, tooth_voxels=len(coords), containment=fraction, domain_truncated=trunc, voxel_mm=float(spacing[0]), pca_tilt_deg=tilt, horn_boundary_min_mm=float(d.min()), horn_boundary_p05_mm=float(np.quantile(d, 0.05)), horn_boundary_median_mm=float(np.median(d)), occlusal_vertical_min_mm=float(min(occl)) if occl else None, axial_cardinal_horn_min_mm=float(min(ax0)) if ax0 else None, axial_cardinal_1p5_below_min_mm=float(min(ax15)) if ax15 else None, digital_two_surface_radius_mm=float(np.linalg.norm(spacing)), horn_boundary_lower_mm=float(max(0, d.min() - np.linalg.norm(spacing))))
    return ((r, check, {'tooth': t, 'pulp': p, 'spacing': np.array(spacing), 'crop_origin_zyx': lo}), None)

def write_csv(path, rows):
    if not rows:
        return
    with path.open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--output', default='R6')
    args = ap.parse_args()
    started = time.time()
    cpu = time.process_time()
    rows = index()
    lookup = {r['name']: r for r in rows}
    source_stat = {str(p): {'size': p.stat().st_size, 'mtime_ns': p.stat().st_mtime_ns} for p in [ARCHIVE, Path(ZIP)]}
    (ROOT / 'raw' / f'{args.output}_source_stat.json').write_text(json.dumps(source_stat, indent=2) + '\n')
    ct = [r for r in rows if r['name'].endswith('/data.nii.gz')]
    ct.sort(key=lambda r: int(r['name'].split('/')[1][1:]))
    pairs = []
    measurements = []
    excluded = []
    controls = []
    manifest = []
    saved = 0
    with zipfile.ZipFile(ZIP) as z:
        names = set(z.namelist())
        for r in ct:
            pid = r['name'].split('/')[1]
            case = 'ToothFairy2P_' + pid[1:].zfill(3)
            iname = f'{TFROOT}/imagesTr/{case}_0000.mha'
            lname = f'{TFROOT}/labelsTr/{case}.mha'
            pname = f'Pulpy3D/{pid}/gt_instance.nii.gz'
            if not r['complete'] or iname not in names or lname not in names or (pname not in lookup):
                pairs.append({'case': pid, 'paired': False, 'reason': 'TRUNCATED_OR_MISSING_CANDIDATE'})
                continue
            t0 = time.time()
            (image, sha) = nifti(r)
            a = np.asanyarray(image.dataobj)
            raw = z.read(iname)
            (b, sp, h) = read_mha_bytes(raw)
            matches = exact_match(a, b)
            record = {'case': pid, 'tf2_case': case, 'paired': bool(matches), 'image_transforms': matches, 'physical_spacing_zyx': sp, 'pulpy_nifti_spacing': [float(x) for x in image.header.get_zooms()], 'pulpy_CT_sha256': sha, 'tf2_CT_sha256': hashlib.sha256(raw).hexdigest(), 'image_voxels': int(a.size), 'identity_wall_s': time.time() - t0}
            if not matches:
                record['reason'] = 'CT_NOT_IDENTICAL'
                pairs.append(record)
                continue
            if not controls:
                corrupted = b.copy()
                corrupted.flat[len(corrupted) // 2] += 1
                record['mutated_CT_rejected'] = not exact_match(a, corrupted)
                del corrupted
            del image, a, b, raw
            labelraw = z.read(lname)
            (tooth, sp_label, lh) = read_mha_bytes(labelraw)
            (pimg, psha) = nifti(lookup[pname])
            pulp = np.asanyarray(pimg.dataobj)
            mapped = None
            if len(matches) == 1:
                pose = transform(pulp, matches[0][0], matches[0][1])
                lut = np.arange(65536, dtype=np.uint16)
                for fdi in FDI:
                    lut[fdi] = fdi + 10 if fdi < 40 else fdi - 10
                mapped = lut[pose]
                coords = np.argwhere(pose > 0)
                v = pose[tuple(coords.T)]
                tv = tooth[tuple(coords.T)]
                swapped = mapped[tuple(coords.T)]
                record['identity_FDI_containment'] = float(np.mean(tv == v)) if len(v) else 0
                record['swapped_FDI_containment'] = float(np.mean(tv == swapped)) if len(v) else 0
                record['aligned'] = bool(record['swapped_FDI_containment'] >= 0.95)
                record['pulp_physical_pose'] = matches[0]
                record['semantic_map'] = 'left_right_swap'
                if not record['aligned']:
                    mapped = None
            else:
                record['aligned'] = False
                record['reason'] = 'NONUNIQUE_IMAGE_POSE'
            record['label_grid_consistent'] = bool(tuple(sp) == tuple(sp_label) and h.get('TransformMatrix') == lh.get('TransformMatrix') and (h.get('Offset') == lh.get('Offset')))
            record['pulpy_pulp_sha256'] = psha
            record['TF2_tooth_sha256'] = hashlib.sha256(labelraw).hexdigest()
            if mapped is not None and record['label_grid_consistent']:
                boxes = ndi.find_objects(tooth, max_label=48)
                for (label, fdi) in FDI.items():
                    if boxes[fdi - 1] is None:
                        excluded.append({'reason': 'MISSING_TOOTH', 'case': pid, 'fdi': fdi})
                        continue
                    (result, exc) = tooth_measure(pid, fdi, tooth, mapped, sp, lh, validate=len(controls) < 10, box=boxes[fdi - 1])
                    if exc:
                        excluded.append(exc)
                        continue
                    (m, control, crop) = result
                    measurements.append(m)
                    if control:
                        controls.append(control)
                    if saved < 3 and fdi % 10 == 6 and (not m['domain_truncated']):
                        file = DATA / f'{args.output}_{pid}_{fdi}_paired.npz'
                        np.savez_compressed(file, **crop)
                        manifest.append({'path': str(file), 'sha256': hashlib.sha256(file.read_bytes()).hexdigest(), 'bytes': file.stat().st_size, 'case': pid, 'fdi': fdi})
                        saved += 1
            pairs.append(record)
            del tooth, pulp, pimg, mapped, labelraw
            if len([p for p in pairs if p.get('paired')]) % 10 == 0 or args.limit:
                print(pid, 'paired', sum((p.get('paired', False) for p in pairs)), 'aligned', sum((p.get('aligned', False) for p in pairs)), 'teeth', len(measurements), flush=True)
            for (name, obj) in [('pairs', pairs), ('excluded', excluded), ('controls', controls), ('manifest', manifest)]:
                (ROOT / 'raw' / f'{args.output}_{name}.json').write_text(json.dumps(obj, indent=1) + '\n')
            write_csv(ROOT / 'raw' / f'{args.output}_teeth.csv', measurements)
            (ROOT / 'CURRENT_WORK_STATE.json').write_text(json.dumps({'lane': 'X12-pulpy3d', 'milestone': 'PAIRED_ATLAS_RUNNING', 'latest_case': pid, 'paired_cases': sum((p.get('paired', False) for p in pairs)), 'aligned_cases': sum((p.get('aligned', False) for p in pairs)), 'measured_teeth': len(measurements), 'last_gate': 'R2 exact CT identity and FDI containment', 'next_operation': 'remaining cases; external comparator and material/heat port', 'elapsed_s': time.time() - started}, indent=2) + '\n')
            if args.limit and sum((p.get('paired', False) for p in pairs)) >= args.limit:
                break
    cost = {'wall_s': time.time() - started, 'cpu_s': time.process_time() - cpu, 'peak_RSS_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'fit_s': 0, 'GPU_s': 0, 'cases_considered': len(pairs), 'paired': sum((p.get('paired', False) for p in pairs)), 'aligned': sum((p.get('aligned', False) for p in pairs)), 'measured_teeth': len(measurements)}
    (ROOT / 'raw' / f'{args.output}_cost.json').write_text(json.dumps(cost, indent=2) + '\n')
    print(json.dumps(cost), flush=True)
if __name__ == '__main__':
    main()
