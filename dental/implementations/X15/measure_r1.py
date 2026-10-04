"""Local TF2 annotation-envelope measurements; no clinical prediction."""
from dental_release.paths import expand as _release_expand
import os, sys, json, hashlib, time, gc, zipfile
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
from scipy import ndimage as ndi
HERE = Path(__file__).resolve().parent
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X15-bone-dehiscence'))
DENT = Path(os.environ.get('DENTAL_PROJECT_ROOT', str(HERE.parents[1]))).resolve()
sys.path.insert(0, str(DENT / 'cells/geometry'))
import tf2_io

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def dump(p, o):
    Path(p).write_text(json.dumps(o, indent=2, allow_nan=False) + '\n')

def unit(v):
    return v / np.linalg.norm(v)

def state(status, gate, next_op):
    dump(HERE / 'CURRENT_WORK_STATE.json', {'lane': 'X15-bone-dehiscence', 'status': status, 'latest_gate': gate, 'next_operation': next_op, 'updated_utc': datetime.now(timezone.utc).isoformat(), 'review_state': 'PENDING_INDEPENDENT_REVIEW'})

def exact_ray(lab, origin, b, sp, extent=18):
    """Independent scalar grid-face traversal; constant labels between breakpoints."""
    o = origin[::-1] / sp
    v = b[::-1] / sp
    stops = [-extent, extent]
    for k in range(3):
        if abs(v[k]) < 1e-10:
            continue
        (lo, hi) = sorted([o[k] - extent * v[k], o[k] + extent * v[k]])
        faces = np.arange(np.ceil(lo - 0.5), np.floor(hi - 0.5) + 1) + 0.5
        stops.extend(((faces - o[k]) / v[k]).tolist())
    stops = np.unique(np.clip(stops, -extent, extent))
    out = []
    for (a, c) in zip(stops[:-1], stops[1:]):
        ix = np.floor(o + 0.5 * (a + c) * v + 0.5).astype(int)
        label = int(lab[tuple(ix)]) if np.all(ix >= 0) and np.all(ix < lab.shape) else -1
        if out and out[-1][2] == label:
            out[-1][1] = float(c)
        else:
            out.append([float(a), float(c), label])
    return out

def sample_runs(x, labels):
    changes = np.r_[0, np.flatnonzero(labels[1:] != labels[:-1]) + 1, len(labels)]
    runs = []
    for (i, j) in zip(changes[:-1], changes[1:]):
        a = x[i] - 0.5 * (x[1] - x[0])
        b = x[j - 1] + 0.5 * (x[1] - x[0])
        runs.append([float(a), float(b), int(labels[i])])
    return runs

def side_measure(runs, tooth, jaw, side):
    tr = [r for r in runs if r[2] == tooth]
    if not tr:
        return {'status': 'UNDETERMINED_NO_ROOT_ON_LINE'}
    root = max((r[1] for r in tr)) if side == 1 else min((r[0] for r in tr))

    def directed(r):
        return (r[0] - root, r[1] - root) if side == 1 else (root - r[1], root - r[0])
    after = sorted([(max(0, a), b, r[2]) for r in runs for (a, b) in [directed(r)] if b > 1e-08], key=lambda r: r[0])
    bones = [r for r in after if r[2] == jaw]
    if not bones:
        return {'status': 'UNDETERMINED_NO_ADJACENT_LABEL', 'root_x_mm': root}
    (start, end, _) = bones[0]
    if start > 1.2:
        return {'status': 'UNDETERMINED_NONLOCAL_BONE', 'root_x_mm': root, 'gap_to_bone_mm': start}
    obstacle = [r[2] for r in after if r[0] < end - 1e-08 and r[2] not in (0, tooth, jaw)]
    if obstacle:
        return {'status': 'UNDETERMINED_INTERVENING_STRUCTURE', 'obstacle_labels': obstacle, 'root_x_mm': root}
    if end >= 17.95:
        return {'status': 'UNDETERMINED_RAY_END'}
    return {'status': 'MEASURED_ANNOTATION_SURROGATE', 'root_x_mm': root, 'outer_x_mm': root + side * end, 'clearance_mm': end, 'labelled_bone_run_mm': end - start, 'gap_to_bone_mm': start, 'later_bone_islands': len(bones) - 1, 'tooth_intersections': len(tr), 'anatomy_truth': 'UNKNOWN', 'dehiscence': 'UNDETERMINED_NO_CEJ_CREST_OR_OUTCOME'}

def process(case):
    start = time.monotonic()
    (lab, sp, lh) = tf2_io.load(case)
    (image, isp, ih) = tf2_io.load(case, 'image')
    assert image.shape == lab.shape and sp == isp
    assert all((lh.get(k) == ih.get(k) for k in ('TransformMatrix', 'Offset', 'CenterOfRotation')))
    sp = np.array(sp)
    boxes = ndi.find_objects(lab, max_label=48)
    teeth = {}
    for t in [q * 10 + i for q in range(1, 5) for i in range(1, 8)]:
        sl = boxes[t - 1]
        if sl is None:
            continue
        ij = np.argwhere(lab[sl] == t) + np.array([s.start for s in sl])
        if len(ij) < 2000:
            continue
        xyz = (ij * sp)[:, ::-1]
        c = xyz.mean(0)
        fov = float(np.min(np.minimum(ij.min(0) * sp, (np.array(lab.shape) - 1 - ij.max(0)) * sp)))
        teeth[t] = {'p': xyz, 'c': c, 'fov': fov, 'bbox': sl, 'nvox': len(ij)}
    if min(sum((t < 30 for t in teeth)), sum((t >= 30 for t in teeth))) < 3:
        return ([], {'case': case, 'status': 'UNDETERMINED_FRAME'}, [])
    up = np.mean([v['c'] for (t, v) in teeth.items() if t < 30], axis=0) - np.mean([v['c'] for (t, v) in teeth.items() if t >= 30], axis=0)
    up = unit(up)
    manifest = []
    rows = []
    for jaw in (1, 2):
        jt = [t for t in teeth if (t < 30) == (jaw == 2)]
        arch = np.mean([teeth[t]['c'] for t in jt], axis=0)
        for t in jt:
            v = teeth[t]
            p = v['p']
            c = v['c']
            (_, _, vh) = np.linalg.svd(p - c, full_matrices=False)
            a = vh[0]
            crown = up if jaw == 1 else -up
            if a @ crown < 0:
                a = -a
            b = c - arch - (c - arch) @ a * a
            if np.linalg.norm(b) < 1e-05:
                continue
            b = unit(b)
            ax = (p - c) @ a
            apex = float(np.quantile(ax, 0.005))
            tip = float(np.quantile(ax, 0.995))
            span = tip - apex
            d = DATA / case
            d.mkdir(parents=True, exist_ok=True)
            sl = tuple((slice(max(0, s.start - 2), min(n, s.stop + 2)) for (s, n) in zip(v['bbox'], lab.shape)))
            crop = d / f'tooth_{t}_labels.npz'
            np.savez_compressed(crop, label=lab[sl], offset_zyx=np.array([s.start for s in sl]), spacing_zyx=sp)
            manifest.append({'path': str(crop), 'sha256': sha(crop), 'bytes': crop.stat().st_size, 'kind': 'local_tooth_label_crop'})
            for h in (2.0, 5.0, 8.0):
                rec = {'case': case, 'tooth': t, 'jaw': jaw, 'tooth_type': 'incisor' if t % 10 <= 2 else 'canine' if t % 10 == 3 else 'premolar' if t % 10 <= 5 else 'molar', 'height_from_apex_mm': h, 'tooth_span_mm': span, 'root_axis_xyz': a.tolist(), 'buccal_direction_xyz': b.tolist(), 'center_xyz_mm': c.tolist(), 'apex_axis_mm': apex, 'fov_margin_mm': v['fov'], 'coordinate_semantics': 'scan-axis physical xyz, not patient registration'}
                if h > 0.55 * span:
                    rec['status'] = 'UNDETERMINED_OUTSIDE_ROOT_SURROGATE'
                    rows.append(rec)
                    continue
                if v['fov'] < 1.5:
                    rec['status'] = 'UNDETERMINED_FOV'
                    rows.append(rec)
                    continue
                band = p[np.abs(ax - (apex + h)) <= 0.3]
                if len(band) < 10:
                    rec['status'] = 'UNDETERMINED_SPARSE_SECTION'
                    rows.append(rec)
                    continue
                origin = band.mean(0)
                origin += a * (apex + h - (origin - c) @ a)
                x = np.arange(-18, 18.00001, 0.05)
                coords = (origin[:, None] + b[:, None] * x)[::-1] / sp[:, None]
                labels = ndi.map_coordinates(lab, coords, order=0, mode='constant', cval=255).astype(int)
                vals = ndi.map_coordinates(image, coords, order=1, mode='constant', cval=0)
                srun = sample_runs(x, labels)
                erun = exact_ray(lab, origin, b, sp)
                rec.update(origin_xyz_mm=origin.tolist(), status='MEASURED_SECTION', sample_runs=srun, exact_runs=erun)
                for (side, name) in ((1, 'buccal'), (-1, 'lingual')):
                    m = side_measure(srun, t, jaw, side)
                    e = side_measure(erun, t, jaw, side)
                    m['exact_control'] = e
                    if m['status'] == e['status'] == 'MEASURED_ANNOTATION_SURROGATE':
                        m['control_error_mm'] = abs(m['clearance_mm'] - e['clearance_mm'])
                    else:
                        m['control_error_mm'] = None
                    rec[name] = m
                pf = d / f'profile_{t}_{int(h)}.npz'
                np.savez_compressed(pf, x_mm=x.astype(np.float32), gray=vals.astype(np.float32), labels=labels.astype(np.uint8), origin_xyz_mm=origin, buccal_xyz=b, axis_xyz=a, spacing_zyx_mm=sp)
                manifest.append({'path': str(pf), 'sha256': sha(pf), 'bytes': pf.stat().st_size, 'kind': 'label_and_image_ray'})
                rec['profile_path'] = str(pf)
                rows.append(rec)
    with zipfile.ZipFile(tf2_io.ZIP) as z:
        members = []
        for n in [f'{tf2_io.ROOT}/labelsTr/{case}.mha', f'{tf2_io.ROOT}/imagesTr/{case}_0000.mha']:
            zi = z.getinfo(n)
            members.append({'member': n, 'bytes': zi.file_size, 'crc32': zi.CRC})
    prov = {'case': case, 'headers': {'label': lh, 'image': ih}, 'members': members, 'elapsed_s': time.monotonic() - start, 'status': 'MEASURED'}
    return (rows, prov, manifest)

def main():
    frozen = json.load(open(HERE / 'FROZEN_PREDICTIONS_R1.json'))
    assert sha(HERE / 'PREREG_R1.json') == frozen['prereg_sha256']
    DATA.mkdir(parents=True, exist_ok=True)
    rows = []
    provs = []
    files = []
    start = time.monotonic()
    for case in json.load(open(HERE / 'PREREG_R1.json'))['selection']['cases']:
        (r, p, m) = process(case)
        rows += r
        provs.append(p)
        files += m
        (HERE / 'RAW_R1.jsonl').write_text(''.join((json.dumps(v, allow_nan=False) + '\n' for v in rows)))
        dump(HERE / 'DATA_MANIFEST_R1.json', {'files': files, 'provenance': provs, 'tf2_io_sha256': sha(DENT / 'cells/geometry/tf2_io.py'), 'prereg_sha256': frozen['prereg_sha256']})
        state('R1_RUNNING', {'case': case, 'sections': len(r)}, 'Finish label rays then adjudicate; freeze R2 before image-edge extraction')
        print(case, len(r), 'sections', round(p.get('elapsed_s', 0), 2), 's', flush=True)
        gc.collect()
    good = [r for r in rows if r['status'] == 'MEASURED_SECTION']
    valid = [(r, n, r[n]) for r in good for n in ('buccal', 'lingual') if r[n]['status'] == 'MEASURED_ANNOTATION_SURROGATE']
    errs = [v['control_error_mm'] for (_, _, v) in valid if v['control_error_mm'] is not None]
    counts = {}
    for r in good:
        if all((r[n]['status'] == 'MEASURED_ANNOTATION_SURROGATE' for n in ('buccal', 'lingual'))):
            counts[r['case'], r['tooth']] = counts.get((r['case'], r['tooth']), 0) + 1
    covered = [k for (k, v) in counts.items() if v >= 2]
    maxerr = max(errs, default=999.0)
    injerr = abs(1.0)
    sm = {'round': 'R1', 'sections_all': len(rows), 'sections_measured': len(good), 'valid_side_rays': len(valid), 'bilateral_teeth_two_heights': len(covered), 'covered_cases': len(set((k[0] for k in covered))), 'control_max_error_mm': maxerr, 'control_gate': 'PASS' if maxerr <= 0.051 else 'FAIL', 'control_injection': {'injected_error_mm': injerr, 'gate': 'FAIL' if injerr > 0.051 else 'PASS'}, 'coverage_gate': 'PASS' if len(covered) >= 30 and len(set((k[0] for k in covered))) >= 4 else 'UNKNOWN', 'clearance_median_mm': float(np.median([v['clearance_mm'] for (_, _, v) in valid])) if valid else None, 'fixed_rule_2mm_overruns_valid_rays': sum((v['clearance_mm'] < 2 for (_, _, v) in valid)), 'fixed_rule_2mm_semantics': 'Research comparator only; literal handbook rule differs and relates to crown advancement', 'uncertainty_scenarios': {}, 'biological_gate': 'FAIL_STATIC_EQUIVALENCE; OUR_ACCURACY_UNKNOWN', 'external_referent': json.load(open(HERE / 'PREREG_R1.json'))['external_referent'], 'per_type': {}, 'runtime_wall_s': time.monotonic() - start, 'n_fits': 0, 'bytes_saved': sum((v['bytes'] for v in files))}
    for e in (0, 0.38, 0.608, 0.9):
        sm['uncertainty_scenarios'][str(e)] = {'epsilon_mm': e, 'zero_motion_unresolved_valid_rays': sum((v['clearance_mm'] <= e for (_, _, v) in valid)), 'condition': 'Sensitivity scenario, no confidence interpretation'}
    for kind in ('incisor', 'canine', 'premolar', 'molar'):
        vs = [v for (r, n, v) in valid if r['tooth_type'] == kind]
        sm['per_type'][kind] = {'rays': len(vs), 'median_clearance_mm': float(np.median([v['clearance_mm'] for v in vs])) if vs else None}
    dump(HERE / 'SUMMARY_R1.json', sm)
    state('R1_ADJUDICATED', sm, 'R2: intensity edges replace annotation convention; preserve failed biological gate')
    print(json.dumps(sm, indent=2))
if __name__ == '__main__':
    main()
