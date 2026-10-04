"""Grid-interval topology and separately vetoed gray-profile edge observations."""
import json, time, sys, gc
import numpy as np
from scipy import ndimage as ndi
from measure_r1 import HERE, DATA, tf2_io, dump, sha, state, side_measure

def slab_control(lab, o_xyz, b_xyz, sp):
    """Enumerate overcomplete neighboring cubes; intersect each cube independently."""
    o = o_xyz[::-1] / sp
    v = b_xyz[::-1] / sp
    n = int(np.ceil(36 * np.max(np.abs(v)) / 0.5)) + 1
    centers = np.floor(o + np.linspace(-18, 18, n)[:, None] * v + 0.5).astype(int)
    neighbors = np.array(np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1])).reshape(3, -1).T
    cubes = np.unique((centers[:, None, :] + neighbors[None, :, :]).reshape(-1, 3), axis=0)
    enter = np.full(len(cubes), -18.0)
    leave = np.full(len(cubes), 18.0)
    for k in range(3):
        if abs(v[k]) < 1e-10:
            good = (o[k] >= cubes[:, k] - 0.5) & (o[k] < cubes[:, k] + 0.5)
            enter[~good] = 99
            leave[~good] = -99
        else:
            a = (cubes[:, k] - 0.5 - o[k]) / v[k]
            c = (cubes[:, k] + 0.5 - o[k]) / v[k]
            enter = np.maximum(enter, np.minimum(a, c))
            leave = np.minimum(leave, np.maximum(a, c))
    good = leave - enter > 1e-10
    cubes = cubes[good]
    enter = enter[good]
    leave = leave[good]
    ix = np.argsort(enter)
    out = []
    for k in ix:
        (a, c) = (float(enter[k]), float(leave[k]))
        ix3 = cubes[k]
        label = int(lab[tuple(ix3)]) if np.all(ix3 >= 0) and np.all(ix3 < np.array(lab.shape)) else -1
        if out and label == out[-1][2] and (abs(a - out[-1][1]) < 1e-08):
            out[-1][1] = c
        else:
            out.append([a, c, label])
    return out

def edge(x, y, label_edge):
    d = x - label_edge
    pin = (d >= -1.2) & (d <= -0.6)
    pout = (d >= 0.6) & (d <= 1.2)
    if min(sum(pin), sum(pout)) < 5:
        return {'status': 'UNDETERMINED_WINDOW'}
    inside = float(np.median(y[pin]))
    outside = float(np.median(y[pout]))
    contrast = inside - outside
    rec = {'contrast_gray': contrast, 'inside_gray': inside, 'outside_gray': outside, 'label_edge_mm': label_edge}
    if contrast < 200:
        rec['status'] = 'UNDETERMINED_CONTRAST'
        return rec
    hm = 0.5 * (inside + outside)
    cross = np.flatnonzero((y[:-1] >= hm) & (y[1:] < hm))
    loc = [float(x[k] + (hm - y[k]) * (x[k + 1] - x[k]) / (y[k + 1] - y[k])) for k in cross if abs(0.5 * (x[k] + x[k + 1]) - label_edge) <= 0.6]
    rec['crossings'] = loc
    if len(loc) != 1:
        rec['status'] = 'UNDETERMINED_MULTIPLE_OR_NO_CROSSING'
        return rec
    s = loc[0]
    smooth = ndi.gaussian_filter1d(y, 0.1 / (x[1] - x[0]))
    grad = np.gradient(smooth, x)
    win = np.flatnonzero(abs(d) <= 0.6)
    gm = float(x[win[np.argmin(grad[win])]])
    rec.update(status='MEASURED_IMAGE_EDGE_SURROGATE', HM_mm=s, gradient_mm=gm, offset_mm=s - label_edge, HM_gradient_discrepancy_mm=abs(s - gm))
    return rec

def main():
    f = json.load(open(HERE / 'FROZEN_PREDICTIONS_R2.json'))
    assert sha(HERE / 'PREREG_R2.json') == f['prereg_sha256']
    rr = [json.loads(l) for l in open(HERE / 'RAW_R1.jsonl')]
    out = []
    start = time.monotonic()
    diff = []
    diffrec = []
    for case in sorted(set((r['case'] for r in rr))):
        (lab, sp, hdr) = tf2_io.load(case)
        sp = np.array(sp)
        for r in [r for r in rr if r['case'] == case and r['status'] == 'MEASURED_SECTION']:
            slab = slab_control(lab, np.array(r['origin_xyz_mm']), np.array(r['buccal_direction_xyz']), sp)
            ex = r['exact_runs']
            if len(slab) != len(ex) or any((a[2] != b[2] for (a, b) in zip(slab, ex))):
                err = 999.0
                detail = 'RUN_SEQUENCE_MISMATCH'
            else:
                err = max((abs(a[j] - b[j]) for (a, b) in zip(slab, ex) for j in (0, 1)))
                detail = 'SAME_LABEL_SEQUENCE'
            diff.append(err)
            diffrec.append({'case': case, 'tooth': r['tooth'], 'h_mm': r['height_from_apex_mm'], 'max_error_mm': err, 'status': detail, 'slab_runs': slab})
            p = np.load(r['profile_path'])
            x = p['x_mm'].astype(float)
            y = p['gray'].astype(float)
            item = {k: r[k] for k in ('case', 'tooth', 'tooth_type', 'jaw', 'height_from_apex_mm', 'root_axis_xyz', 'buccal_direction_xyz', 'origin_xyz_mm', 'center_xyz_mm', 'apex_axis_mm', 'tooth_span_mm', 'profile_path')}
            item['topology_control_error_mm'] = err
            item['exact_runs'] = ex
            for (side, name) in ((1, 'buccal'), (-1, 'lingual')):
                m = side_measure(ex, r['tooth'], r['jaw'], side)
                m['slab_control'] = side_measure(slab, r['tooth'], r['jaw'], side)
                if m['status'] == 'MEASURED_ANNOTATION_SURROGATE':
                    if side == 1:
                        (xx, yy) = (x, y)
                    else:
                        (xx, yy) = (-x[::-1], y[::-1])
                    re = edge(xx, yy, side * m['root_x_mm'])
                    be = edge(xx, yy, side * m['outer_x_mm'])
                    if m['labelled_bone_run_mm'] < 0.6:
                        be = {'status': 'UNDETERMINED_THIN_BONE_NO_PLATEAU', 'labelled_bone_run_mm': m['labelled_bone_run_mm']}
                    m.update(root_image_edge=re, bone_image_edge=be)
                    if re['status'] == be['status'] == 'MEASURED_IMAGE_EDGE_SURROGATE':
                        m['image_clearance_mm'] = be['HM_mm'] - re['HM_mm']
                        m['image_minus_label_mm'] = m['image_clearance_mm'] - m['clearance_mm']
                        m['image_status'] = 'MEASURED_IMAGE_CLEARANCE_SURROGATE'
                        m['image_discrepancy_max_mm'] = max(re['HM_gradient_discrepancy_mm'], be['HM_gradient_discrepancy_mm'])
                    else:
                        m['image_status'] = 'UNDETERMINED_EDGE_IDENTIFICATION'
                item[name] = m
            out.append(item)
        state('R2_RUNNING', {'case': case, 'lines': len(out)}, 'Finish image-edge identifiability and exact slab control')
        print(case, 'lines', len(out), flush=True)
        del lab
        gc.collect()
    (HERE / 'RAW_R2.jsonl').write_text(''.join((json.dumps(r, allow_nan=False) + '\n' for r in out)))
    dump(HERE / 'R2_SLAB_CONTROL.json', diffrec)
    vs = [(r, n, r[n]) for r in out for n in ('buccal', 'lingual') if r[n].get('image_status') == 'MEASURED_IMAGE_CLEARANCE_SURROGATE']
    edges = [m[e] for (_, _, m) in vs for e in ('root_image_edge', 'bone_image_edge')]
    all_edges = [m[e] for r in out for n in ('buccal', 'lingual') for m in [r[n]] for e in ('root_image_edge', 'bone_image_edge') if m.get(e, {}).get('status') == 'MEASURED_IMAGE_EDGE_SURROGATE']
    counts = {}
    for r in out:
        if all((r[n].get('image_status') == 'MEASURED_IMAGE_CLEARANCE_SURROGATE' for n in ('buccal', 'lingual'))):
            counts[r['case'], r['tooth']] = counts.get((r['case'], r['tooth']), 0) + 1
    covered = [k for (k, v) in counts.items() if v >= 2]
    p95 = float(np.quantile([e['HM_gradient_discrepancy_mm'] for e in all_edges], 0.95)) if all_edges else None
    from collections import Counter
    statuses = Counter((m.get('image_status', m['status']) for r in out for n in ('buccal', 'lingual') for m in [r[n]]))
    offsets = [m['image_minus_label_mm'] for (_, _, m) in vs]
    sm = {'round': 'R2', 'exact_interval_control_max_error_mm': max(diff), 'exact_control_gate': 'PASS' if max(diff) <= 1e-08 else 'FAIL', 'exact_control_injection': {'error_mm': 1.0, 'gate': 'FAIL'}, 'image_clearance_rays': len(vs), 'identified_edges_all': len(all_edges), 'bilateral_teeth_two_heights': len(covered), 'covered_cases': len(set((k[0] for k in covered))), 'image_coverage_gate': 'PASS' if len(covered) >= 20 and len(set((k[0] for k in covered))) >= 4 else 'UNKNOWN', 'HM_gradient_p95_mm': p95, 'HM_gradient_gate': 'PASS' if p95 is not None and p95 <= 0.3 else 'FAIL', 'image_minus_label_median_mm': float(np.median(offsets)) if offsets else None, 'image_minus_label_p10_p90_mm': np.quantile(offsets, [0.1, 0.9]).tolist() if offsets else None, 'statuses': dict(statuses), 'physical_anatomy_accuracy': 'UNKNOWN', 'biological_accuracy': 'UNKNOWN', 'runtime_wall_s': time.monotonic() - start, 'external_referent': json.load(open(HERE / 'PREREG_R2.json'))['external_referent']}
    dump(HERE / 'SUMMARY_R2.json', sm)
    state('R2_ADJUDICATED', sm, 'R3: coupled exact translation/rotation limits and minimum edge measurement precision; do not promote image to anatomy truth')
    print(json.dumps(sm, indent=2))
if __name__ == '__main__':
    main()
