import json, time
import numpy as np
from scipy.spatial.transform import Rotation
from measure_r1 import HERE, tf2_io, dump, sha, exact_ray, side_measure, state
from measure_r2 import slab_control
from measure_r3 import matrix_control
from motion_envelope import interval, query

def surface_probe(lab, sp, m, ix, t, th):
    a = np.array(m['axis_xyz'])
    b = np.array(m['buccal_xyz'])
    base = np.array(m['base_xyz_mm'])
    pivot = base + m['pivot_height_mm'] * a
    k = np.cross(a, b)
    k /= np.linalg.norm(k)
    xyz = (np.array(ix) * sp)[:, ::-1]
    moved = pivot + Rotation.from_rotvec(np.deg2rad(th) * k).apply(xyz - pivot) + t * b
    margins = []
    missing = 0
    err = 0.0
    records = []
    for (p, src) in zip(moved, ix):
        ex = exact_ray(lab, p, b, sp)
        co = slab_control(lab, p, b, sp)
        if len(ex) != len(co) or any((e[2] != c[2] for (e, c) in zip(ex, co))):
            err = max(err, 999.0)
        else:
            err = max(err, max((abs(e[j] - c[j]) for (e, c) in zip(ex, co) for j in (0, 1))))
        point = {'source_voxel_zyx': src, 'moved_xyz_mm': p.tolist()}
        for (s, n) in ((1, 'buccal'), (-1, 'lingual')):
            x = side_measure(ex, m['tooth'], 1, s)
            c = side_measure(co, m['tooth'], 1, s)
            if x['status'] == c['status'] == 'MEASURED_ANNOTATION_SURROGATE':
                val = float(s * x['outer_x_mm'])
                err = max(err, abs(val - s * c['outer_x_mm']))
                margins.append(val)
                point[n] = {'margin_mm': val, 'status': 'MEASURED_ANNOTATION_SURROGATE'}
            else:
                missing += 1
                point[n] = {'status': 'UNDETERMINED', 'reason': x['status']}
        records.append(point)
    minimum = min(margins) if margins else None
    return {'points': records, 'min_margin_mm': minimum, 'missing_side_queries': missing, 'control_error_mm': err, 'conditional_feasible': bool(not missing and minimum is not None and (minimum >= 0.608)), 'anatomical_status': 'UNDETERMINED'}

def main():
    f = json.load(open(HERE / 'FROZEN_PREDICTIONS_R5.json'))
    assert sha(HERE / 'PREREG_R5.json') == f['prereg_sha256']
    models = [m for m in json.load(open(HERE / 'MOTION_MODELS.json')) if m['source'] == 'label' and m['tooth'] in (31, 32, 41, 42)]
    sample = {(r['case'], r['tooth']): r['sample_points'] for r in json.load(open(HERE / 'R4_SAMPLE_POINTS.json'))}
    rr = []
    errors = []
    marginerrors = []
    surfaceerrors = []
    angles = np.arange(-15, 15.0001, 0.25)
    start = time.monotonic()
    for case in sorted(set((m['case'] for m in models))):
        (lab, sp, _) = tf2_io.load(case)
        sp = np.array(sp)
        for m in [m for m in models if m['case'] == case]:
            (lo, hi) = interval(m, angles)
            (cl, ch) = matrix_control(m, angles)
            span = m['pivot_height_mm'] / 0.55
            lever = span - m['pivot_height_mm']
            a = np.array(m['axis_xyz'])
            b = np.array(m['buccal_xyz'])
            k = np.cross(a, b)
            k /= np.linalg.norm(k)
            tip = (span - m['pivot_height_mm']) * a
            bodily = surface_probe(lab, sp, m, sample[case, m['tooth']], 5.0, 0.0)
            surfaceerrors.append(bodily['control_error_mm'])
            item = {'case': case, 'tooth': m['tooth'], 'historical_5mm_bodily_probe': bodily, 'targets': []}
            for target in (1.0, 5.0):
                ts = target - lever * np.sin(np.deg2rad(angles))
                margins = np.minimum(ts - lo, hi - ts)
                ref = np.minimum(ts - cl, ch - ts)
                marginerrors.append(float(np.max(abs(margins - ref))))
                crown = np.array([(Rotation.from_rotvec(np.deg2rad(th) * k).apply(tip) - tip) @ b + t for (th, t) in zip(angles, ts)])
                errors.append(float(np.max(abs(crown - target))))
                order = sorted(range(len(angles)), key=lambda j: (-margins[j], abs(angles[j])))
                j = order[0]
                t = float(ts[j])
                th = float(angles[j])
                q = query(m, t, th, 0.608)
                probe = surface_probe(lab, sp, m, sample[case, m['tooth']], t, th)
                surfaceerrors.append(probe['control_error_mm'])
                item['targets'].append({'target_crown_advancement_mm': target, 'chosen_translation_mm': t, 'chosen_angle_deg': th, 'chosen_central_margin_mm': float(margins[j]), 'matrix_margin_mm': float(ref[j]), 'crown_target_error_mm': float(abs(crown[j] - target)), 'conditional_central_status': q['conditional_geometric_status'], 'moved_surface_probe': probe, 'anatomical_status': 'UNDETERMINED'})
            rr.append(item)
        state('R5_RUNNING', {'case': case, 'profiles': len(rr)}, 'Finish crown-target comparison; preserve geometry-only inverse candidates')
        print(case, 'profiles', len(rr), flush=True)
    dump(HERE / 'RAW_R5.json', rr)
    stats = {}
    for target in (1.0, 5.0):
        q = [q for r in rr for q in r['targets'] if q['target_crown_advancement_mm'] == target]
        stats[str(target)] = {'profiles': len(q), 'central_conditional_feasible': sum((q['chosen_central_margin_mm'] >= 0.608 for q in q)), 'surface_conditional_feasible': sum((q['moved_surface_probe']['conditional_feasible'] for q in q)), 'median_chosen_angle_deg': float(np.median([q['chosen_angle_deg'] for q in q])) if q else None}
    sm = {'round': 'R5', 'lower_incisor_profiles': len(models), 'crown_target_control_max_error_mm': max(errors, default=999.0), 'margin_control_max_error_mm': max(marginerrors, default=999.0), 'surface_ray_slab_max_error_mm': max(surfaceerrors, default=999.0), 'crown_gate': 'PASS' if max(errors, default=999.0) <= 1e-09 else 'FAIL', 'margin_gate': 'PASS' if max(marginerrors, default=999.0) <= 1e-09 else 'FAIL', 'surface_control_gate': 'PASS' if max(surfaceerrors, default=999.0) <= 1e-08 else 'FAIL', 'historical_5mm_bodily_surface_conditional_feasible': sum((r['historical_5mm_bodily_probe']['conditional_feasible'] for r in rr)), 'historical_5mm_bodily_negative_measured_margin': sum((r['historical_5mm_bodily_probe']['min_margin_mm'] is not None and r['historical_5mm_bodily_probe']['min_margin_mm'] < 0 for r in rr)), 'targets': stats, 'anatomical_accuracy': 'UNKNOWN', 'clinical_benefit': 'UNKNOWN', 'runtime_wall_s': time.monotonic() - start, 'external_referent': json.load(open(HERE / 'PREREG_R5.json'))['external_referent']}
    dump(HERE / 'SUMMARY_R5.json', sm)
    state('R5_ADJUDICATED', sm, 'Deliver one-command demo; whole-root calibrated anatomy and paired outcomes remain next construction')
    print(json.dumps(sm, indent=2))
if __name__ == '__main__':
    main()
