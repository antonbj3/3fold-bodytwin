import json, time, gc
import numpy as np
from scipy import ndimage as ndi
from scipy.spatial.transform import Rotation
from measure_r1 import HERE, DATA, tf2_io, sha, dump, state, exact_ray, side_measure
from measure_r2 import slab_control
from motion_envelope import query

def main():
    f = json.load(open(HERE / 'FROZEN_PREDICTIONS_R4.json'))
    assert sha(HERE / 'PREREG_R4.json') == f['prereg_sha256']
    models = [m for m in json.load(open(HERE / 'MOTION_MODELS.json')) if m['source'] == 'label']
    motions = json.load(open(HERE / 'PREREG_R4.json'))['motions']
    raw = []
    case_meta = []
    maxerr = 0.0
    mismatches = 0
    start = time.monotonic()
    for case in sorted(set((m['case'] for m in models))):
        (lab, sp, hdr) = tf2_io.load(case)
        sp = np.array(sp)
        boxes = ndi.find_objects(lab, max_label=48)
        for m in [m for m in models if m['case'] == case]:
            tooth = m['tooth']
            sl = boxes[tooth - 1]
            mask = lab[sl] == tooth
            boundary = mask & ~ndi.binary_erosion(mask)
            ix = np.argwhere(boundary) + np.array([s.start for s in sl])
            xyz = (ix * sp)[:, ::-1]
            a = np.array(m['axis_xyz'])
            b = np.array(m['buccal_xyz'])
            base = np.array(m['base_xyz_mm'])
            h = (xyz - base) @ a
            keep = (h >= 2) & (h <= 8) & (h <= m['pivot_height_mm'])
            ix = ix[keep]
            xyz = xyz[keep]
            picks = np.linspace(0, len(ix) - 1, min(32, len(ix)), dtype=int) if len(ix) else np.array([], dtype=int)
            ix = ix[picks]
            xyz = xyz[picks]
            pivot = base + m['pivot_height_mm'] * a
            k = np.cross(a, b)
            k /= np.linalg.norm(k)
            case_meta.append({'case': case, 'tooth': tooth, 'available_boundary_voxels': int(keep.sum()), 'sample_points': ix.tolist(), 'frame_source': 'MOTION_MODELS.json'})
            for motion in motions:
                (t, theta) = (motion['t_mm'], motion['theta_deg'])
                rot = Rotation.from_rotvec(np.deg2rad(theta) * k)
                moved = pivot + rot.apply(xyz - pivot) + t * b
                points = []
                margins = []
                missing = 0
                for (i, p) in enumerate(moved):
                    ex = exact_ray(lab, p, b, sp)
                    ctrl = slab_control(lab, p, b, sp)
                    if len(ex) != len(ctrl) or any((a[2] != b[2] for (a, b) in zip(ex, ctrl))):
                        seq_error = 999.0
                        mismatches += 1
                    else:
                        seq_error = max((abs(a[j] - b[j]) for (a, b) in zip(ex, ctrl) for j in (0, 1)))
                    maxerr = max(maxerr, seq_error)
                    point = {'source_voxel_zyx': ix[i].tolist(), 'source_xyz_mm': xyz[i].tolist(), 'moved_xyz_mm': p.tolist(), 'new_height_from_apex_mm': float((p - base) @ a), 'control_sequence_error_mm': seq_error}
                    for (side, name) in ((1, 'buccal'), (-1, 'lingual')):
                        s = side_measure(ex, tooth, 2 if tooth < 30 else 1, side)
                        c = side_measure(ctrl, tooth, 2 if tooth < 30 else 1, side)
                        if s['status'] == c['status'] == 'MEASURED_ANNOTATION_SURROGATE':
                            margin = float(side * s['outer_x_mm'])
                            control_margin = float(side * c['outer_x_mm'])
                            error = abs(margin - control_margin)
                            maxerr = max(maxerr, error)
                            margins.append(margin)
                            point[name] = {'status': 'MEASURED_ANNOTATION_SURFACE_MARGIN', 'margin_mm': margin, 'control_margin_mm': control_margin, 'error_mm': error, 'old_tooth_intersections': s['tooth_intersections'], 'later_bone_islands': s['later_bone_islands']}
                        else:
                            missing += 1
                            point[name] = {'status': 'UNDETERMINED_WALL', 'candidate_status': s['status'], 'control_status': c['status']}
                    points.append(point)
                old = query(m, t, theta, 0.608)
                minmargin = min(margins) if margins else None
                within = bool(len(points) > 0 and missing == 0 and (minmargin is not None) and (minmargin >= 0.608))
                if within:
                    status = 'WITHIN_SAMPLED_SURFACE_CONSTRAINTS'
                elif minmargin is not None and minmargin < -0.608:
                    status = 'VIOLATES_SAMPLED_SURFACE_CONSTRAINTS'
                elif missing or not len(points):
                    status = 'UNDETERMINED_WALL_COVERAGE'
                else:
                    status = 'UNDETERMINED_EDGE_SCENARIO'
                raw.append({'case': case, 'tooth': tooth, 'tooth_type': m['tooth_type'], 'motion': motion, 'R3_query': old, 'R4_min_measured_margin_mm': minmargin, 'missing_side_queries': missing, 'points_n': len(points), 'conditional_status': status, 'R3_admission_survived': within, 'anatomical_status': 'UNDETERMINED', 'surface_points': points})
        (HERE / 'RAW_R4.jsonl').write_text(''.join((json.dumps(r, allow_nan=False) + '\n' for r in raw)))
        dump(HERE / 'R4_SAMPLE_POINTS.json', case_meta)
        state('R4_RUNNING', {'case': case, 'queries': len(raw), 'control_max_error_mm': maxerr}, 'Finish moved-surface comparison; deliver one-command demo with all failed gates')
        print(case, 'queries', len(raw), 's', round(time.monotonic() - start, 2), flush=True)
        del lab
        gc.collect()
    accepted = [r for r in raw if r['R3_query']['conditional_geometric_status'] == 'WITHIN_SECTION_CONSTRAINTS']
    survived = [r for r in accepted if r['R3_admission_survived']]
    ratio = len(survived) / len(accepted) if accepted else None
    from collections import Counter
    sm = {'round': 'R4', 'tooth_profiles': len(models), 'motion_queries': len(raw), 'surface_points': sum((r['points_n'] for r in raw)), 'ray_vs_slab_max_error_mm': maxerr, 'run_sequence_mismatches': mismatches, 'control_gate': 'PASS' if maxerr <= 1e-08 else 'FAIL', 'central_R3_admitted_queries': len(accepted), 'admissions_surviving_sampled_R4': len(survived), 'survival_fraction': ratio, 'survival_gate': 'PASS' if ratio is not None and ratio >= 0.95 else 'FAIL' if ratio is not None else 'UNKNOWN', 'R4_status_counts': dict(Counter((r['conditional_status'] for r in raw))), 'R3_admitted_but_negative_measured_surface_margin': sum((r['R4_min_measured_margin_mm'] is not None and r['R4_min_measured_margin_mm'] < 0 for r in accepted)), 'R3_admitted_with_missing_walls': sum((r['missing_side_queries'] > 0 for r in accepted)), 'control_injection': {'error_mm': 1.0, 'gate': 'FAIL'}, 'survival_injection': {'wrong_survival_fraction': 0.5, 'gate': 'FAIL'}, 'physical_accuracy': 'UNKNOWN', 'biological_accuracy': 'UNKNOWN', 'runtime_wall_s': time.monotonic() - start, 'external_referent': json.load(open(HERE / 'PREREG_R4.json'))['external_referent']}
    dump(HERE / 'SUMMARY_R4.json', sm)
    state('R4_ADJUDICATED', sm, 'Build complete demo, independent-review feedback and R5 measurement-ready whole-root acquisition plan')
    print(json.dumps(sm, indent=2))
if __name__ == '__main__':
    main()
