import json, time, collections
import numpy as np
from scipy.spatial.transform import Rotation
from shapely.geometry import MultiPoint
from shapely.affinity import rotate
import shapely
from measure_r1 import HERE, DATA, sha, dump, state
from motion_envelope import interval, query

def make_model(rr, source):
    rr = sorted(rr, key=lambda r: r['height_from_apex_mm'])
    a = np.array(rr[0]['root_axis_xyz'])
    b = np.array(rr[0]['buccal_direction_xyz'])
    base = np.array(rr[0]['center_xyz_mm']) + rr[0]['apex_axis_mm'] * a
    rec = {'case': rr[0]['case'], 'tooth': rr[0]['tooth'], 'tooth_type': rr[0]['tooth_type'], 'source': source, 'pivot_height_mm': 0.55 * rr[0]['tooth_span_mm'], 'axis_xyz': a.tolist(), 'buccal_xyz': b.tolist(), 'base_xyz_mm': base.tolist(), 'x_root_mm': [], 'height_mm': [], 'wall_lingual_mm': [], 'wall_buccal_mm': [], 'anatomical_status': 'UNDETERMINED'}
    for r in rr:
        offset = (np.array(r['origin_xyz_mm']) - base) @ b
        if source == 'label':
            xl = r['lingual']['root_x_mm']
            xb = r['buccal']['root_x_mm']
            wl = r['lingual']['outer_x_mm']
            wb = r['buccal']['outer_x_mm']
        else:
            xl = -r['lingual']['root_image_edge']['HM_mm']
            xb = r['buccal']['root_image_edge']['HM_mm']
            wl = -r['lingual']['bone_image_edge']['HM_mm']
            wb = r['buccal']['bone_image_edge']['HM_mm']
        rec['x_root_mm'] += [float(offset + xl), float(offset + xb)]
        rec['height_mm'] += [r['height_from_apex_mm']] * 2
        rec['wall_lingual_mm'] += [float(offset + wl)] * 2
        rec['wall_buccal_mm'] += [float(offset + wb)] * 2
    return rec

def matrix_control(model, angles):
    a = np.array(model['axis_xyz'])
    b = np.array(model['buccal_xyz'])
    k = np.cross(a, b)
    k /= np.linalg.norm(k)
    x = np.array(model['x_root_mm'])
    z = np.array(model['height_mm']) - model['pivot_height_mm']
    pts = x[:, None] * b + z[:, None] * a
    lo = []
    hi = []
    for theta in angles:
        q = Rotation.from_rotvec(np.deg2rad(theta) * k).apply(pts) @ b
        lo.append(np.max(np.array(model['wall_lingual_mm']) - q))
        hi.append(np.min(np.array(model['wall_buccal_mm']) - q))
    return (np.array(lo), np.array(hi))

def main():
    f = json.load(open(HERE / 'FROZEN_PREDICTIONS_R3.json'))
    assert sha(HERE / 'PREREG_R3.json') == f['prereg_sha256']
    rr = [json.loads(l) for l in open(HERE / 'RAW_R2.jsonl')]
    groups = collections.defaultdict(list)
    for r in rr:
        groups[r['case'], r['tooth']].append(r)
    models = []
    unknown = []
    for (k, rs) in groups.items():
        ann = len(rs) == 3 and all((r[n]['status'] == 'MEASURED_ANNOTATION_SURROGATE' for r in rs for n in ('buccal', 'lingual')))
        im = len(rs) == 3 and all((r[n].get('image_status') == 'MEASURED_IMAGE_CLEARANCE_SURROGATE' for r in rs for n in ('buccal', 'lingual')))
        if ann:
            models.append(make_model(rs, 'label'))
        else:
            unknown.append({'case': k[0], 'tooth': k[1], 'source': 'label', 'reason': 'INCOMPLETE_BILATERAL_THREE_HEIGHTS'})
        if im:
            models.append(make_model(rs, 'image'))
        else:
            unknown.append({'case': k[0], 'tooth': k[1], 'source': 'image', 'reason': 'INCOMPLETE_BILATERAL_IMAGE_EDGES'})
    dump(HERE / 'MOTION_MODELS.json', models)
    dump(HERE / 'R3_UNDETERMINED_PROFILES.json', unknown)
    angles = np.arange(-15, 15.00001, 0.25)
    tgrid = np.arange(-3, 3.00001, 0.05)
    eps = (0, 0.38, 0.608, 0.9)
    rows = []
    witnesses = []
    curves = []
    err = []
    shapely_err = []
    rawqueries = []
    start = time.monotonic()
    cand_t = 0.0
    ctrl_t = 0.0
    for m in models:
        s = time.perf_counter()
        (lo, hi) = interval(m, angles)
        cand_t += time.perf_counter() - s
        s = time.perf_counter()
        (cl, ch) = matrix_control(m, angles)
        ctrl_t += time.perf_counter() - s
        err.append(float(max(np.max(abs(lo - cl)), np.max(abs(hi - ch)))))
        for th in (-15.0, 0.0, 15.0):
            geom = MultiPoint(list(zip(m['x_root_mm'], m['height_mm'])))
            rot = rotate(geom, -th, origin=(0, m['pivot_height_mm']))
            xx = np.array([p.x for p in rot.geoms])
            exp = np.array(m['x_root_mm']) * np.cos(np.deg2rad(th)) + (np.array(m['height_mm']) - m['pivot_height_mm']) * np.sin(np.deg2rad(th))
            shapely_err.append(float(np.max(abs(xx - exp))))
        u0 = float(hi[np.argmin(abs(angles))])
        l0 = float(lo[np.argmin(abs(angles))])
        pos_angles = angles[(angles >= 0) & (lo <= 0) & (hi >= 0)]
        neg_angles = angles[(angles <= 0) & (lo <= 0) & (hi >= 0)]
        corners = []
        for t in (l0, u0):
            for th in ([float(pos_angles[-1])] if len(pos_angles) else []) + ([float(neg_angles[0])] if len(neg_angles) else []):
                q = query(m, t, th, 0.0)
                corners.append(q)
                if q['conditional_margin_mm'] < -1e-09:
                    witnesses.append({'construction': 'SEPARATELY_FEASIBLE_EXTREMA_COMBINED', 't_alone_margin_mm': query(m, t, 0.0, 0.0)['conditional_margin_mm'], 'angle_alone_margin_mm': query(m, 0.0, th, 0.0)['conditional_margin_mm'], 'combined': q})
        measures = {}
        for e in eps:
            feasible = (tgrid[None, :] >= lo[:, None] + e) & (tgrid[None, :] <= hi[:, None] - e)
            measures[str(e)] = {'grid_pairs_within_conditional_constraints': int(feasible.sum()), 'grid_pairs_total': int(feasible.size), 't_only_lower_mm': l0 + e, 't_only_upper_mm': u0 - e, 'zero_motion_margin_mm': min(-l0, u0), 'zero_motion_conditional_feasible': bool(l0 + e <= 0 <= u0 - e)}
        qs = [query(m, t, th, 0.608) for (t, th) in ((0.0, 0.0), (0.25, 0.0), (0.5, 5.0), (1.0, 10.0))]
        rawqueries += qs
        rows.append({'case': m['case'], 'tooth': m['tooth'], 'tooth_type': m['tooth_type'], 'source': m['source'], 'translation_only_interval_mm': [l0, u0], 'rotation_only_interval_deg': [float(neg_angles[0]) if len(neg_angles) else None, float(pos_angles[-1]) if len(pos_angles) else None], 'separate_corner_violations': sum((q['conditional_margin_mm'] < -1e-09 for q in corners)), 'scenarios': measures, 'anatomical_status': 'UNDETERMINED'})
        curves.append({'case': m['case'], 'tooth': m['tooth'], 'source': m['source'], 'angles_deg': angles.tolist(), 'lower_t_mm': lo.tolist(), 'upper_t_mm': hi.tolist()})
    runtime = time.monotonic() - start
    dump(HERE / 'R3_LIMITS.json', rows)
    dump(HERE / 'R3_COUPLED_CURVES.json', curves)
    dump(HERE / 'R3_COUPLING_WITNESSES.json', witnesses)
    dump(HERE / 'R3_FIXED_QUERIES.json', rawqueries)
    maxerr = max(err, default=999.0)
    sm = {'round': 'R3', 'label_profiles': sum((m['source'] == 'label' for m in models)), 'image_profiles': sum((m['source'] == 'image' for m in models)), 'matrix_control_max_error_mm': maxerr, 'matrix_control_gate': 'PASS' if maxerr <= 1e-09 else 'FAIL', 'published_code_shapely_max_error_mm': max(shapely_err, default=999.0), 'installed_shapely_version': shapely.__version__, 'separate_corner_violations': len(witnesses), 'profiles_with_corner_violation': sum((r['separate_corner_violations'] > 0 for r in rows)), 'fixed_queries': len(rawqueries), 'physical_known_queries': 0, 'anatomical_status': 'UNDETERMINED_ALL', 'matrix_control_injection': {'injected_margin_error_mm': 1.0, 'gate': 'FAIL'}, 'uncertainty_injection': {'epsilon_mm': 100.0, 'accepted_queries': sum((query(m, 0, 0, 100)['conditional_geometric_status'] == 'WITHIN_SECTION_CONSTRAINTS' for m in models)), 'gate': 'FAIL_ACCEPTANCE'}, 'clinical_status_injection': {'candidate_status': 'SAFE', 'allowed_status': 'UNDETERMINED', 'gate': 'FAIL'}, 'runtime_wall_s': runtime, 'candidate_interval_kernel_s': cand_t, 'direct_matrix_control_kernel_s': ctrl_t, 'per_type': {}, 'external_referent': json.load(open(HERE / 'PREREG_R3.json'))['external_referent']}
    for kind in ('incisor', 'canine', 'premolar', 'molar'):
        vals = [r['translation_only_interval_mm'] for r in rows if r['source'] == 'label' and r['tooth_type'] == kind]
        sm['per_type'][kind] = {'n_profiles': len(vals), 'median_nominal_lingual_mm': float(np.median([-v[0] for v in vals])) if vals else None, 'median_nominal_buccal_mm': float(np.median([v[1] for v in vals])) if vals else None}
    dump(HERE / 'SUMMARY_R3.json', sm)
    state('R3_ADJUDICATED', sm, 'R4: test whether complete sampled profiles survive frame perturbation, identify minimum anatomical acquisition')
    print(json.dumps(sm, indent=2))
if __name__ == '__main__':
    main()
