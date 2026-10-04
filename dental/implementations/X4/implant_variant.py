"""Conditional future-implant envelope inside finite straight grafts.

Geometry scenarios are not clinical implant recommendations. The default plane
is a superior-surface proxy, not a measured occlusal plane. A lab can supply a
three-point plane JSON; predictions are frozen before its physical measurement.
"""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '2'
import argparse, json, time, resource
from pathlib import Path
import numpy as np
from freeze import ROOT, DATA, write, sha, state, now, freeze
from planner import load, rigid_fit, transform, extract_rail, chord_errors, graft_mesh, surface_samples, stl, graft_levelset, planes
from constrained import constrained_plan, cut_feasible

def plane_from_three_points(points):
    p = np.asarray(points, float)
    if p.shape != (3, 3) or not np.isfinite(p).all():
        raise ValueError('plane requires three finite 3D points')
    n = np.cross(p[1] - p[0], p[2] - p[0])
    length = np.linalg.norm(n)
    if length < 1e-08:
        raise ValueError('collinear plane points: UNKNOWN')
    n /= length
    if n[2] < 0:
        n = -n
    return (n, float(p[0] @ n))

def envelope(nodes, radius, plane_normal, plane_offset, prosthetic_space=10.0, length=10.0, implant_radius=2.0):
    """Sufficient sampled containment bound, with angular Lipschitz padding.

Convexity of each cylinder/halfspace intersection places the max over implant
axial positions at its two ends. The 64 ring angles have covering radius pi/64;
all radial and halfspace functions are 1-Lipschitz, so add r_implant*pi/64.
Faceted donor mesh radius is reduced by cos(pi/48) conservatively.
"""
    a = np.asarray(plane_normal, float)
    a /= np.linalg.norm(a)
    ref = np.array([1.0, 0.0, 0.0]) if abs(a[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e = np.cross(a, ref)
    e /= np.linalg.norm(e)
    other = np.cross(a, e)
    theta = np.arange(64) * 2 * np.pi / 64
    (u, nn, _) = planes(nodes)
    records = []
    for i in range(len(u)):
        mid = (nodes[i] + nodes[i + 1]) / 2
        top = mid + (plane_offset - prosthetic_space - mid @ a) * a
        rings = []
        for s in (0.0, length):
            rings.append(top - s * a + implant_radius * (np.cos(theta)[:, None] * e + np.sin(theta)[:, None] * other))
        pp = np.r_[tuple(rings)]
        d = pp - nodes[i]
        r_eff = radius * np.cos(np.pi / 48)
        phi = np.maximum.reduce([np.linalg.norm(d - (d @ u[i])[:, None] * u[i], axis=1) - r_eff, -d @ nn[i], (pp - nodes[i + 1]) @ nn[i + 1]])
        upper = float(phi.max() + implant_radius * np.pi / 64 + 0.0001)
        records.append(dict(segment=i + 1, implant_top_mm=top.tolist(), implant_axis=(-a).tolist(), length_mm=length, radius_mm=implant_radius, containment_upper_bound_mm=upper, contained=upper <= 0, bound_scope='assumed circular donor and cut halfspaces only; no nerve/cortex/pedicle/fixation constraints'))
    return records

def freeze_r3():
    r = json.loads((ROOT / 'RESULTS_R2.json').read_text())
    freeze('R3', dict(prediction=['Pre', 'expert inferior rail', 'upper anterior Post surface plane proxy or explicitly supplied three-point plane'], evaluation=['full Post surface', 'same conditional implant envelope in every arm'], blind_reconstruction=False), 'joint graft radius/vertical placement and thickness-feasible segment search, constrained by future implant cylindrical envelope; freeze before surface comparison', dict(all_cases_available=True, median_primary_p95_mm_max=3.0, fraction_all_segment_envelopes_contained_min=0.8, median_surface_penalty_vs_shape_only_mm_max=2.0, all_miters_non_crossing=True, every_gate_fault_injection=True))
    p = ROOT / 'PREREG_R3.json'
    reg = json.loads(p.read_text())
    if 'implant_contract' not in reg:
        reg['implant_contract'] = dict(radii_mm=[6.0, 8.0, 10.0], vertical_offsets_mm=[-10.0, -8.0, -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0], implant_length_mm=10.0, implant_diameter_mm=4.0, prosthetic_space_mm=10.0, angular_samples=64, containment_padding_mm=float(2 * np.pi / 64), top_plane='q90(z) of anterior abs(theta)<60deg Post jaw proxy; anatomy UNKNOWN', candidate='min radius then abs offset among thickness-feasible plans containing all segment envelopes; fallback minimum envelope violation', straightness_repair='axis transport of one cross-section per segment; legacy R1/R2 independently phased rings had twisted longitudinal edges', straightness_tolerance_mm=1e-08, stl_export_tolerance_mm=0.0001, equal_information_control='for each supplied base plane and nodes, enumerate same radius/height domain independently with same envelope bound', surface_tradeoff_control='corrected straight shape-only circular graft of radius6mm on same expert inferior rail; not the invalid R2 legacy mesh', physical_measurement='NOT_PERFORMED; plane can later be supplied from 3 measured points; donor radius needs matched CT', external_referent_scope='only surface error is judged against published Post; implant plane and dimensions are scenario closures, not published measurement')
        write(p, reg)
        (ROOT / 'PREREG_R3.sha256').write_text(sha(p) + '  ' + p.name + '\n')
        decomp = json.loads((ROOT / 'DECOMPOSITION_R3.json').read_text())
        decomp['equations'].extend(['radial_j=(ring_j-node_i)-dot(ring_j-node_i,u_i)*u_i; next_ring=node_next+radial_j-u_i*dot(radial_j,n_next)/dot(u_i,n_next)', 'containment_upper=max(endpoint/ring64 cylinder-halfspace violation)+r_implant*pi/64+1e-4mm', 'effective circular radius=r_donor*cos(pi/48) for inside faceted-cylinder support'])
        decomp['leaves'].extend([dict(name='cross-section transport', status='DERIVED_UNDER_ASSUMPTIONS', basis='equal-radius cylinder intersection at bisector plane; identical radial vectors give parallel generator edges', stopping_argument='circular straight donor surrogate only; arbitrary donor cortex needs a new sweep operation'), dict(name='implant containment', status='DERIVED_UNDER_ASSUMPTIONS', basis='convex maximum on axial endpoints and disk perimeter; one-Lipschitz angular covering bound', stopping_argument='no cortex/nerve/pedicle/prosthetic clinical suitability'), dict(name='plane and sizes', status='CONSTITUTIVE_CLOSURE', basis='upper anterior Post q90 and declared 10mm length/4mm diameter/10mm prosthetic space/radii6,8,10', stopping_argument='not labelled occlusion or measured donor; replace with independent measurements')])
        write(ROOT / 'DECOMPOSITION_R3.json', decomp)
    (ROOT / 'HANDOFF_R2.md').write_text('# R2 decided\n\n' + r['outcome'] + '; ' + str(r['n_available']) + '/118 available.\n\nChanged target information to a supplied expert inferior curve and moved finite miter extents inside segment search. Actual pinned upstream OsteoOpt simplifier executed on same curves. Full OsteoOpt union BO not available. See RESULTS_R2.json, frozen predictions, raw hashes and source code.\n\nNext adds a conditional implant envelope to graft radius/height search. Plane, donor and implant dimensions remain explicitly uncalibrated scenario inputs.\n')
    state('R3_preregistered', 'freeze implant-driven plans before external surface evaluation', r['outcome'])

def predict_all(plane_path=None):
    reg = json.loads((ROOT / 'PREREG_R3.json').read_text())
    frozen = ROOT / 'FROZEN_PREDICTIONS_R3.json'
    if frozen.exists():
        return
    config = reg['implant_contract']
    records = []
    tick = time.perf_counter()
    source_files = ['implant_variant.py', 'planner.py', 'constrained.py', 'freeze.py', 'competitor.py']
    pinned = []
    for name in source_files:
        dst = ROOT / 'sources/R3_code' / name
        dst.parent.mkdir(exist_ok=True)
        dst.write_bytes((ROOT / name).read_bytes())
        pinned.append(dict(path=str(dst), sha256=sha(dst)))
    supplied = json.loads(Path(plane_path).read_text()) if plane_path else None
    for cid in reg['population']:
        out = DATA / 'R3' / cid
        out.mkdir(parents=True, exist_ok=True)
        start = time.perf_counter()
        row = dict(case=cid, round='R3', created_utc=now(), status='UNKNOWN', code_sha256=sha(__file__), prereg_sha256=sha(ROOT / 'PREREG_R3.json'))
        try:
            (pre, a, pmeta) = load(cid, 'Pre')
            (post, pa, pm) = load(cid, 'Post')
            (T, fit) = rigid_fit(a[:4000], pa[:30000])
            q = transform(post, np.linalg.inv(T))
            (curve, rail) = extract_rail(q, pre, a)
            (E, _) = chord_errors(curve)
            if supplied and cid in supplied:
                (n, offset) = plane_from_three_points(supplied[cid]['points_mm'])
                plane_source = 'supplied three-point plane, data contract; physical measurement provenance must be supplied'
            else:
                origin = np.array(rail['origin_proxy_mm'])
                theta = np.arctan2(q[:, 0] - origin[0], -(q[:, 1] - origin[1]))
                anterior = q[abs(theta) < np.pi / 3]
                if len(anterior) < 100:
                    raise ValueError('UNKNOWN_NO_SUPERIOR_PLANE_PROXY_SUPPORT')
                n = np.array([0.0, 0.0, 1.0])
                offset = float(np.quantile(anterior[:, 2], 0.9))
                plane_source = 'Post anterior superior-surface q90 proxy; not labelled occlusion'
            arms = []
            for radius in config['radii_mm']:
                base = curve.copy()
                base[:, 2] += radius - 6
                try:
                    (nodes, plan) = constrained_plan(base, E, radius)
                except ValueError:
                    continue
                for dz in config['vertical_offsets_mm']:
                    moved = nodes + np.array([0.0, 0.0, dz])
                    env = envelope(moved, radius, n, offset)
                    penalty = max((r['containment_upper_bound_mm'] for r in env))
                    arms.append(dict(radius_mm=radius, offset_mm=dz, nodes=moved, plan=plan, envelopes=env, violation_mm=penalty, all_contained=penalty <= 0))
            if not arms:
                raise ValueError('UNKNOWN_NO_RADIUS_PARTITION')
            feasible = [v for v in arms if v['all_contained']]
            selected = min(feasible, key=lambda v: (v['radius_mm'], abs(v['offset_mm']), v['offset_mm'])) if feasible else min(arms, key=lambda v: (v['violation_mm'], v['radius_mm'], abs(v['offset_mm'])))
            controlled = sorted(arms, key=lambda v: (0, v['radius_mm'], abs(v['offset_mm']), v['offset_mm']) if v['all_contained'] else (1, v['violation_mm'], v['radius_mm'], abs(v['offset_mm'])))[0]
            agree = bool(selected['radius_mm'] == controlled['radius_mm'] and selected['offset_mm'] == controlled['offset_mm'])
            (meshes, detail) = graft_mesh(selected['nodes'], selected['radius_mm'])
            p = out / 'prediction.npz'
            (base_nodes, base_plan) = constrained_plan(curve, E, 6.0)
            (base_meshes, base_detail) = graft_mesh(base_nodes, 6.0)
            np.savez_compressed(p, dp_nodes=selected['nodes'].astype('f4'), dp_points=surface_samples(meshes, 3000, 20261002 + int(cid)).astype('f4'), dp_points_dense=surface_samples(meshes, 12000, 20261002 + int(cid)).astype('f4'), curve=curve.astype('f8'), shape_only_nodes=base_nodes.astype('f4'), shape_only_points=surface_samples(base_meshes, 3000, 20261002 + int(cid)).astype('f4'), shape_only_points_dense=surface_samples(base_meshes, 12000, 20261002 + int(cid)).astype('f4'))
            for (si, mesh) in enumerate(meshes):
                stl(out / f'dp_segment_{si + 1}.stl', mesh)
            row.update(status='PREDICTED', rail=rail, dp=selected['plan'], grafts={'dp': detail, 'shape_only': base_detail}, arrays=dict(path=str(p), sha256=sha(p)), shape_only_envelopes=envelope(base_nodes, 6.0, n, offset), information='expert-assisted geometry; scenario plane and donor radius', inputs=[pmeta, pm], plane=dict(normal=n.tolist(), offset_mm=offset, source=plane_source, anatomical_occlusal_status='UNKNOWN'), envelopes=selected['envelopes'], all_envelopes_contained=selected['all_contained'], equal_information_control_agrees=agree, selected_radius_mm=selected['radius_mm'], selected_vertical_offset_mm=selected['offset_mm'], search_arms=[{k: v for (k, v) in arm.items() if k not in ('nodes', 'plan')} for arm in arms], registration=dict(T_Pre_to_Post=T.tolist(), fit=fit))
        except ValueError as exc:
            row['error'] = str(exc)
        row['prediction_seconds'] = time.perf_counter() - start
        path = out / 'plan.json'
        write(path, row)
        records.append(dict(case=cid, path=str(path), sha256=sha(path)))
        state('R3_predicted_' + cid, 'freeze all implant-driven predictions before surface evaluation')
        if len(records) % 20 == 0:
            print('R3 PREDICTED', len(records), flush=True)
    write(frozen, dict(round='R3', frozen_utc=now(), prereg_sha256=sha(ROOT / 'PREREG_R3.json'), predictor_code_sha256=sha(__file__), supplied_plane_input=dict(path=str(plane_path), sha256=sha(plane_path)) if plane_path else None, predictions=records, source_files=pinned, wall_seconds=time.perf_counter() - tick, peak_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    (ROOT / 'FROZEN_PREDICTIONS_R3.sha256').write_text(sha(frozen) + '\n')
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('stage', choices=['freeze', 'predict'])
    ap.add_argument('--plane-json')
    a = ap.parse_args()
    if a.stage == 'freeze':
        freeze_r3()
    else:
        predict_all(a.plane_json)
