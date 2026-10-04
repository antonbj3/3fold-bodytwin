"""R3: same measured crown, existing full-prescan operator, preserved failure or specimen."""
from dental_release.paths import expand as _release_expand
import datetime, hashlib, importlib.util, io, json, math, os, sys, zipfile
from pathlib import Path
import numpy as np
import trimesh
sys.dont_write_bytecode = True
R = Path(__file__).resolve().parents[1]
D = R.parent.parent
FC = D / 'results/PROOF_LANE_FULL_CROWN/code'
sys.path[:0] = [str(FC), _release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_GENCAD_V5/deps'), str(D / 'results/LANE_X49_DESIGN_GATE/code')]
import prescan_local_holes as builder

def read(p):
    return json.loads(Path(p).read_text())

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def clean(x):
    if isinstance(x, np.ndarray):
        return clean(x.tolist())
    if isinstance(x, np.generic):
        return clean(x.item())
    if isinstance(x, dict):
        return {k: clean(v) for (k, v) in x.items()}
    if isinstance(x, (tuple, list)):
        return [clean(v) for v in x]
    if isinstance(x, Path):
        return str(x)
    return x

def write(p, x):
    Path(p).write_text(json.dumps(clean(x), indent=2, sort_keys=True) + '\n')

def npz(p, **arrays):
    with zipfile.ZipFile(p, 'w', compression=zipfile.ZIP_DEFLATED) as z:
        for (n, a) in sorted(arrays.items()):
            b = io.BytesIO()
            np.lib.format.write_array(b, np.asarray(a), allow_pickle=False)
            zi = zipfile.ZipInfo(n + '.npy', (1980, 1, 1, 0, 0, 0))
            zi.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(zi, b.getvalue())

def sample(m, n=1024):
    a = m.area_faces
    ids = np.searchsorted(np.cumsum(a), (np.arange(n) + 0.5) * a.sum() / n)
    return m.triangles_center[ids]

def local_cap(m):
    loops = builder.loops_of(m)
    v = m.vertices
    f = m.faces
    verts = v.tolist()
    new = []
    info = []
    for loop in loops:
        pp = v[loop]
        center = pp.mean(0)
        k = len(verts)
        verts.append(center.tolist())
        for (a, b) in zip(loop, np.roll(loop, -1)):
            new.append([b, a, k])
        info.append(dict(center_mm=center.tolist(), n_vertices=len(loop), height_range_mm=[float(pp[:, 2].min()), float(pp[:, 2].max())], measurement_status='UNMEASURED_COMPLETION_HYPOTHESIS'))
    mesh = trimesh.Trimesh(np.array(verts), np.vstack([f, np.array(new, int).reshape(-1, 3)]), process=False)
    trimesh.repair.fix_normals(mesh, multibody=True)
    return (mesh, dict(cap_faces=len(new), loops=info, watertight=bool(mesh.is_watertight), winding=bool(mesh.is_winding_consistent), completion='R4 local centroid fan instead of basal point; no measured boundary added'))

def main(out, round_='R3'):
    p = read(R / ('PREREG_R3_FULL_PRESCAN.json' if round_ == 'R3' else 'PREREG_R4_LOCAL_COMPLETION.json'))
    sel = read(R / 'PATIENT_SELECTION.json')
    source = np.load(out / 'scene.npz')
    tri = source['source_tooth_triangles_mm']
    margin = float(tri[:, :, 2].min() + 0.75)
    if round_ == 'R4':
        builder.cap = local_cap
    result = dict(round=round_, patient_id=sel['patient_id'], fdi=36, source_scene_sha256=sha(out / 'scene.npz'), operator_sha256=sha(FC / 'prescan_local_holes.py'), virtual_margin_mm=margin, frame='Demo1_local_translated_CBCT', status='FAILED', external_referent=p['external_referent'], full_requested_capability=False, physical_pulp_nerve_force_film='UNKNOWN', complete_axial_specimen=False)
    dest = out / ('full_crown' if round_ == 'R3' else 'full_crown_R4')
    dest.mkdir()
    try:
        (m, roles, die, info) = builder.build(tri, margin)
        src = trimesh.Trimesh(tri.reshape(-1, 3), np.arange(len(tri) * 3).reshape(-1, 3), process=True)
        src = src.slice_plane([0, 0, margin], [0, 0, 1], cap=False)
        ext = trimesh.Trimesh(m.vertices, m.faces[roles == 0], process=False)
        (ps, pe) = (sample(src), sample(ext))
        ds = trimesh.proximity.closest_point(ext, ps)[1]
        de = trimesh.proximity.closest_point(src, pe)[1]
        p95 = max(float(np.percentile(ds, 95)), float(np.percentile(de, 95)))
        control_points = ps[::128]
        dd = trimesh.proximity.closest_point(ext, control_points)[1]
        oracle = trimesh.proximity.closest_point_naive(ext, control_points)[1]
        control_error = float(np.max(abs(dd - oracle)))
        npz(dest / 'mesh.npz', vertices=m.vertices, faces=m.faces, face_roles=roles, die_vertices=die.vertices, die_faces=die.faces, source_preservation_distances_mm=ds, reverse_preservation_distances_mm=de)
        m.export(dest / 'crown.stl')
        die.export(dest / 'prep.stl')
        mm = trimesh.load_mesh(dest / 'crown.stl', process=True)
        rounding = float(np.max(abs(m.vertices.astype(np.float32).astype(float) - m.vertices)))
        report = dict(schema='X77-full-prescan-conditional-gate', patient_id=sel['patient_id'], declared_units='mm', material='KATANA_HTML_PLUS', verdict='UNKNOWN' if p95 <= 0.35 else 'FAIL', inputs=dict(crown=dict(sha256=sha(dest / 'crown.stl'), diagnostics=dict(watertight=bool(mm.is_watertight), winding_consistent=bool(mm.is_winding_consistent)))), rules=dict(mesh_health=dict(status='PASS' if mm.is_watertight else 'FAIL', witnesses=[]), scale=dict(status='PASS', physical_scale_calibration='UNKNOWN'), material_wall=dict(status='CONDITIONAL', certificate=info['certified_wall_lower_mm'], scope='parent complete exterior triangles; no physical/scanner uncertainty'), source_preservation=dict(status='PASS' if p95 <= 0.35 else 'FAIL', p95_mm=p95, gate_mm=0.35), insertion=dict(status='CONDITIONAL', ideal=info['ideal_insertion']), cement_physical=dict(status='UNKNOWN'), force=dict(status='UNKNOWN'), pulp=dict(status='UNKNOWN'), nerve=dict(status='UNKNOWN')))
        write(dest / 'DESIGN_GATE.json', report)
        contract = dict(units='mm', common_frame_confirmed=True, regions=dict(exterior=np.flatnonzero(roles == 0).tolist(), intaglio=np.flatnonzero(roles == 1).tolist(), margin=np.flatnonzero(roles == 2).tolist()), region_semantics='Observed/regularized crown exterior, co-designed virtual intaglio and annular closure, not clinical margin', input_sha256=dict(crown=sha(dest / 'crown.stl'), prep=sha(dest / 'prep.stl'), antagonist=sha(out / 'antagonist.stl')))
        write(dest / 'CONTRACT.json', contract)
        from export_bridge import export, check_native
        paths = dict(crown=str(dest / 'crown.stl'), prep=str(dest / 'prep.stl'), antagonist=str(out / 'antagonist.stl'))
        ex = export(paths, report, dest / 'export', dest / 'CONTRACT.json', 1.25)
        if ex['status'] == 'PASS':
            path = dest / 'export/model.3mf'
            with zipfile.ZipFile(path) as zz:
                members = {n: zz.read(n) for n in zz.namelist()}
            with zipfile.ZipFile(path, 'w') as zz:
                for (n, b) in sorted(members.items()):
                    zi = zipfile.ZipInfo(n, (1980, 1, 1, 0, 0, 0))
                    zi.compress_type = zipfile.ZIP_DEFLATED
                    zz.writestr(zi, b)
            side = dest / 'export/model.3mf.json'
            s = read(side)
            s['asset_sha256'] = sha(path)
            write(side, s)
            (dest / 'export/model.3mf.json.sha256').write_text(sha(side) + '\n')
            check = check_native(path, side, sha(side), sha(dest / 'crown.stl'), sha(dest / 'export/DESIGN_GATE.json'), 1.25)
            check.pop('regions_in_imported_face_order', None)
            ex['receipt'].update(asset_sha256=sha(path), sidecar_sha256=sha(side), roundtrip=check)
            write(dest / 'export/EXPORT_RECEIPT.json', ex['receipt'])
        bad = json.loads(json.dumps(report))
        bad['inputs']['crown']['sha256'] = '0' * 64
        try:
            export(paths, bad, dest / 'must_not_export', dest / 'CONTRACT.json', 1.25)
            hash_fault = False
        except ValueError:
            hash_fault = True
        npz(dest / 'observed_fields.npz', source_points=ps, exterior_points=pe, source_distances=ds, exterior_distances=de)
        bad_radius = float(np.linalg.norm(m.triangles[roles == 0] - np.array(info['axis_top']), axis=2).max() + 1.0)
        (bad_wall, _) = builder.wall_bound(m.triangles[roles == 0], np.array(info['axis_top']), bad_radius)
        faults = dict(wrong_distance_rejected=float(np.max(abs(dd + 1 - oracle))) > 1e-09, wrong_export_hash_rejected=hash_fault, injected_radius_mm=bad_radius, injected_wall_bound_mm=bad_wall, wall_below_requirement_rejected=bad_wall < 0.5, control_error_pass=control_error <= 1e-09)
        result.update(status='COMPLETE_AXIAL_RESEARCH_SPECIMEN' if p95 <= 0.35 and mm.is_watertight else 'AXIAL_SPECIMEN_FAILED_SOURCE_OR_SERIALIZATION_GATE', complete_axial_specimen=bool(mm.is_watertight and mm.is_winding_consistent and (info['certified_wall_lower_mm'] >= 0.5)), surface_preservation_gate=p95 <= 0.35, observed_exterior_p95_mm=p95, source_to_ext_p95_mm=float(np.percentile(ds, 95)), ext_to_source_p95_mm=float(np.percentile(de, 95)), source_samples=len(ps), exterior_samples=len(pe), sampled_surface_enclosure='MISSING; p95 is finite area sampling, not uniform geometry error', certificate=info, STL_coordinate_rounding_max_mm=rounding, serialized_closed=bool(mm.is_watertight), vertices=len(m.vertices), faces=len(m.faces), source_control_error_mm=control_error, export=ex, controls=faults, physical_calibration='UNKNOWN', full_CAM='UNKNOWN_HOLDER_FIXTURE_MACHINE', correlation='All Euclidean specimen self-checks share one lower-jaw transform; pulse/canal distances still blocked by missing tissue and calibration')
    except Exception as e:
        result['reason'] = str(e)
    write(out / ('FULL_CROWN.json' if round_ == 'R3' else 'FULL_CROWN_R4.json'), result)
    print(result['status'], result.get('reason', result.get('observed_exterior_p95_mm')), flush=True)
if __name__ == '__main__':
    main(Path(sys.argv[2]), 'R4' if sys.argv[1] == 'full_crown_R4' else 'R3')
