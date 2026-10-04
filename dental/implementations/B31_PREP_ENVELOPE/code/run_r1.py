from geometry import *
from carriers import *
from shapely.prepared import prep
from shapely.geometry import box
from column_sign import column_sign

def admitted_field(t, closed, base):
    h = float(H)
    v = closed.vertices
    axes = [np.arange(np.floor(v[:, i].min() / h) * h + h / 2, np.ceil(v[:, i].max() / h) * h, h) for i in range(3)]
    axes[2] = np.arange(base + h / 2, np.ceil(t[:, :, 2].max() / h) * h, h)
    origin = np.array([a[0] for a in axes])
    pts = np.stack(np.meshgrid(*axes, indexing='ij'), -1).reshape(-1, 3)
    (inside, parity_control) = column_sign(closed, axes)
    dist = np.zeros(len(pts))
    dist[inside] = distances(source_mesh(t), pts[inside])
    boundary_nonbase = closed.triangles[~np.all(abs(closed.triangles[:, :, 2] - base) < 1e-09, axis=1)]
    dclosed = np.zeros(len(pts))
    dclosed[inside] = distances(source_mesh(boundary_nonbase), pts[inside])
    rad = np.sqrt(3) * h / 2
    admitted = inside & (dist >= 0.6 + rad) & (dclosed >= rad + 1e-10)
    T = admitted.reshape(tuple((len(a) for a in axes)))
    good = np.flatnonzero(T.any(axis=(0, 1)))
    if not len(good):
        raise ValueError('NO_FULL_ERODED_BOXES')
    start = int(good[0])
    T = T[:, :, start:]
    origin[2] += h * start
    carrier_base = origin[2] - h / 2
    return (T, origin, carrier_base, {'source_grid_points': len(pts), 'inside_samples': int(inside.sum()), 'admitted_boxes': int(T.sum()), 'wall_guard_mm': 0.6, 'box_radius_mm': rad, 'source_closure_sign': 'vertical-line parity; source closure self intersection UNKNOWN', 'parity_control': parity_control, 'rigorous_roundoff_enclosure': 'MISSING'})

def solve(r, closure):
    tic = time.perf_counter()
    d = load_np(closure['path'])
    t = d['original']
    closed = trimesh.Trimesh(d['vertices'], d['faces'], process=False)
    (T, origin, base, meta) = admitted_field(t, closed, float(d['base']))
    (F, C, lower, upper) = dictionary(T)
    check = explicit_check(T, F, C)
    rr = {'key': r['key'], 'case_key': r['case_key'], 'family': r['family'], 'resolution': 'PER_TOOTH', 'field': meta, 'carrier_count': int(F.sum()), 'undominated_count': len(C), 'upper_outside_admitted_boxes': int((upper & ~T).sum()), 'independent_checks': len(check), 'independent_failures': sum((bool(x['outside_boxes']) for x in check)), 'retained_volume_interval_mm3': [float(lower.sum()) * float(H) ** 3, float(upper.sum()) * float(H) ** 3], 'cone_canonical': 'EXACT_RATIONAL_PRIMITIVE_PROPERTY; source admission arithmetic not enclosed', 'direction': [0, 0, 1], 'full_source_certification': 'UNKNOWN_SELF_INTERSECTION_AND_ROUNDOFF'}
    path = DATA / 'r1' / (r['key'] + '_field.npz')
    path.parent.mkdir(exist_ok=True)
    np.savez_compressed(path, T=T, F=F, carriers=C, origin=origin, base=base, union_lower=lower, union_upper=upper)
    rr.update(field_path=path, field_sha256=sha(path))
    dump(ROOT / 'raw/r1' / (r['key'] + '_checks.json'), check)
    if not len(C):
        rr.update(status='UNKNOWN_NO_CARRIER', missing='No feasible dictionary carrier; not a global impossibility proof')
        return rr
    (inner, coords, info) = union_surface(C, origin, base)
    rr['surface'] = info
    ii = ~np.all(abs(inner.triangles[:, :, 2] - base) < 2e-06, axis=1)
    n = inner.face_normals[ii]
    values = n[:, 2] - float(S) * np.abs(n[:, :2]).sum(1)
    rr['export_insertion'] = {'axial_min_normal_z': float(n[:, 2].min()), 'square_cone_min_margin': float(values.min()), 'square_cone_failing_facets': int((values < -1e-08).sum()), 'axial_failing_facets': int((n[:, 2] < -1e-08).sum()), 'resolution': 'PER_SURFACE_REGION', 'fullpath_certificate': 'UNKNOWN_ROUNDOFF_AND_OPENING_BRIDGE'}
    try:
        (shell, roles, opened) = crown_shell(t, inner, base)
        wm = wall_bound(t, opened.triangles)
        q = DATA / 'r1' / (r['key'] + '_candidate.npz')
        np.savez_compressed(q, vertices=shell.vertices, faces=shell.faces, roles=roles, original=t, carrier_coords=coords, base=base)
        stl = DATA / 'r1' / (r['key'] + '_candidate_UNQUALIFIED.stl')
        shell.export(stl)
        ext = shell.vertices[shell.faces[roles == 0]]
        sort = lambda a: np.sort(a, axis=1)
        rr.update(status='CANDIDATE_WITH_UNKNOWN', wall=wm, closed_shell=bool(shell.is_watertight and shell.is_winding_consistent and (shell.volume > 0)), shell_components=len(shell.split(only_watertight=False)), exterior_coordinate_error_mm=float(np.max(abs(sort(ext) - sort(t)))), margin_coordinate_error_mm=0.0, shell_self_intersection='UNKNOWN_NOT_YET_CHECKED', candidate_path=q, candidate_sha256=sha(q), stl_path=stl, stl_sha256=sha(stl), complete_digital_status='UNKNOWN', missing=['source closure self-intersection', 'whole shell/sweep collision certificate', 'rigorous floating/source enclosure', 'actual preparation/finish line and product qualification'])
    except Exception as e:
        rr.update(status='UNKNOWN_SURFACE_RECONSTRUCTION', reason=str(e), complete_digital_status='UNKNOWN')
    rr['seconds'] = time.perf_counter() - tic
    return rr

def main():
    rows = []
    tic = time.perf_counter()
    cohort = read(ROOT / 'FROZEN_COHORT.json')['records']
    closures = {r['key']: r for r in read(ROOT / 'raw/SOURCE_INSPECTION.json')}
    for r in cohort:
        c = closures[r['key']]
        if c['status'] != 'VIRTUAL_SOURCE_CLOSED':
            rr = {'key': r['key'], 'family': r['family'], 'status': 'UNKNOWN_SOURCE_CLOSURE', 'reason': c['reason'], 'complete_digital_status': 'UNKNOWN', 'resolution': 'PER_TOOTH'}
        else:
            try:
                rr = solve(r, c)
            except Exception as e:
                rr = {'key': r['key'], 'family': r['family'], 'status': 'UNKNOWN_CONSTRUCTION', 'reason': repr(e), 'complete_digital_status': 'UNKNOWN', 'resolution': 'PER_TOOTH'}
        rows.append(rr)
        dump(ROOT / 'raw/R1_RESULTS.json', rows)
        state('R1_RUNNING', {'completed': len(rows), 'requested': 18}, 'Finish same 18 then change the binding basal closure/surface representation')
        print(rr['key'], rr['status'], rr.get('reason', ''), rr.get('carrier_count'), rr.get('export_insertion'), flush=True)
    dump(ROOT / 'raw/R1_COST.json', {'seconds': time.perf_counter() - tic, **usage()})
    state('R1_DECIDED', {'complete': 0, 'requested': 18, 'unknown': 18}, 'Preserve R1; change basal closure to a spatial triangle fan, retain original native rim and new prereg')
if __name__ == '__main__':
    main()
