from run_r1 import admitted_field
from geometry import *
from carriers import *
from exact_union import *
from controls import triangle_lp

def fan_closure(t, margin, base):
    m = source_mesh(t)
    ids = boundary(m)
    rim = m.vertices[ids]
    apex = np.r_[rim[:, :2].mean(0), base]
    cap = np.array([[rim[i], rim[(i + 1) % len(rim)], apex] for i in range(len(rim))])
    closed = source_mesh(np.r_[t, cap])
    trimesh.repair.fix_normals(closed)
    if closed.volume < 0:
        closed.invert()
    if not (closed.is_watertight and closed.is_winding_consistent):
        raise ValueError('SPATIAL_FAN_TOPOLOGY_FAILURE')
    return (closed, cap)

def separating_axis(A, B):
    ea = np.roll(A, -1, axis=0) - A
    eb = np.roll(B, -1, axis=0) - B
    na = np.cross(ea[0], ea[1])
    nb = np.cross(eb[0], eb[1])
    axes = [na, nb] + [np.cross(a, b) for a in ea for b in eb] + [np.cross(na, e) for e in ea] + [np.cross(nb, e) for e in eb]
    for ax in axes:
        mag = np.linalg.norm(ax)
        if mag < 1e-15:
            continue
        x = A @ (ax / mag)
        y = B @ (ax / mag)
        gap = max(float(x.min() - y.max()), float(y.min() - x.max()))
        if gap > 1e-09:
            return {'separated': True, 'gap_mm': gap, 'roundoff_enclosure': 'MISSING'}
    return {'separated': False}

def check_added_surfaces(original, added, max_lp=256):
    m = source_mesh(original)
    tree = m.triangles_tree
    checked = 0
    separated = 0
    shared = 0
    unknown = 0
    hits = []
    LP = []
    for (k, A) in enumerate(added):
        bb = np.r_[A.min(0) - 1e-09, A.max(0) + 1e-09]
        for j in tree.intersection(bb):
            B = original[j]
            checked += 1
            if any((np.array_equal(a, b) for a in A for b in B)):
                shared += 1
                continue
            sat = separating_axis(A, B)
            if sat['separated']:
                separated += 1
                continue
            if len(LP) >= max_lp:
                unknown += 1
                continue
            w = triangle_lp(A, B)
            w.update(added_facet=k, source_facet=int(j))
            LP.append(w)
            if w['status'] == 'NUMERIC_INTERSECTION_WITNESS':
                hits.append(w)
            else:
                unknown += 1
    return {'AABB_pairs': checked, 'floating_SAT_separations': separated, 'shared_vertex_pairs_excluded': shared, 'LP_queries': len(LP), 'numeric_intersection_witnesses': len(hits), 'unresolved_pairs': unknown, 'witnesses': hits[:32], 'LP_raw': LP, 'whole_surface_status': 'COUNTERWITNESS_FLOATING' if hits else 'UNKNOWN_SHARED_ADJACENCY_AND_ROUNDOFF', 'scope': 'added faces vs original only; original self-intersections and added/added still unresolved'}

def solve2(r, field_cache=None):
    tic = time.perf_counter()
    t = load_np(r['private_path'])['target']
    p = load_np(r['public_path'])
    base0 = np.floor((p['margin_curve'][:, 2].min() - 0.2) / 0.2) * 0.2
    (closed, cap) = fan_closure(t, p['margin_curve'], base0)
    if field_cache is None:
        (T, origin, base, meta) = admitted_field(t, closed, base0)
        (F, C, lower, upper) = dictionary(T)
    else:
        assert sha(field_cache['field_path']) == field_cache['field_sha256']
        fd = load_np(field_cache['field_path'])
        assert np.array_equal(fd['source'], t) and np.array_equal(fd['source_cap'], cap)
        (T, origin, base, meta) = (fd['T'], fd['origin'], float(fd['base']), field_cache['field'])
        (F, C, lower, upper) = dictionary(T)
        assert np.array_equal(F, fd['F']) and np.array_equal(C, fd['carriers'])
    check = explicit_check(T, F, C)
    out = DATA / 'r2'
    out.mkdir(exist_ok=True)
    fp = out / (r['key'] + '_field.npz')
    np.savez_compressed(fp, T=T, F=F, carriers=C, origin=origin, base=base, source=t, source_cap=cap)
    dump(ROOT / 'raw/r2' / (r['key'] + '_FIELD_CHECKPOINT.json'), {'key': r['key'], 'field_path': fp, 'field_sha256': sha(fp), 'field': meta, 'carrier_count': int(F.sum()), 'undominated_count': len(C), 'upper_outside_admitted_boxes': int((upper & ~T).sum()), 'independent_checks': len(check), 'independent_failures': sum((bool(q['outside_boxes']) for q in check)), 'source_private_sha256': r['private_sha256']})
    coords = origin + float(H) * C
    (pair, selection) = choose_pair(coords, base)
    coords = coords[pair]
    r0 = H / 2 - S * H / 2
    P = [primitive(tuple((Q(float(x)) for x in c[:2])), Q(float(c[2])), Q(float(base)), r0) for c in coords]
    (inner, cert, rawpts, face_roles) = mesh_exact(P)
    PP = [primitive(tuple((Q(float(x)) for x in c[:2])), Q(float(c[2])) - Q(1, 20), Q(float(base)), r0 + S * Q(1, 20) - Q(51, 1000)) for c in coords]
    (prep, pcert, _, prep_roles) = mesh_exact(PP)
    rr = {'key': r['key'], 'case_key': r['case_key'], 'family': r['family'], 'resolution': 'PER_TOOTH', 'field': meta, 'native_margin_z_span_mm': np.ptp(p['margin_curve'][:, 2]), 'carrier_count': int(F.sum()), 'undominated_count': len(C), 'upper_outside_admitted_boxes': int((upper & ~T).sum()), 'independent_checks': len(check), 'independent_failures': sum((bool(x['outside_boxes']) for x in check)), 'maximum_dictionary_retained_interval_mm3': [float(lower.sum()) * float(H) ** 3, float(upper.sum()) * float(H) ** 3], 'pair_selection': selection, 'canonical_intaglio': cert, 'canonical_preparation': pcert, 'cement_min_clearance_model_mm': 0.05, 'cement_side_support_shift_mm': 0.051, 'original_exterior_coordinate_error_mm': 0.0, 'native_margin_coordinate_error_mm': 0.0, 'actual_source_containment': 'UNKNOWN_SPATIAL_FAN_AND_SOURCE_TRUTH', 'complete_digital_status': 'UNKNOWN', 'source_closure': 'CONSTITUTIVE_SPATIAL_FAN', 'direction': [0, 0, 1]}
    out = DATA / 'r2'
    out.mkdir(exist_ok=True)
    fp = out / (r['key'] + '_field.npz')
    np.savez_compressed(fp, T=T, F=F, carriers=C, origin=origin, base=base, pair_ids=pair, source=t, source_cap=cap)
    rr.update(field_path=fp, field_sha256=sha(fp))
    dump(ROOT / 'raw/r2' / (r['key'] + '_carrier_checks.json'), check)
    if not cert['watertight'] or cert['exact_cone_failures'] or cert['exact_volume_mm3'] != selection['selected_pair_exact_volume_mm3']:
        raise ValueError('EXACT_UNION_CERTIFICATE_OR_VOLUME_FAILED')
    try:
        (shell, roles, opened) = crown_shell(t, inner, base)
        wm = wall_bound(t, opened.triangles)
        stl = out / (r['key'] + '_crown_UNQUALIFIED.stl')
        shell.export(stl)
        ps = out / (r['key'] + '_prep_UNQUALIFIED.stl')
        prep.export(ps)
        q = out / (r['key'] + '_candidate.npz')
        np.savez_compressed(q, vertices=shell.vertices, faces=shell.faces, roles=roles, intaglio_vertices=inner.vertices, intaglio_faces=inner.faces, prep_vertices=prep.vertices, prep_faces=prep.faces, carrier_coords=coords, base=base, original=t)
        added = shell.triangles[roles == 2]
        cc = check_added_surfaces(t, added)
        dump(ROOT / 'raw/r2' / (r['key'] + '_surface_control.json'), cc)
        rr.update(status='CANDIDATE_WITH_UNKNOWN', closed_shell=bool(shell.is_watertight and shell.is_winding_consistent and (shell.volume > 0)), shell_components=len(shell.split(only_watertight=False)), wall=wm, candidate_path=q, candidate_sha256=sha(q), stl_path=stl, stl_sha256=sha(stl), prep_stl_path=ps, prep_stl_sha256=sha(ps), added_source_control={k: v for (k, v) in cc.items() if k not in ['witnesses', 'LP_raw']}, missing=['source fan self intersection and true inside', 'whole shoulder/preparation swept-path certificate', 'source/original self intersection and float arithmetic enclosure', 'actual margin/preparation/product qualification'])
    except Exception as e:
        rr.update(status='UNKNOWN_SHELL_RECONSTRUCTION', reason=str(e))
    source = next((s for s in read(ROOT / 'FROZEN_COHORT.json')['sources'] if s['case_key'] == r['case_key']))
    iq = out / (r['key'] + '_design_spatial.json')
    dump(iq, {'center_apex_mm': coords, 'base_mm': base, 'r0_mm': str(r0), 'slope': str(S), 'exact_union_certificate': cert, 'exact_preparation_certificate': pcert, 'integer_rational_vertices': [[str(x) for x in p] for p in rawpts], 'surface_faces': inner.faces, 'basal_face_roles': face_roles, 'source_private_sha256': r['private_sha256'], 'frame': 'R4 tooth local mm', 'coordinate_contract': {'case_key': r['case_key'], 'tooth_fdi': r['fdi'], 'units': 'mm', 'origin_world_rotated_mm': p['origin_world'], 'arch_rotation_rows': source['frame'], 'to_source_obj_equation': 'row point_source_OBJ = (point_tooth_local + origin_world_rotated) @ arch_rotation_rows', 'pose': 'unchanged source-derived orientation; verified bite UNKNOWN', 'region': 'original exterior, virtual cervical shoulder, two-core intaglio and preparation', 'time': 'SIMULTANEOUS digital design snapshot', 'resolution': 'PER_POINT'}, 'whole_crown': 'UNKNOWN'})
    rr.update(carrier_design_path=iq, carrier_design_sha256=sha(iq))
    rr['seconds'] = time.perf_counter() - tic
    return rr

def main():
    tic = time.perf_counter()
    rows = []
    for r in read(ROOT / 'FROZEN_COHORT.json')['records']:
        try:
            assert sha(r['private_path']) == r['private_sha256']
            assert sha(r['public_path']) == r['public_sha256']
            rr = solve2(r)
        except Exception as e:
            rr = {'key': r['key'], 'family': r['family'], 'status': 'UNKNOWN_CONSTRUCTION', 'reason': repr(e), 'complete_digital_status': 'UNKNOWN', 'resolution': 'PER_TOOTH'}
        rows.append(rr)
        dump(ROOT / 'raw/R2_RESULTS.json', rows)
        state('R2_RUNNING', {'completed': len(rows), 'requested': 18}, 'Finish all18 source-bound exact carrier constructions and whole-surface independent checks')
        print(rr['key'], rr['status'], rr.get('reason'), rr.get('canonical_intaglio'), rr.get('wall', {}).get('geometric_lower_mm'), flush=True)
    dump(ROOT / 'raw/R2_COST.json', {'seconds': time.perf_counter() - tic, **usage()})
    state('R2_DECIDED', {'complete': 0, 'requested': 18, 'unknown': 18}, 'Preserve two-core result; next construction must certify shoulder or remove unmeasured marginal volume closure')
if __name__ == '__main__':
    main()
