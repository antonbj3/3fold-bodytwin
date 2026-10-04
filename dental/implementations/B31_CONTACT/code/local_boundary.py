from common import *
from geometry import *
from optimize import solve
from bounds import box
from shapely import polygons, STRtree
import time, resource

def overlap_bary(q, upper, lower):
    xy = lower[:, :2]
    det = (xy[1, 0] - xy[0, 0]) * (xy[2, 1] - xy[0, 1]) - (xy[1, 1] - xy[0, 1]) * (xy[2, 0] - xy[0, 0])
    if det < 0:
        lower = lower[[0, 2, 1]]
        xy = lower[:, :2]
    p = np.eye(3)
    for (a, b) in zip(xy, np.roll(xy, -1, axis=0)):
        u = upper[:, :2] - a
        edge = b - a
        values = edge[0] * u[:, 1] - edge[1] * u[:, 0]
        p = q.Q.clip(p, values)
        if len(p) < 3:
            return (None, None)
    points = p @ upper[:, :2]
    (ar, _) = q.measure(p, upper[:, :2])
    if ar <= 1e-14:
        return (None, None)
    mat = np.c_[lower[:, :2], np.ones(3)]
    coefficients = np.linalg.solve(mat, lower[:, 2])
    zl = np.c_[points, np.ones(len(points))] @ coefficients
    p = np.maximum(p, 0)
    p /= p.sum(1)[:, None]
    return (p, p @ upper[:, 2] - zl)

def build(q, t, r):
    start = time.perf_counter()
    v = t['vertices']
    f = t['faces']
    tri = v[f]
    nz = np.cross(tri[:, 1, :2] - tri[:, 0, :2], tri[:, 2, :2] - tri[:, 0, :2])
    upids = np.flatnonzero((t['face_roles'] == 0) & (nz > 1e-12) & np.any(t['taper'][f] > 0, axis=1))
    downids = np.flatnonzero(nz < -1e-12)
    tree = STRtree(polygons(tri[downids, :, :2]))
    rows = []
    cols = []
    vals = []
    rhs = []
    pairs = 0
    vertices = 0
    crossings = []
    pruned = 0
    pruned_lb = np.inf
    fault = None
    for chunk in np.array_split(upids, max(1, int(np.ceil(len(upids) / 512)))):
        candidates = tree.query(polygons(tri[chunk, :, :2]))
        for (ui, di) in candidates.T:
            uf = int(chunk[ui])
            df = int(downids[di])
            (p, separation) = overlap_bary(q, tri[uf], tri[df])
            if p is None:
                continue
            pairs += 1
            vertices += len(p)
            if separation.min() < -1e-09 and separation.max() > 1e-09:
                crossings.append(dict(upper_face=uf, lower_face=df, separation_min_mm=float(separation.min()), separation_max_mm=float(separation.max())))
                continue
            if separation.max() <= 1e-09 and separation.min() < -1e-09:
                continue
            for (a, s) in zip(p, separation):
                s = float(s)
                maxmotion = r['relief_cap_mm'] * float(a @ t['taper'][f[uf]])
                if fault is None and a.max() > 0.1 and (s > 1e-06) and np.any(t['taper'][f[uf]] > 0):
                    fault = dict(upper_face=uf, lower_face=df, upper_vertex_ids=f[uf], barycentric=a, separation_mm=s)
                if s >= maxmotion:
                    pruned += 1
                    pruned_lb = min(pruned_lb, s - maxmotion)
                    continue
                rid = len(rhs)
                rows.extend([rid] * 3)
                cols.extend(f[uf])
                vals.extend(a)
                rhs.append(max(0.0, s))
    A = csr_matrix((vals, (rows, cols)), shape=(len(rhs), len(v)))
    b = np.array(rhs)
    stats = dict(upward_potentially_moved_faces=len(upids), fixed_downward_faces=len(downids), positive_area_overlap_pairs=pairs, overlap_polygon_vertices_tested=vertices, cap_implied_constraints=pruned, cap_implied_minimum_slack_lower_mm=None if not np.isfinite(pruned_lb) else float(pruned_lb), explicit_local_constraints=len(rhs), original_mixed_order_overlap_pairs=len(crossings), original_crossing_witnesses=crossings[:5], seconds=time.perf_counter() - start, full_original_self_intersection_proof='UNKNOWN', rigorous_machine_enclosure='MISSING')
    return (A, b, stats, fault)

def run():
    check_lock()
    start = time.perf_counter()
    (q, reuse) = parents()
    pr = read(ROOT / 'PREREG_R3.json')
    old = read(ROOT / 'raw/ROUND1.json')['rows']
    rows = []
    predictions = []
    cache = {}
    for (r, rr) in zip(pr['cohort'], old):
        t = load_case(r)
        (A, b, stats, fault) = build(q, t, r)
        cp = DATA / 'R3' / (r['key'] + '_constraints.npz')
        cp.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(cp, A_data=A.data, A_indices=A.indices, A_indptr=A.indptr, A_shape=A.shape, b=b)
        row = dict(key=r['key'], case_key=r['case_key'], local_polytope=stats, constraints_path=str(cp), constraints_sha256=sha(cp), outputs={})
        pred = dict(key=r['key'], outputs={})
        for (name, reg) in [('regional', True), ('ordinary_mesh', False)]:
            (d, opt) = solve(t, r, reg, extra=(A, b))
            (v, g) = wall_and_removal(t, d, r)
            predg = gap(t, d)
            predcontact = q.compare(t['xy'], t['grid_faces'], predg, t['reference_gap'])
            p = DATA / 'R3' / r['key'] / (name + '.npz')
            p.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(p, vertices=v, faces=t['faces'], face_roles=t['face_roles'], relief_mm=d, predicted_gap=predg)
            export_stl(ROOT / 'exports' / r['key'] / ('R3_' + name + '.stl'), v, t['faces'])
            with np.load(rr['outputs'][name]['path']) as a:
                identical = bool(np.array_equal(v, a['vertices']))
                difference = float(abs(v - a['vertices']).max())
            slack = b - A @ d
            minslack = float(slack.min()) if len(slack) else None
            pruned_lower = stats['cap_implied_minimum_slack_lower_mm']
            conditional_order_pass = stats['original_mixed_order_overlap_pairs'] == 0 and (minslack is None or minslack >= -1e-08) and (g['protected_identity_error_mm'] == 0.0)
            g['local_vertical_interval_order_pass'] = conditional_order_pass
            g['local_interval_slack_lower_mm'] = minslack
            g['cap_implied_pruned_slack_lower_mm'] = pruned_lower
            g['removal_status_R3'] = 'CONDITIONAL_ON_ORIGINAL_VALID_SOLID' if conditional_order_pass else 'UNKNOWN_OR_ORIGINAL_ORDER_FAIL'
            row['outputs'][name] = dict(path=str(p), sha256=sha(p), optimization=opt, wall_and_removal=g, R1_output_bit_identical=identical, R1_max_vertex_difference_mm=difference, nominal_contact_prediction=predcontact)
            pred['outputs'][name] = dict(artifact_sha256=sha(p), nominal_symdiff_mm2=predcontact['symdiff_mm2'], wall_lower_mm=g['wall_lower_mm'], area_integral_mm3=g['area_integral_abs_displacement_mm3'], force_N=None, physical_adjustment_time_min=None)
        if fault is not None:
            ids = np.asarray(fault['upper_vertex_ids'])
            a = np.asarray(fault['barycentric'])
            eligible = t['taper'][ids] > 0
            weights = np.where(eligible, a, -1)
            k = int(weights.argmax())
            injected_relief = (fault['separation_mm'] + 0.01) / a[k]
            violation = a[k] * injected_relief - fault['separation_mm']
            fault.update(injected_vertex=int(ids[k]), injected_relief_mm=float(injected_relief), local_violation_mm=float(violation), rejected=violation > 1e-08, scope='Local-boundary fault may also violate inherited wall cap; tested independently of that cap')
        row['downward_crossing_fault'] = fault
        rows.append(row)
        predictions.append(pred)
        cache[r['key']] = (t, A, b)
        dump(ROOT / 'raw/R3_CONSTRUCTION_PROGRESS.json', dict(rows=rows))
        print('local constructed', r['key'], 'vertices', stats['overlap_polygon_vertices_tested'], 'constraints', len(b), 'old crossings', stats['original_mixed_order_overlap_pairs'], flush=True)
    fp = ROOT / 'FROZEN_PREDICTIONS_R3.json'
    if not fp.exists():
        freeze(fp, dict(prereg_sha256=sha(ROOT / 'PREREG_R3.json'), predictions=predictions, physical_measurement='NOT_RUN', source_used_for_design=True, next='Actual mesh direct-height validation; box reuse only on bit-identical whole vertices'))
    else:
        for (a, b) in zip(read(fp)['predictions'], predictions):
            for n in a['outputs']:
                if a['outputs'][n]['artifact_sha256'] != b['outputs'][n]['artifact_sha256']:
                    raise ValueError('R3 frozen output drift')
    for (r, row, rr) in zip(pr['cohort'], rows, old):
        (t, A, b) = cache[r['key']]
        for (name, x) in row['outputs'].items():
            with np.load(x['path']) as a:
                v = a['vertices']
                d = a['relief_mm']
            actual = t['ceiling'] - reuse.VF.height(v[t['faces'][t['face_roles'] == 0]], t['xy'])
            contact = q.compare(t['xy'], t['grid_faces'], actual, t['reference_gap'])
            x['contact'] = contact
            x['prediction_symdiff_error_mm2'] = abs(contact['symdiff_mm2'] - x['nominal_contact_prediction']['symdiff_mm2'])
            if x['R1_output_bit_identical']:
                x['boxes'] = rr['outputs'][name]['boxes']
                x['box_binding'] = 'Bit-identical complete vertices/faces and unchanged source gap/grid; inherited R1 raw/ROUND1.json'
            else:
                x['boxes'] = {str(b): box(q, t, actual, t['reference_gap'], b, pr['parameters']) for b in [0.01, 0.05]}
                x['box_binding'] = 'R3 direct continuous interval recomputation'
            x['local_bound_injected_plus001_rejected'] = row['downward_crossing_fault'] is not None and row['downward_crossing_fault']['rejected']
            x['physical_eligibility'] = 'UNKNOWN_ORIGINAL_SELF_INTERSECTION_ROUNDOFF_PROCESS_AND_MEASUREMENT'
        print('local validated', r['key'], 'nominal', row['outputs']['regional']['contact']['symdiff_mm2'], flush=True)
    summary = dict(claim_type='capability', retained=10, case_clusters=4, total_overlap_vertices_tested=sum((r['local_polytope']['overlap_polygon_vertices_tested'] for r in rows)), original_mixed_order_cases=sum((r['local_polytope']['original_mixed_order_overlap_pairs'] > 0 for r in rows)), conditional_local_order_regional_passes=sum((r['outputs']['regional']['wall_and_removal']['local_vertical_interval_order_pass'] for r in rows)), R1_regional_outputs_bit_identical=sum((r['outputs']['regional']['R1_output_bit_identical'] for r in rows)), physical_passes=0, inherited_digital_gate='FAIL_PRESERVED; R2 all-family counterwitness still applies', seconds=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    dump(ROOT / 'raw/ROUND3.json', dict(summary=summary, rows=rows))
    state('R3_DECIDED', summary, 'Seal one-command demo and isolated replay; missing observation: valid source solid/loaded pose plus fabricated before/after scan and active time')
    print(json.dumps(summary), flush=True)
if __name__ == '__main__':
    run()
