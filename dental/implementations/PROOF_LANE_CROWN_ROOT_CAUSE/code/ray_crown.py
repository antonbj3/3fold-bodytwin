from joint_solid import *

def capsule(c, H, r, z0, n=128, nr=32):
    theta = np.arange(n) * 2 * np.pi / n
    circle = np.c_[np.cos(theta), np.sin(theta)]
    verts = [np.c_[c + r * circle, np.full(n, z0)]]
    for j in range(nr):
        a = j * np.pi / (2 * nr)
        verts.append(np.c_[c + r * np.cos(a) * circle, np.full(n, H + r * np.sin(a))])
    v = np.vstack(verts)
    bot = len(v)
    top = bot + 1
    v = np.r_[v, [[*c, z0], [*c, H + r]]]
    faces = []
    for j in range(nr):
        for i in range(n):
            k = (i + 1) % n
            a = j * n + i
            b = j * n + k
            u = (j + 1) * n + i
            w = (j + 1) * n + k
            faces.extend([[a, b, w], [a, w, u]])
    for i in range(n):
        k = (i + 1) % n
        faces.extend([[bot, k, i], [nr * n + i, nr * n + k, top]])
    m = trimesh.Trimesh(v, np.array(faces), process=False)
    trimesh.repair.fix_normals(m, multibody=True)
    if m.volume < 0:
        m.invert()
    return (m.vertices, m.faces)

def select_ray(v, f, native, key):
    pr = read(R / 'PREREG_C.json')
    S = v[f[native]].copy()
    zmin = float(v[:, 2].min())
    zmax = float(v[:, 2].max())
    rows = []
    for Hrel in pr['search']['top_z_relative_max_mm']:
        H = zmax + Hrel
        zs = np.linspace(zmin, H, max(2, int(np.ceil((H - zmin) / 0.05)) + 1))
        spacing = float(np.diff(zs).max())
        for x in pr['search']['center_offsets_mm']:
            for y in pr['search']['center_offsets_mm']:
                c = np.array([x, y])
                pts = np.c_[np.full(len(zs), x), np.full(len(zs), y), zs]
                distance = fast_nearest(S, pts)[1]
                lb = float(distance.min() - spacing / 2)
                radius = np.floor((lb - 0.65) * 1000000) / 1000000
                wn = float(igl.winding_number(v, f, np.array([[x, y, H]]))[0])
                inside = abs(wn) > 0.5
                feasible = inside and radius >= 0.5
                volume_proxy = np.pi * radius ** 2 * (H - zmin) + 2 * np.pi * radius ** 3 / 3 if feasible else None
                rows.append(dict(center_xy=c.tolist(), axis_top_mm=H, radius_mm=float(radius), axis_clearance_lower_mm=lb, axis_sample_min_mm=float(distance.min()), axis_cover_mm=spacing / 2, inside_top_winding=wn, feasible=bool(feasible), volume_proxy_mm3=volume_proxy, reject_reason=None if feasible else 'AXIS_TOP_OUTSIDE' if not inside else 'RADIUS_LT0P5'))
    dump(R / 'raw' / ('C_SEARCH_' + key + '.json'), rows)
    good = [x for x in rows if x['feasible']]
    if not good:
        raise ValueError('NO_FEASIBLE_RAY_IN100_CANDIDATES')
    return (max(good, key=lambda x: x['volume_proxy_mm3']), rows)

def generate_c(rec):
    st = time.perf_counter()
    (v, f, native, ci) = close(rec)
    (chosen, search) = select_ray(v, f, native, rec['key'])
    c = np.array(chosen['center_xy'])
    H = chosen['axis_top_mm']
    r = chosen['radius_mm']
    z0 = float(v[:, 2].min() - 2)
    (qv, qf) = capsule(c, H, r, z0)
    (kv, kf) = capsule(c, H, r + 0.05, z0)
    support = D / (rec['key'] + '_C_support.mesh')
    cavity = D / (rec['key'] + '_C_cavity.mesh')
    out = D / (rec['key'] + '_C_crown.mesh')
    meshwrite(support, v, f)
    meshwrite(cavity, kv, kf)
    proc = subprocess.run([str(D / 'boolean'), str(support), str(cavity), str(out)], capture_output=True, text=True)
    if proc.returncode:
        raise ValueError('BOOLEAN_FAILURE ' + str(proc.returncode) + ' ' + proc.stdout[:200] + ' ' + proc.stderr[:200])
    boolean = json.loads(proc.stdout)
    (cv, cf) = meshread(out)
    S = v[f[native]].copy()
    A = v[f[~native]].copy()
    source = {tuple(sorted((tuple(x) for x in t))) for t in S}
    kept = np.array([tuple(sorted((tuple(x) for x in t))) in source for t in cv[cf]])
    roles = np.zeros(len(cf), int)
    pts = cv[cf[~kept]].mean(1)
    da = fast_nearest(A, pts)[1]
    dk = fast_nearest(kv[kf], pts)[1]
    roles[~kept] = np.where(dk < da, 1, 2)
    path = D / (rec['key'] + '_C.npz')
    np.savez_compressed(path, vertices=cv, faces=cf, roles=roles, prep_vertices=qv, prep_faces=qf, cavity_vertices=kv, cavity_faces=kf, native_triangles=S, closure_triangles=A, support_vertices=v, support_faces=f)
    result = dict(key=rec['key'], family=rec['family'], status='GENERATED', mesh_path=str(path), mesh_sha256=sha(path), chosen=chosen, search_requested=len(search), search_rejected=sum((not x['feasible'] for x in search)), native_facets=len(S), retained_native_facets=int(kept.sum()), boolean=boolean, closure=ci, wall_analytic_lower_mm=chosen['axis_clearance_lower_mm'] - r - 0.05, seconds=time.perf_counter() - st)
    print('C', rec['key'], 'r', r, 'wallLB', result['wall_analytic_lower_mm'], 'native', int(kept.sum()), len(S), 'seconds', result['seconds'], flush=True)
    return result

def run():
    records = read(B / 'PROOF_LANE_FULL_CROWN_R4/FROZEN_INPUTS_B.json')['records']
    keys = read(R / 'PREREG_A.json')['selection']['keys']
    records = [r for r in records if '--all' in sys.argv or r['key'] in keys]
    rows = []
    st = time.perf_counter()
    for rec in records:
        try:
            out = generate_c(rec)
        except Exception as e:
            out = dict(key=rec['key'], family=rec['family'], status='REJECTED', reason=repr(e))
            print(out, flush=True)
        rows.append(out)
        dump(R / 'raw/C_GENERATION.json', rows)
        dump(R / 'CURRENT_WORK_STATE.json', dict(lane='PROOF_LANE-crown-root-cause', phase='C_CONSTRUCTING', latest_gate=out, next_operation='Freeze geometry then independent wall/gap and reduction checks', review_state='PENDING_INDEPENDENT_REVIEW'))
    p = R / ('FROZEN_PREDICTIONS_C_ALL.json' if '--all' in sys.argv else 'FROZEN_PREDICTIONS_C.json')
    assert not p.exists()
    dump(p, dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), prereg_sha256=sha(R / 'PREREG_C.json'), code_sha256=sha(Path(__file__)), rows=rows, seconds=time.perf_counter() - st, physical_predictions=None, scope='Geometric virtual construction with native input; no clinical efficacy prediction'))
    p.with_suffix('.sha256').write_text(sha(p) + '  ' + p.name + '\n')
if __name__ == '__main__':
    run()
