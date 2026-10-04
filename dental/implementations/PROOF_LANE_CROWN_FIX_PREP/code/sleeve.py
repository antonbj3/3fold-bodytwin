from core import *
from generate import construct

def surfaces(rec):
    info = read(R / 'raw' / (rec['key'] + '_R2_SUPPORTS.json'))
    T = (*meshread(info['native']['path']), info['native'])
    E = (*meshread(info['outer']['path']), info['outer'])
    S = load(rec['private_path'])['target']
    x = load(rec['outer']['mesh_path'])
    (ev, ef) = (x['vertices'], x['faces'])
    if 'roles' in x:
        (ev, ef) = compact(ev, ef[x['roles'] == 0])
    (oldv, oldf) = meshread(D / (rec['key'] + '_outer_closed.mesh'))
    original = {tuple(sorted((tuple(y) for y in t))) for t in ev[ef]}
    mask = np.array([tuple(sorted((tuple(y) for y in t))) in original for t in oldv[oldf]])
    artificial = oldv[oldf[~mask]]
    (av, af) = triangles_mesh(artificial)
    pts = E[0][E[1]].mean(1)
    de = Distance(ev, ef).query(pts)[0]
    dc = Distance(av, af).query(pts)[0]
    ext = de <= dc
    outside = E[0][E[1][ext]]
    return (T, E, S, outside, ext)

def sep(path, queries, key, surface=True):
    p = D / (key + '_queries.txt')
    with p.open('w') as h:
        h.write(str(len(queries)) + '\n')
        np.savetxt(h, queries, fmt='%.17g')
    cmd = [str(D / 'exact_surface'), str(path), str(p)] + (['surface'] if surface else [])
    z = subprocess.run(cmd, capture_output=True, text=True, timeout=90)
    if z.returncode:
        raise ValueError('SURFACE_EXACT_EXIT ' + str(z.returncode) + ' ' + z.stderr[:100])
    out = json.loads(z.stdout)
    ok = all((x['pass'] and (surface or x['inside']) for x in out['rows']))
    out.update(all_pass=ok, inside_checked=not surface, support_path=str(path), support_sha256=sha(path), query_path=str(p), query_sha256=sha(p))
    dest = R / 'raw' / (key + '_EXACT.json')
    dump(dest, out)
    return dict(path=str(dest), sha256=sha(dest), all_pass=ok, inside_checked=not surface, balls=len(queries), scope='Exact surface separation ONLY, not containment' if surface else 'Exact surface separation and inside for complete closed support')

def search3(rec, mat, T, E, S, X):
    pr = read(R / 'PREREG_R3.json')
    s = pr['search']
    (sv, sf) = triangles_mesh(S)
    (xv, xf) = triangles_mesh(X)
    sd = Distance(sv, sf)
    xd = Distance(xv, xf)
    td = Distance(T[0], T[1])
    ed = Distance(E[0], E[1])
    best = None
    rows = []
    for (mi, M) in enumerate(frame_candidates(rec['key'])):
        st = S @ M.T
        xt = X @ M.T
        low = min(st[:, :, 2].min(), xt[:, :, 2].min()) - 2.0
        hi = min(st[:, :, 2].max(), xt[:, :, 2].max())
        center = (st.reshape(-1, 3).min(0) + st.reshape(-1, 3).max(0)) / 2
        for extra in s['top_reduction_extra_mm']:
            H = hi - mat['top_reduction_mm'] - extra
            if H - low < 0.75:
                continue
            zs = np.linspace(low, H, max(2, int(np.ceil((H - low) / 0.1)) + 1))
            mid = (zs[:-1] + zs[1:]) / 2
            cover = np.diff(zs) / 2
            off = offset(zs[:-1], H)
            for dx in s['center_offsets_mm']:
                for dy in s['center_offsets_mm']:
                    c = center[:2] + [dx, dy]
                    pts = np.c_[np.tile(c, (len(mid), 1)), mid] @ M
                    ar = float(min(np.min(sd.query(pts)[0] - cover - off - mat['source_reduction_mm'] - 0.025), np.min(xd.query(pts)[0] - cover - off - mat['wall_mm'] - G - 0.025)))
                    a = np.floor(ar * 1000000.0) / 1000000.0
                    inside = bool(td.inside(np.array([[*c, H]]) @ M)[0] and ed.inside(np.array([[*c, H]]) @ M)[0]) if a >= 0.15 else False
                    good = a >= 0.15 and inside
                    proxy = float(np.pi * np.trapezoid((a + offset(zs, H)) ** 2, zs)) if good else None
                    r = dict(axis_id=mi, M=M, c=c, b=low, H=H, flat_radius_mm=a, feasible=good, volume_proxy_mm3=proxy, reject_reason=None if good else 'RADIUS_TOO_SMALL' if a < 0.15 else 'TOP_OUTSIDE')
                    rows.append(r)
                    if good and (best is None or proxy > best['volume_proxy_mm3']):
                        best = r
    path = D / (rec['key'] + '_' + mat['id'] + '_R3_SEARCH.json')
    dump(path, rows)
    return (best, dict(path=str(path), sha256=sha(path), requested=len(rows), rejected=sum((not x['feasible'] for x in rows)), reasons={k: sum((x['reject_reason'] == k for x in rows)) for k in ['RADIUS_TOO_SMALL', 'TOP_OUTSIDE']}))

def build(rec, mat, T, E, S, X, emask, ch):
    key = rec['key'] + '_' + mat['id'] + '_R3'
    M = np.array(ch['M'])
    c = np.array(ch['c'])
    a = ch['flat_radius_mm']
    b = ch['b']
    H = ch['H']
    (qv, qf, zq, _) = body(c, a, b, H, M)
    (kv, kf, zk, _) = body(c, a, b, H + G, M, rho=RHO + G)
    (sv, sf) = triangles_mesh(S)
    xp = D / (key + '_exterior.mesh')
    sp = D / (key + '_native.mesh')
    (xv, xf) = triangles_mesh(X)
    meshwrite(xp, xv, xf)
    meshwrite(sp, sv, sf)
    qs = sep(sp, query_ball_rows(qv, zq, w=mat['source_reduction_mm']), key + '_Q_native_surface')
    we = sep(xp, query_ball_rows(kv, zk, w=mat['wall_mm']), key + '_K_outer_surface')
    globalQ = sep(T[2]['path'], query_ball_rows(qv, zq, w=mat['wall_mm'] + G), key + '_Q_global_native', False)
    if not qs['all_pass'] or not we['all_pass']:
        raise ValueError('SURFACE_COVER_FAILED')
    cav = D / (key + '_cutter.mesh')
    meshwrite(cav, kv, kf)
    dest = D / (key + '_crown.mesh')
    bo = boolean(E[2]['path'], cav, dest)
    (v, f) = meshread(dest)
    pts = v[f].mean(1)
    cd = Distance(kv, kf).query(pts)[0]
    (_, idx, _) = Distance(E[0], E[1]).query(pts)
    roles = np.where(cd < 1e-07, 1, np.where(emask[idx], 0, 2))
    path = D / (key + '.npz')
    np.savez_compressed(path, vertices=v, faces=f, roles=roles, prep_vertices=qv, prep_faces=qf, cavity_vertices=kv, cavity_faces=kf, target=S, M=M, base=b, top=H, center=c, flat_radius=a)
    qpath = D / (key + '_prep_RESEARCH.stl')
    trimesh.Trimesh(qv, qf, process=False).export(qpath)
    return dict(status='GENERATED', mesh_path=str(path), mesh_sha256=sha(path), prep_stl=str(qpath), prep_stl_sha256=sha(qpath), construction=ch, certificates=dict(Q_separated_from_native_surface=qs, cavity_wall_inside_outer=we, Q_dilated_inside_native=globalQ), boolean=bo, source_closure=T[2], outer_closure=E[2], prediction=dict(nominal_gap_mm=0.05, wall_min_mm=mat['wall_mm'], global_native_containment='PASS' if globalQ['all_pass'] else 'FAIL', shape='UNMEASURED', cervical_shoulder='MISSING'), scope='Open-ended cap and model cervical port. Surface separation is not the requested literal full native containment; that separate result is recorded.')

def run():
    rows = []
    pr = read(R / 'PREREG_R3.json')
    st = time.perf_counter()
    for rec in read(R / 'INPUTS.json')['records']:
        try:
            (T, E, S, X, mask) = surfaces(rec)
            err = None
        except Exception as e:
            err = repr(e)
        for mat in pr['materials']:
            t = time.perf_counter()
            r = dict(key=rec['key'], case_key=rec['case_key'], family=rec['family'], fdi=rec['fdi'], material=mat['id'], material_contract=mat, round='R3', outer_route=rec['outer'], resolution='PER_TOOTH')
            try:
                if err:
                    raise ValueError(err)
                (ch, sr) = search3(rec, mat, T, E, S, X)
                r['search'] = sr
                if ch is None:
                    raise ValueError('NO_FEASIBLE_TAPERED_SLEEVE_IN_FROZEN_FAMILY')
                r.update(build(rec, mat, T, E, S, X, mask, ch))
            except Exception as e:
                r.update(status='REJECTED', reason=repr(e))
            r['generation_seconds'] = time.perf_counter() - t
            rows.append(r)
            dump(R / 'raw/R3_GENERATION.json', rows)
            state('R3_GENERATING', dict(key=rec['key'], material=mat['id'], status=r['status'], reason=r.get('reason')), 'Freeze whole-form candidates, evaluate every original link and full native containment')
            print(rec['key'], mat['id'], r['status'], r.get('reason', '')[:110], flush=True)
    lock('FROZEN_PREDICTIONS_R3.json', dict(rows=rows, prereg_sha256=sha(R / 'PREREG_R3.json'), code_sha256=sha(Path(__file__)), seconds=time.perf_counter() - st, physical_predictions=None, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
if __name__ == '__main__':
    run()
