from common import *
from scipy.spatial import ConvexHull

def frames(rec):
    rows = read(PREV / 'raw/CERVICAL_FRAME_DIAGNOSTIC.json')
    out = [np.array(next((r['M'] for r in rows if r['key'] == rec['key']))), np.eye(3)]
    old = read(B / 'PROOF_LANE_CROWN_ROOT_CAUSE/RESULTS_D.json')['rows']
    q = next((r for r in old if r['key'] == rec['key']))
    if q['status'] == 'GENERATED':
        out.append(np.array(q['chosen']['world_to_local']))
    return out

def profile(a, b, H, g=0.0, step=0.15):
    rho = 0.25
    zt = b + rho * (1 - SIN)
    ztop = H - 0.5 + 0.5 * SIN
    if zt >= ztop:
        raise ValueError('FILLETS_OVERLAP')
    R0 = a + offset(np.array([b]), H)[0]
    rc = R0 + rho * (SEC - K)
    zc = b + rho
    rr = rho - g
    angs = np.linspace(-np.pi / 2, -np.arctan(K), 13)
    zb = zc + rr * np.sin(angs)
    rb = rc - rr * np.cos(angs)
    ztopg = H - 0.5 + (0.5 + g) * SIN
    zl = np.linspace(zb[-1], ztopg, max(2, int(np.ceil((ztopg - zb[-1]) / step)) + 1))[1:]
    rl = a + offset(zl, H + g, 0.5 + g)
    angt = np.linspace(np.arctan(K), np.pi / 2, 17)[1:]
    zu = H - 0.5 + (0.5 + g) * np.sin(angt)
    ru = a + (0.5 + g) * np.cos(angt)
    return (np.r_[zb, zl, zu], np.r_[rb, rl, ru])

def make_body(c, a, b, H, M, g=0.0, n=64, extend=False):
    (zs, rs) = profile(a, b, H, g)
    if extend:
        zs = np.r_[b - 10, zs]
        rs = np.r_[rs[0], rs]
    theta = np.arange(n) * 2 * np.pi / n
    xy = np.c_[np.cos(theta), np.sin(theta)]
    v = np.concatenate([np.c_[c + rr * xy, np.full(n, z)] for (rr, z) in zip(rs, zs)]) @ M
    bot = len(v)
    top = bot + 1
    v = np.r_[v, np.array([[*c, zs[0]], [*c, zs[-1]]]) @ M]
    f = []
    for j in range(len(zs) - 1):
        for i in range(n):
            k = (i + 1) % n
            f.extend([[j * n + i, j * n + k, (j + 1) * n + k], [j * n + i, (j + 1) * n + k, (j + 1) * n + i]])
    for i in range(n):
        k = (i + 1) % n
        f.extend([[bot, k, i], [(len(zs) - 1) * n + i, (len(zs) - 1) * n + k, top]])
    return (v, np.array(f), zs, rs)

def propose(rec, mat):
    (v, f, info) = native(rec)
    td = Distance(v, f)
    best = None
    rows = []
    count = 0
    for (mi, M) in enumerate(frames(rec)):
        tv = v @ M.T
        (lo, hi) = (tv.min(0), tv.max(0))
        center = (lo + hi) / 2
        for ex in [0, 0.5, 1]:
            H = hi[2] - mat['top_reduction_mm'] - ex
            for b in np.arange(lo[2] + mat['wall_mm'] + G + 0.05, H - 0.75, 0.25):
                for dx in [-1, 0, 1]:
                    for dy in [-1, 0, 1]:
                        c = center[:2] + [dx, dy]
                        l = 0.15
                        u = min(np.ptp(tv[:, 0]), np.ptp(tv[:, 1])) / 2

                        def feasible(a):
                            try:
                                (qv, qf, zs, rs) = make_body(c, a, b, H, M, n=24)
                            except ValueError:
                                return False
                            xx = np.linspace(-1, 1, 7)
                            disk = np.array([(x, y) for x in xx for y in xx if x * x + y * y <= 1])
                            caps = np.r_[np.c_[c + rs[0] * disk, np.full(len(disk), b)], np.c_[c + a * disk, np.full(len(disk), H)]] @ M
                            pp = np.r_[qv, caps]
                            d = td.query(pp)[0]
                            if d.min() < mat['wall_mm'] + G + 0.025:
                                return False
                            return bool(td.inside(pp).all())
                        ok = feasible(l)
                        count += 1
                        if ok:
                            for _ in range(10):
                                a = (l + u) / 2
                                if feasible(a):
                                    l = a
                                else:
                                    u = a
                            a = np.floor(l * 1000000.0) / 1000000.0
                            (zs, rs) = profile(a, b, H)
                            vol = float(np.pi * np.trapezoid(rs * rs, zs))
                            q = dict(axis_id=mi, M=M, c=c, b=b, H=H, flat_radius_mm=a, volume_proxy_mm3=vol)
                            if best is None or vol > best['volume_proxy_mm3']:
                                best = q
                        rows.append(dict(axis=mi, extra=ex, base=b, dx=dx, dy=dy, feasible=ok, flat_radius=l if ok else None))
    path = D / (rec['key'] + '_' + mat['id'] + '_A_SEARCH.json')
    dump(path, rows)
    return (best, dict(path=str(path), sha256=sha(path), requested=count, rejected=sum((not x['feasible'] for x in rows)), reject_reason='Sampled signed clearance or overlapping fillets; exact acceptance separate'))

def exact_mesh(support, qpath, w, key, surface=False):
    cmd = [str(D / 'exact_mesh_clearance'), str(support), str(qpath), str(w)] + (['surface'] if surface else [])
    z = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
    if z.returncode:
        raise ValueError('EXACT_MESH_EXIT ' + str(z.returncode) + ' ' + z.stderr[:500])
    x = json.loads(z.stdout)
    x.update(support_path=str(support), support_sha256=sha(support), query_path=str(qpath), query_sha256=sha(qpath), required_mm=w, surface_only=surface, command=cmd)
    path = R / 'raw' / (key + '_CERT.json')
    dump(path, x)
    return dict(path=str(path), sha256=sha(path), **x)

def build(rec, mat, ch):
    key = rec['key'] + '_' + mat['id'] + '_A'
    M = np.array(ch['M'])
    c = np.array(ch['c'])
    a = ch['flat_radius_mm']
    b = ch['b']
    H = ch['H']
    (qv, qf, zq, rq) = make_body(c, a, b, H, M)
    qp = D / (key + '_Q.mesh')
    meshwrite(qp, qv, qf)
    T = native(rec)
    cert = exact_mesh(T[2]['path'], qp, mat['wall_mm'] + G, key + '_native')
    if not cert['all_pass']:
        return dict(status='REJECTED_EXACT', construction=ch, containment=cert)
    (kv, kf, zk, rk) = make_body(c, a, b, H, M, g=G, extend=True)
    kp = D / (key + '_cutter.mesh')
    meshwrite(kp, kv, kf)
    info = read(PREV / 'raw' / (rec['key'] + '_R2_SUPPORTS.json'))['outer']
    box = trimesh.creation.box(extents=[100, 100, 100])
    box.vertices += np.array([c[0], c[1], b + G - 50])
    box.vertices = box.vertices @ M
    bp = D / (key + '_cut.mesh')
    meshwrite(bp, box.vertices, box.faces)
    trim = D / (key + '_trim.mesh')
    clip = boolean(info['path'], bp, trim)
    cp = D / (key + '_crown.mesh')
    bo = boolean(trim, kp, cp)
    (v, f) = meshread(cp)
    pts = v[f].mean(1)
    ci = Distance(kv, kf).query(pts)[0]
    roles = np.where(abs(pts @ M[2] - b - G) < 1e-07, 2, np.where(ci < 1e-07, 1, 0))
    path = D / (key + '.npz')
    S = load(rec['private_path'])['target']
    np.savez_compressed(path, vertices=v, faces=f, roles=roles, prep_vertices=qv, prep_faces=qf, cavity_vertices=kv, cavity_faces=kf, target=S, M=M, base=b, top=H, center=c, flat_radius=a)
    return dict(status='GENERATED', mesh_path=str(path), mesh_sha256=sha(path), construction=ch, containment=cert, boolean=bo, clip=clip, source_closure=T[2], outer_closure=info)

def run(limit=None):
    rows = []
    pr = read(R / 'PREREG_A.json')
    st = time.perf_counter()
    for rec in inputs()[:limit]:
        for mat in pr['materials']:
            tick = time.perf_counter()
            r = dict(key=rec['key'], case_key=rec['case_key'], fdi=rec['fdi'], family=rec['family'], material=mat['id'], material_contract=mat, round='A', resolution='PER_TOOTH')
            try:
                (ch, search) = propose(rec, mat)
                r['search'] = search
                if ch is None:
                    r.update(status='REJECTED', reason='NO_CORE_IN_FROZEN_FAMILY')
                else:
                    r.update(build(rec, mat, ch))
            except Exception as e:
                r.update(status='ERROR', error=repr(e))
            r['generation_seconds'] = time.perf_counter() - tick
            rows.append(r)
            dump(R / 'raw/A_GENERATION.json', rows)
            state('A_GENERATING', dict(key=r['key'], material=r['material'], status=r['status'], error=r.get('error')), 'Freeze and assess subset core/cervical crown')
            print(r['key'], r['material'], r['status'], r.get('error', ''), r.get('construction'), flush=True)
    lock('FROZEN_PREDICTIONS_A' + ('_PILOT' if limit else '') + '.json', dict(rows=rows, prereg_sha256=sha(R / 'PREREG_A.json'), prereg_a1_sha256=sha(R / 'PREREG_A1.json'), code_sha256=sha(Path(__file__)), seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
if __name__ == '__main__':
    run(3 if '--pilot' in sys.argv else None)
