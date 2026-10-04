from common import *
from construction import profile, frames, exact_mesh
from evaluate import boundary_edges
from scipy.ndimage import map_coordinates

def footprint(a, ratio, n):
    t = np.arange(n) * 2 * np.pi / n
    p = np.c_[a * ratio * np.cos(t), a * np.sin(t)]
    nn = np.c_[np.cos(t) / ratio, np.sin(t)]
    nn /= np.linalg.norm(nn, axis=1)[:, None]
    return (p, nn)

def geometry(c, a, ratio, b, H, M, g=0.0, n=96, extend=False, points_only=False):
    (zs, rs) = profile(0, b, H, g, step=0.15)
    if extend:
        zs = np.r_[b - 10, zs]
        rs = np.r_[rs[0], rs]
    (p, nn) = footprint(a, ratio, n)
    vv = np.concatenate([np.c_[c + p + rr * nn, np.full(n, z)] for (z, rr) in zip(zs, rs)]) @ M
    if points_only:
        rr = np.array([0, 0.25, 0.5, 0.75])
        extra = np.concatenate([np.c_[c + factor * (p + rs[0] * nn), np.full(n, b)] for factor in rr] + [np.c_[c + factor * p, np.full(n, H)] for factor in rr]) @ M
        return np.r_[vv, extra]
    bot = len(vv)
    top = bot + 1
    v = np.r_[vv, np.array([[*c, zs[0]], [*c, zs[-1]]]) @ M]
    f = []
    for j in range(len(zs) - 1):
        for i in range(n):
            k = (i + 1) % n
            f.extend([[j * n + i, j * n + k, (j + 1) * n + k], [j * n + i, (j + 1) * n + k, (j + 1) * n + i]])
    for i in range(n):
        k = (i + 1) % n
        f.extend([[bot, k, i], [(len(zs) - 1) * n + i, (len(zs) - 1) * n + k, top]])
    return (v, np.array(f), zs, rs)

class Grid:

    def __init__(self, rec):
        (v, f, info) = native(rec)
        h = 0.15
        lo = v.min(0) - h
        dims = np.ceil((v.max(0) - lo) / h).astype(int) + 2
        N = int(np.prod(dims))
        assert N < 5000000
        path = D / (rec['key'] + '_B_discovery_grid.npz')
        if path.exists():
            x = load(path)
            self.field = x['field']
            self.lo = x['lo']
            self.h = float(x['h'])
        else:
            dd = Distance(v, f)
            val = np.empty(N, np.float32)
            for i in range(0, N, 32768):
                idx = np.arange(i, min(i + 32768, N))
                p = lo + np.stack(np.unravel_index(idx, dims), axis=1) * h
                d = dd.query(p)[0]
                val[idx] = np.where(dd.inside(p), d, -d)
            self.field = val.reshape(dims)
            self.lo = lo
            self.h = h
            np.savez_compressed(path, field=self.field, lo=lo, h=h)
        self.info = dict(path=str(path), sha256=sha(path), nodes=N, pitch_mm=h, scope='Discovery only. Fast winding/interpolation are never inclusion acceptance.')

    def query(self, p):
        return map_coordinates(self.field, ((p - self.lo) / self.h).T, order=1, mode='constant', cval=-100, prefilter=False)

def propose(rec, mat, grid):
    (v, f, _) = native(rec)
    S = load(rec['private_path'])['target']
    (sv, sf) = triangles_mesh(S)
    bd = sv[np.unique(boundary_edges(sf))]
    best = None
    rows = []
    rejects = {}
    count = 0
    for (mi, M0) in enumerate(frames(rec)):
        vv = v @ M0.T
        (evals, evecs) = np.linalg.eigh(np.cov(vv[:, :2].T))
        u = evecs[:, -1]
        rot = np.array([[u[0], u[1], 0], [-u[1], u[0], 0], [0, 0, 1]])
        M = rot @ M0
        tv = v @ M.T
        (lo, hi) = (tv.min(0), tv.max(0))
        center = (lo + hi) / 2
        start = float((bd @ M.T)[:, 2].max() + 0.025)
        for ex in [0, 0.5, 1, 1.5, 2, 2.5, 3]:
            H = hi[2] - mat['top_reduction_mm'] - ex
            for dx in [-1, 0, 1]:
                for dy in [-1, 0, 1]:
                    c = center[:2] + [dx, dy]
                    if grid.query(np.array([[*c, H]]) @ M)[0] < mat['wall_mm'] + G + 0.14:
                        rejects['TOP_POINT_CLEARANCE'] = rejects.get('TOP_POINT_CLEARANCE', 0) + 1
                        continue
                    for b in np.arange(start, H - 0.75, 0.4):
                        for ratio in [1, 2, 4]:
                            count += 1
                            l = 0.05
                            u = min(np.ptp(tv[:, 0]) / ratio, np.ptp(tv[:, 1])) / 2

                            def feasible(a):
                                try:
                                    p = geometry(c, a, ratio, b, H, M, n=24, points_only=True)
                                except ValueError:
                                    return False
                                return grid.query(p).min() >= mat['wall_mm'] + G + 0.14
                            ok = feasible(l)
                            if ok:
                                for _ in range(10):
                                    a = (l + u) / 2
                                    if feasible(a):
                                        l = a
                                    else:
                                        u = a
                                a = np.floor(l * 1000000.0) / 1000000.0
                                (zs, ss) = profile(0, b, H)
                                (p, _) = footprint(a, ratio, 96)
                                per = np.linalg.norm(p - np.roll(p, 1, axis=0), axis=1).sum()
                                vol = float(np.trapezoid(np.pi * a * a * ratio + per * ss + np.pi * ss * ss, zs))
                                ch = dict(axis_id=mi, M=M, c=c, b=b, H=H, flat_radius_mm=a, ellipse_ratio=ratio, volume_proxy_mm3=vol, native_boundary_max_z_mm=start - 0.025)
                                if best is None or vol > best['volume_proxy_mm3']:
                                    best = ch
                            rows.append(dict(axis=mi, extra=ex, dx=dx, dy=dy, b=b, ratio=ratio, feasible=ok, minor_radius_mm=l if ok else None))
    path = D / (rec['key'] + '_' + mat['id'] + '_B_SEARCH.json')
    dump(path, dict(rows=rows, early_rejections=rejects))
    return (best, dict(path=str(path), sha256=sha(path), explicit_candidates=count, rejected=sum((not x['feasible'] for x in rows)), early_rejections=rejects, grid=grid.info))

def segments(tri, z, M):
    t = tri @ M.T
    out = []
    for tt in t:
        p = []
        for (i, j) in [(0, 1), (1, 2), (2, 0)]:
            (a, b) = (tt[i], tt[j])
            if (a[2] - z) * (b[2] - z) < 0:
                p.append(a + (b - a) * ((z - a[2]) / (b[2] - a[2])))
        if len(p) == 2:
            out.append(p)
    return np.array(out)

def ray_radii(seg, c, angles):
    p = seg[:, 0, :2] - c
    e = seg[:, 1, :2] - seg[:, 0, :2]
    out = []
    mult = []
    for angle in angles:
        d = np.array([np.cos(angle), np.sin(angle)])
        den = d[0] * e[:, 1] - d[1] * e[:, 0]
        safe = np.abs(den) > 1e-12
        r = np.full(len(e), np.nan)
        t = r.copy()
        r[safe] = (p[safe, 0] * e[safe, 1] - p[safe, 1] * e[safe, 0]) / den[safe]
        t[safe] = (p[safe, 0] * d[1] - p[safe, 1] * d[0]) / den[safe]
        vals = np.sort(r[safe & (r > 0) & (t >= -1e-08) & (t <= 1 + 1e-08)])
        unique = vals[np.r_[True, np.diff(vals) > 1e-07]] if len(vals) else vals
        mult.append(len(unique))
        out.append(unique[0] if len(unique) else np.nan)
    return (np.array(out), np.array(mult))

def morph_collar(v, f, roles, Touter, S, ch):
    M = np.array(ch['M'])
    c = np.array(ch['c'])
    b = ch['b']
    x = v @ M.T
    old = x.copy()
    newseg = segments(S, b, M)
    oldseg = segments(Touter, b + G, M)
    if len(newseg) == 0 or len(oldseg) == 0:
        raise ValueError('NO_FINISH_LINE_SECTION')
    extids = np.unique(f[roles == 0])
    baseids = np.unique(f[roles == 2])
    ids = np.union1d(extids[x[extids, 2] < b + G + 1 + 1e-07], baseids)
    d = x[ids, :2] - c
    rad = np.linalg.norm(d, axis=1)
    ang = np.arctan2(d[:, 1], d[:, 0])
    (rt, ct) = ray_radii(newseg, c, ang)
    (re, ce) = ray_radii(oldseg, c, ang)
    if np.any(ct != 1) or np.any(ce != 1):
        raise ValueError('FINISH_LINE_NOT_STAR_SHAPED ' + str((np.unique(ct, return_counts=True), np.unique(ce, return_counts=True))))
    basemask = np.isin(ids, baseids)
    outermask = np.isin(ids, extids)
    disp = np.zeros_like(rad)
    disp[outermask] = (rt - re)[outermask] * np.clip(1 - (x[ids[outermask], 2] - b - G), 0, 1)
    innedges = boundary_edges(f[roles == 1])
    ipt = x[np.unique(innedges)]
    ipt = ipt[abs(ipt[:, 2] - b - G) < 1e-06]
    if len(ipt) < 3:
        raise ValueError('NO_INNER_SHOULDER_SEAM')
    ia = np.arctan2(ipt[:, 1] - c[1], ipt[:, 0] - c[0])
    ipt = ipt[np.argsort(ia)]
    iseg = np.stack([ipt, np.roll(ipt, -1, axis=0)], axis=1)
    (ri, ci) = ray_radii(iseg, c, ang)
    if np.any(ci != 1):
        raise ValueError('INNER_SHOULDER_RADIAL_AMBIGUITY')
    u = np.clip((rad - ri) / (re - ri), 0, 1)
    disp[basemask] = u[basemask] * (rt - re)[basemask]
    x[ids, :2] += disp[:, None] * d / rad[:, None]
    x[ids[basemask], 2] = b + G * (1 - u[basemask])
    rim = ids[basemask & outermask]
    x[rim, 2] = b
    return (x @ M, dict(outer_finish_vertices=rim, expected_shoulder_gap_mm=x[baseids, 2] - b, base_vertex_ids=baseids, source_section=newseg @ M, scope='Piecewise-planar native surface intersection; nominal shoulder gap ramps0..50um. Independent distance/mesh checks follow.'))

def build(rec, mat, ch):
    key = rec['key'] + '_' + mat['id'] + '_B'
    M = np.array(ch['M'])
    c = np.array(ch['c'])
    a = ch['flat_radius_mm']
    ratio = ch['ellipse_ratio']
    b = ch['b']
    H = ch['H']
    (qv, qf, zq, rq) = geometry(c, a, ratio, b, H, M)
    qp = D / (key + '_Q.mesh')
    meshwrite(qp, qv, qf)
    T = native(rec)
    cert = exact_mesh(T[2]['path'], qp, mat['wall_mm'] + G, key + '_native')
    if not cert['all_pass']:
        return dict(status='REJECTED_EXACT', construction=ch, containment=cert)
    (kv, kf, zk, rk) = geometry(c, a, ratio, b, H, M, g=G, extend=True)
    kp = D / (key + '_cutter.mesh')
    meshwrite(kp, kv, kf)
    info = read(PREV / 'raw' / (rec['key'] + '_R2_SUPPORTS.json'))['outer']
    (ev, ef) = meshread(info['path'])
    box = trimesh.creation.box(extents=[100, 100, 100])
    box.vertices += np.array([c[0], c[1], b + G - 50])
    box.vertices = box.vertices @ M
    bp = D / (key + '_cut.mesh')
    meshwrite(bp, box.vertices, box.faces)
    trim = D / (key + '_trim.mesh')
    clip = boolean(info['path'], bp, trim)
    cp = D / (key + '_unmorphed.mesh')
    bo = boolean(trim, kp, cp)
    (v, f) = meshread(cp)
    pts = v[f].mean(1)
    ci = Distance(kv, kf).query(pts)[0]
    roles = np.where(abs(pts @ M[2] - b - G) < 1e-07, 2, np.where(ci < 1e-07, 1, 0))
    S = load(rec['private_path'])['target']
    (v, margin) = morph_collar(v, f, roles, ev[ef], S, ch)
    path = D / (key + '.npz')
    np.savez_compressed(path, vertices=v, faces=f, roles=roles, prep_vertices=qv, prep_faces=qf, cavity_vertices=kv, cavity_faces=kf, target=S, M=M, base=b, top=H, center=c, flat_radius=a, ellipse_ratio=ratio, margin_vertex_ids=margin.pop('outer_finish_vertices'), shoulder_vertex_ids=margin.pop('base_vertex_ids'), shoulder_gap_mm=margin.pop('expected_shoulder_gap_mm'), native_section=margin.pop('source_section'))
    return dict(status='GENERATED', mesh_path=str(path), mesh_sha256=sha(path), construction=ch, containment=cert, boolean=bo, clip=clip, source_closure=T[2], outer_closure=info, margin_construction=margin)

def run(limit=None):
    pr = read(R / 'PREREG_A.json')
    rows = []
    st = time.perf_counter()
    pilot = R / 'FROZEN_PREDICTIONS_B_PILOT.json'
    reused = read(pilot)['rows'] if pilot.exists() and limit is None else []
    for rec in inputs()[:limit]:
        old = [x for x in reused if x['key'] == rec['key']]
        if len(old) == 2:
            rows.extend(old)
            print('REUSE', rec['key'], flush=True)
            continue
        tick = time.perf_counter()
        grid = Grid(rec)
        gridseconds = time.perf_counter() - tick
        for mat in pr['materials']:
            tick = time.perf_counter()
            r = dict(key=rec['key'], case_key=rec['case_key'], fdi=rec['fdi'], family=rec['family'], material=mat['id'], material_contract=mat, round='B', resolution='PER_TOOTH', grid_seconds=gridseconds)
            try:
                (ch, sr) = propose(rec, mat, grid)
                r['search'] = sr
                if ch is None:
                    r.update(status='REJECTED', reason='NO_CORE_ABOVE_OBSERVED_FINISH_LINE')
                else:
                    r.update(build(rec, mat, ch))
            except Exception as e:
                r.update(status='ERROR', error=repr(e))
                __import__('traceback').print_exc()
            r['generation_seconds'] = time.perf_counter() - tick
            rows.append(r)
            dump(R / 'raw/B_GENERATION.json', rows)
            state('B_GENERATING', dict(key=r['key'], material=r['material'], status=r['status'], error=r.get('error')), 'Freeze whole assembly candidates then score every gate')
            print(r['key'], r['material'], r['status'], r.get('error'), r.get('construction'), flush=True)
    lock('FROZEN_PREDICTIONS_B' + ('_PILOT' if limit else '') + '.json', dict(rows=rows, prereg_sha256=sha(R / 'PREREG_B.json'), code_sha256=sha(Path(__file__)), seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
if __name__ == '__main__':
    run(3 if '--pilot' in sys.argv else None)
