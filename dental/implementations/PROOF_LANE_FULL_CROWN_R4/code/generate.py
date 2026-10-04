"""Generator imports no scorer and reads only public arrays/template bank."""
from geometry import *
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import spsolve
from scipy.optimize import minimize_scalar

def curve_parameter(cv):
    ell = np.linalg.norm(np.roll(cv, -1, axis=0) - cv, axis=1)
    return (np.r_[0, np.cumsum(ell)], ell.sum())

def curve_eval(cv, u):
    (s, L) = curve_parameter(cv)
    q = np.asarray(u) % 1 * L
    idx = np.minimum(np.searchsorted(s, q, side='right') - 1, len(cv) - 1)
    w = (q - s[idx]) / np.maximum(s[idx + 1] - s[idx], 1e-20)
    return cv[idx] + w[:, None] * (cv[(idx + 1) % len(cv)] - cv[idx])

def correspond(a, b):
    (s, L) = curve_parameter(a)
    u = s[:-1] / L
    best = None
    for flip in [False, True]:
        bb = b[::-1] if flip else b

        def fun(offset):
            return float(np.mean(np.sum((a - curve_eval(bb, u + offset)) ** 2, axis=1)))
        shifts = np.arange(128) / 128
        values = [fun(t) for t in shifts]
        x = shifts[np.argmin(values)]
        opt = minimize_scalar(fun, bounds=(x - 1 / 128, x + 1 / 128), method='bounded', options={'xatol': 1e-12})
        cand = (opt.fun, curve_eval(bb, u + opt.x), dict(flipped=flip, offset=float(opt.x), rms_mm=float(np.sqrt(opt.fun))))
        if best is None or cand[0] < best[0]:
            best = cand
    return (best[1], best[2])

def deform(m, public):
    v = m.vertices.copy()
    f = m.faces
    ll = loops(m)
    if len(ll) != 1:
        raise ValueError('donor boundary count%d' % len(ll))
    ids = ll[0]
    (target, info) = correspond(v[ids], public['margin_curve'])
    d = cKDTree(v[ids]).query(v, workers=1)[0]
    e = m.edges_unique
    rows = np.r_[e[:, 0], e[:, 1]]
    cols = np.r_[e[:, 1], e[:, 0]]
    A = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(v), len(v))).tocsr()
    L = diags(np.asarray(A.sum(1)).ravel()) - A
    fixed = np.unique(np.r_[ids, np.flatnonzero(d >= 1.5)])
    u = np.zeros_like(v)
    u[ids] = target - v[ids]
    free = np.setdiff1d(np.arange(len(v)), fixed)
    mat = (L + diags(np.full(len(v), 0.1))).tocsr()
    for j in range(3):
        u[free, j] = spsolve(mat[free][:, free], -mat[free][:, fixed] @ u[fixed, j])
    out = trimesh.Trimesh(v + u, f, process=False)
    info.update(max_displacement_mm=float(np.linalg.norm(u, axis=1).max()), fixed_high_vertices=int((d >= 1.5).sum()), operation='screened boundary displacement; no injectivity proof')
    return (out, info)

def shell(ext, p):
    ids = loops(ext)
    if len(ids) != 1:
        raise ValueError('cannot stitch multiple boundary')
    ids = ids[0]
    v = ext.vertices.copy()
    N = len(v)
    rim = v[ids]
    (target, ci) = correspond(rim, p['margin_curve'])
    inner = target.copy()
    c = p['cavity_center']
    inner[:, :2] = c + float(p['cavity_scale']) * (inner[:, :2] - c)
    apex = p['cavity_apex']
    verts = np.r_[v, inner, apex[None]]
    faces = ext.faces.tolist()
    roles = [0] * len(faces)
    for (j, a) in enumerate(ids):
        k = (j + 1) % len(ids)
        b = ids[k]
        faces.extend([[int(a), int(b), N + k], [int(a), N + k, N + j]])
        roles.extend([2, 2])
        faces.append([N + j, N + k, N + len(ids)])
        roles.append(1)
    out = trimesh.Trimesh(verts, np.array(faces), process=False)
    trimesh.repair.fix_normals(out, multibody=True)
    if out.volume < 0:
        out.invert()
    return (out, np.array(roles), dict(cavity='fixed nonplanar ruled cone; no measured preparation', apex_above_all_rim=bool(apex[2] > inner[:, 2].max()), inner_boundary_correspondence=ci))

def run(tag='B', only=None):
    st = time.perf_counter()
    lock = read(ROOT / 'FROZEN_INPUTS_B.json')
    rows = []
    for rec in lock['records']:
        if only and rec['key'] != only:
            continue
        for method in ['rigid_control', 'boundary_deform', 'missing_donor']:
            r = dict(key=rec['key'], family=rec['family'], case_key=rec['case_key'], method=method, status='FAILED')
            t0 = time.perf_counter()
            try:
                if rec['status'] != 'PREPARED':
                    raise ValueError(rec['reason'])
                assert sha(rec['public_path']) == rec['public_sha256']
                p = npz(rec['public_path'])
                if method == 'missing_donor':
                    options = []
                    for b in lock['templates']:
                        if b.get('family') != rec['family'] or b['status'] != 'AVAILABLE':
                            continue
                        dat = npz(b['path'])
                        (R, u, errs) = fit(sample(dat['collar'], 1024), sample(p['collar'], 1024))
                        options.append((min(errs), b, dat, R, u, errs))
                    if not options:
                        raise ValueError('no same-FDI templates')
                    (_, b, dat, R, u, errs) = min(options, key=lambda x: x[0])
                    tri = dat['triangles'] @ R.T + u
                    r['template_key'] = b['key']
                    r['template_candidates'] = len(options)
                else:
                    dc = p['donor_collar'].copy()
                    dc[:, :, 0] *= -1
                    dt = p['donor'].copy()
                    dt[:, :, 0] *= -1
                    dt = dt[:, ::-1]
                    (R, u, errs) = fit(sample(dc, 1024), sample(p['collar'], 1024))
                    tri = dt @ R.T + u
                m = mesh(tri, True)
                r['fit'] = dict(R=R, translation_mm=u, start_rms_mm=errs)
                if method != 'rigid_control':
                    (m, r['deformation']) = deform(m, p)
                (out, roles, info) = shell(m, p)
                r['build'] = info
                path = DATA / tag / method / (rec['key'] + '.npz')
                path.parent.mkdir(exist_ok=True, parents=True)
                np.savez_compressed(path, vertices=out.vertices, faces=out.faces, roles=roles)
                r.update(status='GENERATED', mesh_path=path, mesh_sha256=sha(path))
            except Exception as e:
                r['reason'] = repr(e)
            r['seconds'] = time.perf_counter() - t0
            rows.append(r)
            dump(ROOT / 'raw' / f'{tag}_GENERATION.json', rows)
            print(tag, r['key'], method, r['status'], r.get('reason', ''), flush=True)
            state(tag + '_GENERATING', str(len(rows)) + ' predictions written', 'Hash all outputs before scoring held target')
    if only:
        return
    freeze(ROOT / f'FROZEN_PREDICTIONS_{tag}.json', dict(claim_type='capability', input_sha256=sha(ROOT / 'FROZEN_INPUTS_B.json'), prereg_sha256=sha(ROOT / 'PREREG_B.json'), code={str(x): sha(x) for x in sorted((ROOT / 'code').glob('*.py'))}, rows=rows, seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
if __name__ == '__main__':
    run()
