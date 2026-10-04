from generate import curve_parameter, curve_eval, correspond
from geometry import *
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import spsolve

def exact_trace(m, curve):
    v = m.vertices.copy()
    f = m.faces.copy()
    ids = loops(m)[0]
    a = v[ids]
    (_, info) = correspond(a, curve)
    bb = curve[::-1] if info['flipped'] else curve
    offset = info['offset']
    (s, L) = curve_parameter(a)
    u = s[:-1] / L
    (tb, LB) = curve_parameter(bb)
    need = (tb[:-1] / LB - offset) % 1
    verts = v.tolist()
    desired = {int(i): q for (i, q) in zip(ids, curve_eval(bb, u + offset))}
    add = {}
    inserted = []
    for (j, un) in enumerate(need):
        k = min(int(np.searchsorted(u, un, side='right') - 1), len(ids) - 1)
        k = max(k, 0)
        uu = un
        end = u[k + 1] if k + 1 < len(u) else 1.0
        frac = (uu - u[k]) / (end - u[k])
        aa = int(ids[k])
        b = int(ids[(k + 1) % len(ids)])
        idx = len(verts)
        verts.append(((1 - frac) * v[aa] + frac * v[b]).tolist())
        desired[idx] = bb[j].copy()
        add.setdefault((aa, b), []).append((frac, idx))
        inserted.append(idx)
    for edge in add:
        add[edge].sort()
    faces = []
    for tri in f:
        poly = []
        changed = False
        for (aa, b) in zip(tri, np.roll(tri, -1)):
            poly.append(int(aa))
            ls = add.get((int(aa), int(b)))
            if ls is not None:
                poly.extend((i for (_, i) in ls))
                changed = True
            elif (int(b), int(aa)) in add:
                poly.extend((i for (_, i) in add[int(b), int(aa)][::-1]))
                changed = True
        if not changed:
            faces.append(tri.tolist())
            continue
        mid = len(verts)
        verts.append(v[tri].mean(0).tolist())
        faces.extend([[mid, aa, b] for (aa, b) in zip(poly, np.roll(poly, -1))])
    vv = np.array(verts)
    ff = np.array(faces)
    mm = trimesh.Trimesh(vv, ff, process=False)
    e = mm.edges_unique
    ii = np.r_[e[:, 0], e[:, 1]]
    jj = np.r_[e[:, 1], e[:, 0]]
    A = coo_matrix((np.ones(len(ii)), (ii, jj)), shape=(len(vv), len(vv))).tocsr()
    L = diags(np.asarray(A.sum(1)).ravel()) - A
    bound = np.array(sorted(desired))
    d = cKDTree(vv[bound]).query(vv, workers=1)[0]
    fixed = np.unique(np.r_[bound, np.flatnonzero(d >= 1.5)])
    free = np.setdiff1d(np.arange(len(vv)), fixed)
    disp = np.zeros_like(vv)
    for (i, q) in desired.items():
        disp[i] = q - vv[i]
    mat = (L + diags(np.full(len(vv), 0.1))).tocsr()
    for j in range(3):
        disp[free, j] = spsolve(mat[free][:, free], -mat[free][:, fixed] @ disp[fixed, j])
    out = trimesh.Trimesh(vv + disp, ff, process=False)
    for (i, q) in desired.items():
        out.vertices[i] = q
    return (out, dict(**info, inserted_target_knots=len(inserted), target_knot_identity_error_mm=float(np.max(np.linalg.norm(out.vertices[inserted] - bb, axis=1))), curve_error_mm=curve_error(out.vertices[loops(out)[0]], curve), max_displacement_mm=float(np.linalg.norm(disp, axis=1).max())))

def exact_shell(ext, p):
    ids = loops(ext)[0]
    v = ext.vertices.copy()
    N = len(v)
    rim = v[ids]
    inner = rim.copy()
    c = p['cavity_center']
    inner[:, :2] = c + float(p['cavity_scale']) * (inner[:, :2] - c)
    apex = p['cavity_apex']
    verts = np.r_[v, inner, apex[None]]
    faces = ext.faces.tolist()
    roles = [0] * len(faces)
    for (j, a) in enumerate(ids):
        k = (j + 1) % len(ids)
        b = ids[k]
        faces.extend([[int(a), int(b), N + k], [int(a), N + k, N + j], [N + j, N + k, N + len(ids)]])
        roles.extend([2, 2, 1])
    out = trimesh.Trimesh(verts, np.array(faces), process=False)
    trimesh.repair.fix_normals(out, multibody=True)
    if out.volume < 0:
        out.invert()
    return (out, np.array(roles))

def wall_repair(ext, p):
    base = ext.vertices.copy()
    (m, roles) = exact_shell(ext, p)
    inner = m.triangles[roles == 1]
    (q, d, j) = closest(inner, p['margin_curve'])
    k = int(np.argmin(d))
    witness = dict(margin_point=p['margin_curve'][k], inner_closest_point=q[k], distance_mm=float(d[k]), largest_allowed_wall_mm=float(d[k] + 0.025), infeasible=bool(d[k] + 0.025 < 0.5), proof='distance1-Lipschitz necessary condition; floating evaluation, roundoff enclosure missing')
    if witness['infeasible']:
        return (ext, dict(status='ABSTAIN_INFEASIBLE_FIXED_MARGIN', witness=witness))
    ids = loops(ext)[0]
    free = np.ones(len(base), bool)
    free[ids] = False
    history = []
    out = ext.copy()
    norm = np.cross(inner[:, 1] - inner[:, 0], inner[:, 2] - inner[:, 0])
    norm /= np.maximum(np.linalg.norm(norm, axis=1)[:, None], 1e-20)
    norm[norm[:, 2] < 0] *= -1
    for _ in range(12):
        (q, d, j) = closest(inner, out.vertices)
        signed = np.einsum('ij,ij->i', out.vertices - q, norm[j])
        active = free & (signed < 0.65)
        step = np.zeros_like(base)
        step[active] = (0.65 - signed[active, None]) * norm[j[active]]
        desired = out.vertices + step
        total = desired - base
        length = np.linalg.norm(total, axis=1)
        total *= np.minimum(1, 0.5 / np.maximum(length, 1e-20))[:, None]
        out.vertices = base + total
        history.append(float(d.min()))
        if not active.any():
            break
    return (out, dict(status='BOUNDED_LOCAL_PROJECTION', witness=witness, minimum_distance_history_mm=history, additional_max_displacement_mm=float(np.linalg.norm(out.vertices - base, axis=1).max())))

def run():
    st = time.perf_counter()
    lock = {r['key']: r for r in read(ROOT / 'FROZEN_INPUTS_B.json')['records']}
    prior = read(ROOT / 'FROZEN_PREDICTIONS_B.json')
    rows = []
    for rec in prior['rows']:
        if rec['method'] != 'rigid_control':
            continue
        p = npz(lock[rec['key']]['public_path'])
        raw = npz(rec['mesh_path'])
        ext = mesh(raw['vertices'][raw['faces'][raw['roles'] == 0]], True)
        for method in ['exact_margin', 'fixed_prep_projection']:
            r = {k: rec[k] for k in ['key', 'case_key', 'family']}
            r.update(method=method, status='FAILED')
            t0 = time.perf_counter()
            try:
                (e, info) = exact_trace(ext, p['margin_curve'])
                r['trace'] = info
                if method == 'fixed_prep_projection':
                    (e, r['repair']) = wall_repair(e, p)
                (m, roles) = exact_shell(e, p)
                path = DATA / 'C' / method / (r['key'] + '.npz')
                path.parent.mkdir(parents=True, exist_ok=True)
                np.savez_compressed(path, vertices=m.vertices, faces=m.faces, roles=roles)
                r.update(status='GENERATED', mesh_path=path, mesh_sha256=sha(path), build=dict(apex_above_all_rim=bool(p['cavity_apex'][2] > p['margin_curve'][:, 2].max())))
            except Exception as e:
                r['reason'] = repr(e)
            r['seconds'] = time.perf_counter() - t0
            rows.append(r)
            dump(ROOT / 'raw/C_GENERATION.json', rows)
            print('C', r['key'], method, r['status'], r.get('reason', ''), flush=True)
    freeze(ROOT / 'FROZEN_PREDICTIONS_C.json', dict(claim_type='capability', prereg_sha256=sha(ROOT / 'PREREG_C.json'), input_sha256=sha(ROOT / 'FROZEN_INPUTS_B.json'), code_sha256=sha(Path(__file__)), rows=rows, seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
if __name__ == '__main__':
    run()
