from common import *
from elliptic import Grid, propose, geometry
from construction import exact_mesh
from evaluate import boundary_edges
import shapely
from shapely.geometry import Polygon

def clip(v, f, z, M, above=True):
    xx = v @ M.T
    vv = list(v)
    out = []
    cache = {}
    sign = 1 if above else -1

    def crossing(i, j):
        key = tuple(sorted((int(i), int(j))))
        if key not in cache:
            (a, b) = key
            t = (z - xx[a, 2]) / (xx[b, 2] - xx[a, 2])
            cache[key] = len(vv)
            vv.append(v[a] + t * (v[b] - v[a]))
        return cache[key]
    for face in f:
        poly = []
        for (i, j) in zip(face, np.roll(face, -1)):
            ai = sign * (xx[i, 2] - z) >= 0
            aj = sign * (xx[j, 2] - z) >= 0
            if ai:
                poly.append(int(i))
            if ai != aj:
                poly.append(crossing(i, j))
        for k in range(1, len(poly) - 1):
            out.append([poly[0], poly[k], poly[k + 1]])
    if not out:
        raise ValueError('EMPTY_CLIP')
    return compact(np.array(vv), np.array(out))

def loops(v, f, M):
    e = boundary_edges(f)
    adj = {}
    for (a, b) in e:
        adj.setdefault(int(a), []).append(int(b))
        adj.setdefault(int(b), []).append(int(a))
    if any((len(a) != 2 for a in adj.values())):
        raise ValueError('SECTION_NOT_TWO_REGULAR')
    unseen = set(adj)
    out = []
    x = v @ M.T
    while unseen:
        first = min(unseen)
        cur = first
        prev = None
        seq = []
        while True:
            seq.append(cur)
            unseen.discard(cur)
            nxt = adj[cur][0] if adj[cur][0] != prev else adj[cur][1]
            (prev, cur) = (cur, nxt)
            if cur == first:
                break
            if cur not in unseen:
                raise ValueError('SECTION_NON_SIMPLE_INCIDENCE')
        ids = np.array(seq)
        p = x[ids, :2]
        area = np.sum(p[:, 0] * np.roll(p[:, 1], -1) - p[:, 1] * np.roll(p[:, 0], -1)) / 2
        if area < 0:
            ids = ids[::-1]
        out.append(ids)
    return out

def arclen(p):
    d = np.linalg.norm(np.roll(p, -1, axis=0) - p, axis=1)
    return np.r_[0, np.cumsum(d)] / d.sum()

def interp(p, t):
    arc = arclen(p)
    return np.stack([np.interp(t, arc, np.r_[p[:, k], p[0, k]]) for k in range(p.shape[1])], axis=1)

def align(a, b):
    sampleb = interp(b, np.arange(128) / 128)
    best = None
    for i in np.unique(np.linspace(0, len(a) - 1, min(len(a), 128)).astype(int)):
        aa = np.roll(a, -i, axis=0)
        err = np.sum((interp(aa, np.arange(128) / 128)[:, :2] - sampleb[:, :2]) ** 2)
        if best is None or err < best[0]:
            best = (err, i)
    return int(best[1])

def zipper(aidx, bidx, v, M):
    a = v[aidx] @ M.T
    b = v[bidx] @ M.T
    ta = arclen(a[:, :2])
    tb = arclen(b[:, :2])
    i = j = 0
    out = []
    na = len(aidx)
    nb = len(bidx)
    while i < na or j < nb:
        if j == nb or (i < na and ta[i + 1] <= tb[j + 1]):
            out.append([aidx[i % na], aidx[(i + 1) % na], bidx[j % nb]])
            i += 1
        else:
            out.append([aidx[i % na], bidx[(j + 1) % nb], bidx[j % nb]])
            j += 1
    return np.array(out)

def annulus(v, outer, inner, M):
    x = v @ M.T
    p = Polygon(x[outer, :2], holes=[x[inner, :2]])
    if not p.is_valid:
        raise ValueError('INVALID_SHOULDER_POLYGON ' + shapely.is_valid_reason(p))
    mapping = {tuple(x[i, :2]): int(i) for i in np.r_[outer, inner]}
    tris = shapely.constrained_delaunay_triangles(p)
    out = []
    for tri in tris.geoms:
        q = np.asarray(tri.exterior.coords)[:3]
        out.append([mapping[tuple(y)] for y in q])
    out = np.array(out)
    area = float(sum((t.area for t in tris.geoms)))
    if abs(area - p.area) > 1e-08:
        raise ValueError('SHOULDER_AREA_PARTITION')
    return (out, dict(triangles=len(out), polygon_area_mm2=p.area, triangulated_area_mm2=area))

def build(rec, mat, ch):
    key = rec['key'] + '_' + mat['id'] + '_C'
    M = np.array(ch['M'])
    c = np.array(ch['c'])
    a = ch['flat_radius_mm']
    ratio = ch['ellipse_ratio']
    b = ch['b']
    H = ch['H']
    (qv, qf, _, _) = geometry(c, a, ratio, b, H, M)
    qp = D / (key + '_Q.mesh')
    meshwrite(qp, qv, qf)
    T = native(rec)
    cert = exact_mesh(T[2]['path'], qp, mat['wall_mm'] + G, key + '_native')
    if not cert['all_pass']:
        return dict(status='REJECTED_EXACT', construction=ch, containment=cert)
    S = load(rec['private_path'])['target']
    (sv, sf) = triangles_mesh(S)
    (sv, sf) = clip(sv, sf, b, M)
    sl = loops(sv, sf, M)
    if len(sl) != 1:
        raise ValueError('NATIVE_FINISH_LOOPS ' + str(len(sl)))
    contour = sv[sl[0]]
    info = read(PREV / 'raw' / (rec['key'] + '_R2_SUPPORTS.json'))['outer']
    (ev0, ef0) = meshread(info['path'])
    selected = None
    trials = []
    for h in [1, 0.5, 0.25]:
        try:
            (ev, ef) = clip(ev0, ef0, b + G + h, M)
            el = loops(ev, ef, M)
            nc = len(trimesh.graph.connected_components(trimesh.Trimesh(ev, ef, process=False).face_adjacency, nodes=np.arange(len(ef))))
            trials.append(dict(height_mm=h, boundary_loops=len(el), components=nc))
            if len(el) == 1 and nc == 1:
                selected = (ev, ef, el[0], h)
                break
        except Exception as e:
            trials.append(dict(height_mm=h, error=repr(e)))
    if selected is None:
        raise ValueError('NO_SINGLE_UPPER_CAP ' + str(trials))
    (ev, ef, eids, h) = selected
    phase = align(contour @ M.T, ev[eids] @ M.T)
    contour = np.roll(contour, -phase, axis=0)
    (kv, kf, _, _) = geometry(c, a, ratio, b, H, M, g=G)
    klocal = kv[kf] @ M.T
    basal = abs(klocal[:, :, 2] - b - G).max(1) < 1e-07
    kf = kf[~basal]
    il = loops(kv, kf, M)
    if len(il) != 1:
        raise ValueError('INNER_LOOPS ' + str(len(il)))
    v = np.r_[ev, contour, kv]
    co = np.arange(len(ev), len(ev) + len(contour))
    ii = il[0] + len(ev) + len(contour)
    innerf = kf + len(ev) + len(contour)
    collar = zipper(co, eids, v, M)
    (shoulder, shinfo) = annulus(v, co, ii, M)
    f = np.r_[ef, collar, innerf, shoulder]
    roles = np.r_[np.zeros(len(ef) + len(collar), int), np.ones(len(innerf), int), np.full(len(shoulder), 2, int)]
    mm = trimesh.Trimesh(v, f, process=False)
    trimesh.repair.fix_normals(mm, multibody=True)
    if mm.volume < 0:
        mm.invert()
    (v, f) = (mm.vertices.copy(), mm.faces.copy())
    path = D / (key + '.npz')
    np.savez_compressed(path, vertices=v, faces=f, roles=roles, prep_vertices=qv, prep_faces=qf, cavity_vertices=kv, cavity_faces=kf, target=S, M=M, base=b, top=H, center=c, flat_radius=a, ellipse_ratio=ratio, margin_vertex_ids=co, shoulder_vertex_ids=np.r_[co, ii], shoulder_gap_mm=(v[np.r_[co, ii]] @ M.T)[:, 2] - b, native_section=contour)
    return dict(status='GENERATED', mesh_path=str(path), mesh_sha256=sha(path), construction=ch, containment=cert, source_closure=T[2], outer_closure=info, collar=dict(height_mm=h, trials=trials, phase=phase, native_vertices=len(co), cap_vertices=len(eids)), shoulder_triangulation=shinfo)

def run(limit=None):
    pr = read(R / 'PREREG_A.json')
    rows = []
    st = time.perf_counter()
    pilot = R / 'FROZEN_PREDICTIONS_C_PILOT.json'
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
            r = dict(key=rec['key'], case_key=rec['case_key'], fdi=rec['fdi'], family=rec['family'], material=mat['id'], material_contract=mat, round='C', resolution='PER_TOOTH', grid_seconds=gridseconds)
            try:
                (ch, sr) = propose(rec, mat, grid)
                r['search'] = sr
                r['construction'] = ch
                if ch is None:
                    r.update(status='REJECTED', reason='NO_CORE_ABOVE_OBSERVED_FINISH_LINE')
                else:
                    r.update(build(rec, mat, ch))
            except Exception as e:
                r.update(status='ERROR', error=repr(e))
                __import__('traceback').print_exc()
            r['generation_seconds'] = time.perf_counter() - tick
            rows.append(r)
            dump(R / 'raw/C_GENERATION.json', rows)
            state('C_GENERATING', dict(key=r['key'], material=r['material'], status=r['status'], error=r.get('error')), 'Freeze then assess source-contour assembly')
            print(r['key'], r['material'], r['status'], r.get('error'), flush=True)
    lock('FROZEN_PREDICTIONS_C' + ('_PILOT' if limit else '') + '.json', dict(rows=rows, prereg_sha256=sha(R / 'PREREG_C.json'), code_sha256=sha(Path(__file__)), seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
if __name__ == '__main__':
    run(3 if '--pilot' in sys.argv else None)
