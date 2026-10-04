from common import *
from construction import frames, exact_mesh
from loop_crown import clip, loops, zipper, align
from scipy.ndimage import gaussian_filter
from scipy.spatial.distance import cdist
import shapely
from shapely.geometry import Polygon

def outer(rec):
    S = load(rec['private_path'])['target']
    (sv0, sf0) = triangles_mesh(S)
    info = read(PREV / 'raw' / (rec['key'] + '_R2_SUPPORTS.json'))['outer']
    (ev0, ef0) = meshread(info['path'])
    bd = sv0[np.unique(__import__('evaluate').boundary_edges(sf0))]
    trials = []
    for (mi, M) in enumerate(frames(rec)):
        for extra in [0.025, 0.1, 0.25, 0.5]:
            b = float((bd @ M.T)[:, 2].max() + extra)
            try:
                (sv, sf) = clip(sv0, sf0, b, M)
                ll = loops(sv, sf, M)
                if len(ll) != 1:
                    raise ValueError('NATIVE_LOOPS_' + str(len(ll)))
                contour = sv[ll[0]]
                poly = Polygon((contour @ M.T)[:, :2])
                if not poly.is_valid:
                    raise ValueError('INVALID_NATIVE_POLYGON')
                for band in [1, 0.5, 0.25]:
                    try:
                        (ev, ef) = clip(ev0, ef0, b + G + band, M)
                        el = loops(ev, ef, M)
                        if len(el) != 1:
                            continue
                        phase = align(contour @ M.T, ev[el[0]] @ M.T)
                        contour = np.roll(contour, -phase, axis=0)
                        v = np.r_[ev, contour]
                        co = np.arange(len(ev), len(v))
                        ff = np.r_[ef, zipper(co, el[0], v, M)]
                        mm = trimesh.Trimesh(v, ff, process=False)
                        if len(trimesh.graph.connected_components(mm.face_adjacency, nodes=np.arange(len(ff)))) != 1:
                            continue
                        return dict(M=M, b=b, ev=ev, ef=ef, eids=el[0], contour=contour, outer_vertices=v, outer_faces=ff, axis_id=mi, band=band, source=S, info=info, trials=trials)
                    except Exception as e:
                        trials.append(dict(axis=mi, extra=extra, band=band, error=repr(e)))
            except Exception as e:
                trials.append(dict(axis=mi, extra=extra, error=repr(e)))
    raise ValueError('NO_SIMPLE_SOURCE_AND_UPPER_CAP ' + str(trials))

def field(rec, mat, o):
    M = o['M']
    b = o['b']
    poly = Polygon((o['contour'] @ M.T)[:, :2])
    lo = np.array(poly.bounds[:2]) - 0.4
    hi = np.array(poly.bounds[2:]) + 0.4
    h = 0.2
    (nx, ny) = np.ceil((hi - lo) / h).astype(int) + 1
    (xx, yy) = np.meshgrid(lo[0] + h * np.arange(nx), lo[1] + h * np.arange(ny), indexing='ij')
    xy = np.c_[xx.ravel(), yy.ravel()]
    inside = shapely.contains_xy(poly, xy[:, 0], xy[:, 1])
    (tv, tf, _) = native(rec)
    td = Distance(tv, tf)
    od = Distance(o['outer_vertices'], o['outer_faces'])
    top = min((tv @ M.T)[:, 2].max(), (o['ev'] @ M.T)[:, 2].max())
    zs = np.arange(b, top + 0.1, 0.1)
    raw = np.zeros(len(xy))
    active = inside.copy()
    required = mat['wall_mm'] + G + 0.2
    for z in zs:
        ids = np.flatnonzero(active)
        if len(ids) == 0:
            break
        p = np.c_[xy[ids], np.full(len(ids), z)] @ M
        good = np.minimum(td.query(p)[0], od.query(p)[0]) >= required
        good &= td.inside(p)
        active[ids[~good]] = False
        raw[ids[good]] = z - b
    L = 1 / (np.sqrt(2) * 0.12)
    env = np.empty_like(raw)
    for i in range(0, len(xy), 128):
        env[i:i + 128] = np.min(raw[None, :] + L * cdist(xy[i:i + 128], xy), axis=1)
    env = env.reshape(nx, ny)
    for _ in range(3):
        env = np.minimum(env, gaussian_filter(env, 1, mode='constant', cval=0))
    env[env < 1e-06] = 0
    path = D / (rec['key'] + '_' + mat['id'] + '_D_field.npz')
    np.savez_compressed(path, xy=xy, height=env, raw_height=raw, lo=lo, pitch=h)
    if env.max() < 0.75:
        raise ValueError('NATIVE_FIELD_PEAK_BELOW_0.75 ' + str(env.max()))
    verts = []
    faces = []
    vmap = {}
    heights = []

    def vertex(p, z):
        key = tuple(np.round(p, 10))
        if key not in vmap:
            vmap[key] = len(verts)
            verts.append([p[0], p[1], b + z])
            heights.append(z)
        return vmap[key]
    flat = env.ravel()
    for i in range(nx - 1):
        for j in range(ny - 1):
            a = i * ny + j
            bb = (i + 1) * ny + j
            c = bb + 1
            d = a + 1
            for ids in [[a, bb, c], [a, c, d]]:
                p = xy[ids]
                tri = Polygon(p)
                if not tri.intersects(poly):
                    continue
                part = tri.intersection(poly)
                geoms = [part] if part.geom_type == 'Polygon' else list(part.geoms) if part.geom_type == 'MultiPolygon' else []
                coeff = np.linalg.solve(np.c_[p, np.ones(3)], flat[ids])
                for piece in geoms:
                    if piece.area < 1e-15:
                        continue
                    for t in shapely.constrained_delaunay_triangles(piece).geoms:
                        pts = np.asarray(t.exterior.coords)[:3]
                        zz = np.maximum(0, np.c_[pts, np.ones(3)] @ coeff)
                        faces.append([vertex(q, z) for (q, z) in zip(pts, zz)])
    v = np.array(verts) @ M
    f = np.array(faces)
    hh = np.array(heights)
    mm = trimesh.Trimesh(v, f, process=False)
    trimesh.repair.fix_normals(mm, multibody=True)
    if np.median(mm.face_normals @ M[2]) < 0:
        mm.faces = mm.faces[:, ::-1]
    f = mm.faces.copy()
    positive = hh[f].max(1) > 1e-08
    qf = f[positive]
    used = np.unique(qf)
    qv = v[used]
    remap = np.full(len(v), -1, int)
    remap[used] = np.arange(len(used))
    qf = remap[qf]
    qheight = hh[used]
    bv = qv.copy()
    bv = bv @ M.T
    bv[:, 2] = b
    bv = bv @ M
    bottom = np.empty(len(qv), int)
    vv = list(qv)
    for k in range(len(qv)):
        if qheight[k] <= 1e-08:
            bottom[k] = k
        else:
            bottom[k] = len(vv)
            vv.append(bv[k])
    allf = np.r_[qf, bottom[qf[:, ::-1]]]
    qm = trimesh.Trimesh(np.array(vv), allf, process=False)
    trimesh.repair.fix_normals(qm, multibody=True)
    if qm.volume < 0:
        qm.invert()
    return (v, f, hh, qm, dict(path=str(path), sha256=sha(path), nodes=len(xy), raw_max_height_mm=raw.max(), max_height_mm=env.max(), slope_L=L, grid_mm=h, retained_Q_volume_mm3=float(qm.volume), scope='Native-adaptive PL graph, discovery buffer0.2mm; exact mesh clearance follows'))

def cavity(pv, pf, hh, o):
    v = pv.copy()
    f = pf.copy()
    for _ in range(2):
        (v, f) = trimesh.remesh.subdivide(v, f)
    M = o['M']
    b = o['b']
    local = v @ M.T
    poly = Polygon((o['contour'] @ M.T)[:, :2])
    distmargin = shapely.distance(shapely.points(local[:, :2]), poly.boundary)
    g = 0.05 * np.minimum(1, distmargin / 0.3)
    g[g < 1e-09] = 0
    pd = Distance(pv, pf)
    low = np.zeros(len(v))
    high = np.full(len(v), 2.0)
    for _ in range(16):
        mid = (low + high) / 2
        p = local.copy()
        p[:, 2] += mid
        dd = pd.query(p @ M)[0]
        sel = dd < g
        low[sel] = mid[sel]
        high[~sel] = mid[~sel]
    local[:, 2] += np.where(g > 0, (low + high) / 2, 0)
    kv = local @ M
    center = (v[f].mean(1) @ M.T)[:, 2] - b
    roles = np.where(center > 1e-05, 1, 2)
    return (kv, f, roles, v, g)

def build(rec, mat, o):
    key = rec['key'] + '_' + mat['id'] + '_D'
    M = o['M']
    b = o['b']
    (pv, pf, hh, qm, fi) = field(rec, mat, o)
    certs = []
    for (k, component) in enumerate(qm.split(only_watertight=False)):
        if not component.is_watertight:
            raise ValueError('Q_COMPONENT_OPEN')
        qp = D / (key + '_Q_component' + str(k) + '.mesh')
        meshwrite(qp, component.vertices, component.faces)
        certs.append(exact_mesh(native(rec)[2]['path'], qp, mat['wall_mm'] + G, key + '_native_' + str(k)))
    containment = dict(all_pass=all((x['all_pass'] for x in certs)), components=certs, scope='Union of all certified Q components; no component dropped')
    if not containment['all_pass']:
        return dict(status='REJECTED_EXACT', containment=containment, field=fi)
    (kv, kf, kr, nominalv, g) = cavity(pv, pf, hh, o)
    kl = loops(kv, kf, M)
    if len(kl) != 1:
        raise ValueError('CAVITY_BOUNDARY_LOOPS ' + str(len(kl)))
    contour = kv[kl[0]]
    phase = align(contour @ M.T, o['ev'][o['eids']] @ M.T)
    ki = np.roll(kl[0], -phase)
    v = np.r_[o['ev'], kv]
    ii = ki + len(o['ev'])
    innerf = kf + len(o['ev'])
    collar = zipper(ii, o['eids'], v, M)
    f = np.r_[o['ef'], collar, innerf]
    roles = np.r_[np.zeros(len(o['ef']) + len(collar), int), kr]
    mm = trimesh.Trimesh(v, f, process=False)
    trimesh.repair.fix_normals(mm, multibody=True)
    if mm.volume < 0:
        mm.invert()
    path = D / (key + '.npz')
    np.savez_compressed(path, vertices=mm.vertices, faces=mm.faces, roles=roles, prep_vertices=qm.vertices, prep_faces=qm.faces, prep_cap_vertices=pv, prep_cap_faces=pf, cavity_vertices=kv, cavity_faces=kf, target=o['source'], M=M, base=b, top=b + fi['max_height_mm'], center=np.mean((pv @ M.T)[:, :2], axis=0), flat_radius=np.nan, margin_vertex_ids=ii, shoulder_vertex_ids=np.unique(innerf[kr == 2]), shoulder_gap_mm=(mm.vertices[np.unique(innerf[kr == 2])] @ M.T)[:, 2] - b, native_section=contour, cavity_nominal_vertices=nominalv, cavity_spacer_target_mm=g)
    return dict(status='GENERATED', mesh_path=str(path), mesh_sha256=sha(path), containment=containment, field=fi, source_closure=native(rec)[2], outer_closure=o['info'], construction=dict(M=M, b=b, H=b + fi['max_height_mm'], axis_id=o['axis_id'], retained_Q_volume_mm3=fi['retained_Q_volume_mm3'], representation='Native-adaptive slope-bounded height field'), collar=dict(height_mm=o['band'], trials=o['trials']))

def run(limit=None):
    rows = []
    st = time.perf_counter()
    pilot = R / 'FROZEN_PREDICTIONS_D_PILOT.json'
    reuse = read(pilot)['rows'] if pilot.exists() and limit is None else []
    for rec in inputs()[:limit]:
        old = [x for x in reuse if x['key'] == rec['key']]
        if len(old) == 2:
            rows.extend(old)
            print('REUSE', rec['key'], flush=True)
            continue
        try:
            o = outer(rec)
            outererr = None
        except Exception as e:
            outererr = repr(e)
        for mat in read(R / 'PREREG_A.json')['materials']:
            tick = time.perf_counter()
            r = dict(key=rec['key'], case_key=rec['case_key'], fdi=rec['fdi'], family=rec['family'], material=mat['id'], material_contract=mat, round='D', resolution='PER_TOOTH')
            try:
                if outererr:
                    raise ValueError(outererr)
                r.update(build(rec, mat, o))
            except Exception as e:
                r.update(status='REJECTED', error=repr(e))
                __import__('traceback').print_exc()
            r['generation_seconds'] = time.perf_counter() - tick
            rows.append(r)
            dump(R / 'raw/D_GENERATION.json', rows)
            state('D_GENERATING', dict(key=r['key'], material=r['material'], status=r['status'], error=r.get('error')), 'Freeze adaptive-cap outputs; exact full-chain assessment')
            print(r['key'], r['material'], r['status'], r.get('error'), r.get('field'), flush=True)
    lock('FROZEN_PREDICTIONS_D' + ('_PILOT' if limit else '') + '.json', dict(rows=rows, prereg_sha256=sha(R / 'PREREG_D.json'), code_sha256=sha(Path(__file__)), seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
if __name__ == '__main__':
    run(3 if '--pilot' in sys.argv else None)
