from common import *
from height_crown import outer
from loop_crown import loops, align, zipper
from shapely.geometry import Polygon
import shapely

def build(row):
    rec = next((x for x in inputs() if x['key'] == row['key']))
    m = load(row['mesh_path'])
    o = outer(rec)
    M = m['M']
    b = float(m['base'])
    assert np.array_equal(M, o['M']) and b == o['b']
    pv = m['prep_cap_vertices']
    pf = m['prep_cap_faces']
    v = pv.copy()
    f = pf.copy()
    for _ in range(3):
        (v, f) = trimesh.remesh.subdivide(v, f)
    xy = v @ M.T
    poly = Polygon((o['contour'] @ M.T)[:, :2])
    dd = shapely.distance(shapely.points(xy[:, :2]), poly.boundary)
    g = 0.05 * np.minimum(1, dd / 0.3)
    ll = loops(v, f, M)
    assert len(ll) == 1
    bid = ll[0]
    g[bid] = 0
    g[g < 1e-12] = 0
    pd = Distance(pv, pf)
    low = np.zeros(len(v))
    high = np.full(len(v), 2.0)
    for _ in range(32):
        mid = (low + high) / 2
        p = xy.copy()
        p[:, 2] += mid
        dist = pd.query(p @ M)[0]
        sel = dist < g
        low[sel] = mid[sel]
        high[~sel] = mid[~sel]
    local = xy.copy()
    local[:, 2] += np.where(g > 0, low, 0)
    kv = local @ M
    kv[bid] = v[bid]
    center = (v[f].mean(1) @ M.T)[:, 2] - b
    kr = np.where(center > 1e-05, 1, 2)
    phase = align(kv[bid] @ M.T, o['ev'][o['eids']] @ M.T)
    order = np.roll(bid, -phase)
    neck = kv[order] + 0.1 * M[2]
    vv = np.r_[o['ev'], kv, neck]
    ii = order + len(o['ev'])
    ni = np.arange(len(o['ev']) + len(kv), len(vv))
    N = len(ii)
    strip = []
    for i in range(N):
        j = (i + 1) % N
        strip.extend([[ii[i], ii[j], ni[j]], [ii[i], ni[j], ni[i]]])
    collar = zipper(ni, o['eids'], vv, M)
    inner = f + len(o['ev'])
    ff = np.r_[o['ef'], np.array(strip), collar, inner]
    roles = np.r_[np.zeros(len(o['ef']) + len(strip) + len(collar), int), kr]
    mm = trimesh.Trimesh(vv, ff, process=False)
    trimesh.repair.fix_normals(mm, multibody=True)
    if mm.volume < 0:
        mm.invert()
    path = D / (row['key'] + '_' + row['material'] + '_E.npz')
    np.savez_compressed(path, vertices=mm.vertices, faces=mm.faces, roles=roles, prep_vertices=m['prep_vertices'], prep_faces=m['prep_faces'], prep_cap_vertices=pv, prep_cap_faces=pf, cavity_vertices=kv, cavity_faces=f, target=m['target'], M=M, base=b, top=m['top'], center=m['center'], flat_radius=np.nan, margin_vertex_ids=ii, shoulder_vertex_ids=np.unique(inner[kr == 2]), shoulder_gap_mm=(mm.vertices[np.unique(inner[kr == 2])] @ M.T)[:, 2] - b, native_section=kv[order], cavity_nominal_vertices=v, cavity_spacer_target_mm=g)
    return dict(**{k: x for (k, x) in row.items() if k not in ['round', 'mesh_path', 'mesh_sha256', 'generation_seconds']}, round='E', mesh_path=str(path), mesh_sha256=sha(path), parent_mesh_path=row['mesh_path'], parent_mesh_sha256=row['mesh_sha256'], offset=dict(root_iterations=32, root_bracket_max_width_mm=float((high - low).max()), endpoint='lower', refinements=3, neck_height_mm=0.1))

def run():
    rows = []
    st = time.perf_counter()
    for r in read(R / 'FROZEN_PREDICTIONS_D.json')['rows']:
        t = time.perf_counter()
        if r['status'] == 'GENERATED':
            try:
                r = build(r)
            except Exception as e:
                r = dict(**{k: v for (k, v) in r.items() if k not in ['status', 'round']}, round='E', status='REJECTED', offset_error=repr(e))
                __import__('traceback').print_exc()
        else:
            r = dict(**{k: v for (k, v) in r.items() if k != 'round'}, round='E', inherited_D_rejection=True)
        r['generation_seconds'] = time.perf_counter() - t
        rows.append(r)
        dump(R / 'raw/E_GENERATION.json', rows)
        state('E_GENERATING', dict(key=r['key'], material=r['material'], status=r['status']), 'Freeze refinement; exact continuous film and chain checks')
        print(r['key'], r['material'], r['status'], r.get('offset_error'), flush=True)
    lock('FROZEN_PREDICTIONS_E.json', dict(rows=rows, prereg_sha256=sha(R / 'PREREG_E.json'), parent_frozen_sha256=sha(R / 'FROZEN_PREDICTIONS_D.json'), code_sha256=sha(Path(__file__)), seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
if __name__ == '__main__':
    run()
