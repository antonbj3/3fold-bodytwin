from core import *
from scipy.spatial import cKDTree

def membership(tri, queries, wall):
    scale = 10 ** 9
    if max(abs(tri).max(), abs(queries[:, :3]).max()) > 1000000.0:
        raise ValueError('INTEGER_ENCLOSURE_DOMAIN')
    ids = cKDTree(queries[:, :3]).query(tri.mean(1), k=min(5, len(queries)))[1]
    if ids.ndim == 1:
        ids = ids[:, None]
    lo = np.floor(tri * scale).astype(np.int64) - 1
    hi = np.ceil(tri * scale).astype(np.int64) + 1
    clo = np.floor(queries[:, :3] * scale).astype(np.int64) - 1
    chi = np.ceil(queries[:, :3] * scale).astype(np.int64) + 1
    radii = np.array([int((F(float(x)) - F(str(wall))) * scale) for x in queries[:, 3]], dtype=object)
    ok = np.zeros(len(tri), bool)
    chosen = np.full(len(tri), -1, int)
    for k in range(ids.shape[1]):
        jj = ids[:, k]
        upper = np.maximum(abs(lo - chi[jj, None, :]), abs(hi - clo[jj, None, :])).astype(object)
        d2 = (upper * upper).sum(2).max(1)
        fits = np.asarray(d2 <= radii[jj] ** 2, dtype=bool)
        chosen[~ok & fits] = jj[~ok & fits]
        ok |= fits
    return dict(all_covered=bool(ok.all()), faces=len(tri), uncovered=int((~ok).sum()), cover_index=chosen, coordinate_interval_tick_mm=1e-09, scope='Rigorous integer coordinate enclosure under explicit|coordinate|<=1e6mm; squared norms use unbounded integers, no floating sensitivity estimate')

def prove(ext, inn, row):
    from sleeve import sep
    old = read(row['certificates']['cavity_wall_inside_outer']['path'])
    qq = np.loadtxt(old['query_path'], skiprows=1)
    cov = membership(inn, qq, row['material_contract']['wall_mm'])
    (v, f) = triangles_mesh(ext)
    key = Path(row['mesh_path']).stem + '_ACTUAL_WALL'
    p = D / (key + '_outer.mesh')
    meshwrite(p, v, f)
    try:
        cert = sep(p, qq, key, surface=True)
    except Exception as e:
        cert = dict(all_pass=False, error=repr(e))
    out = dict(all_pass=bool(cov['all_covered'] and cert['all_pass']), actual_outer_separation=cert, actual_inner_triangle_cover=cov, wall_mm=row['material_contract']['wall_mm'], scope='Every emitted inner triangle in a proven ball; every emitted outer triangle separated from that ball by wall requirement. Cervical rim excluded explicitly.')
    dump(R / 'raw' / (key + '.json'), out)
    return dict(all_pass=out['all_pass'], path=str(R / 'raw' / (key + '.json')), sha256=sha(R / 'raw' / (key + '.json')), uncovered_faces=cov['uncovered'], exact_outer_separation=cert)
