from generate import *
import generate
from skimage.measure import marching_cubes

def regularize(v, f, key):
    pr = read(R / 'PREREG_R2.json')['regularization']
    h = pr['grid_mm']
    lo = v.min(0) - pr['padding_mm']
    hi = v.max(0) + pr['padding_mm']
    dims = np.ceil((hi - lo) / h).astype(int) + 1
    N = int(np.prod(dims))
    if N > pr['max_nodes']:
        raise ValueError('GRID_BUDGET ' + str(N))
    ds = Distance(v, f)
    field = np.empty(N, np.float32)
    st = time.perf_counter()
    for start in range(0, N, 32768):
        idx = np.arange(start, min(start + 32768, N))
        pts = lo + np.stack(np.unravel_index(idx, dims), axis=1) * h
        dd = ds.query(pts)[0]
        wn = igl.fast_winding_number(v, f, pts)
        field[idx] = np.where(abs(wn) > 0.5, -dd, dd)
    vol = field.reshape(dims)
    if min(vol[0].min(), vol[-1].min(), vol[:, 0].min(), vol[:, -1].min(), vol[:, :, 0].min(), vol[:, :, -1].min()) <= 0:
        raise ValueError('ISO_TOUCHES_BOX')
    (vv, ff, _, _) = marching_cubes(vol, level=0, spacing=(h, h, h), allow_degenerate=False)
    vv = vv.astype(float) + lo
    m = trimesh.Trimesh(vv, ff, process=False)
    trimesh.repair.fix_normals(m, multibody=True)
    comps = m.split(only_watertight=False)
    volumes = [abs(x.volume) for x in comps]
    j = int(np.argmax(volumes))
    drop = 1 - volumes[j] / sum(volumes)
    m = comps[j]
    if drop > pr['discarded_component_volume_fraction_max']:
        raise ValueError('COMPONENT_DROP ' + str(drop))
    if m.volume < 0:
        m.invert()
    si = intersections(m.vertices, m.faces, key + '_R2_regularized')
    info = dict(method='winding_isosurface', grid_mm=h, grid_shape=dims, nodes=N, components=len(comps), discarded_volume_fraction=drop, seconds=time.perf_counter() - st, self_intersections=si.get('count'), closed=bool(m.is_watertight), original_facets=len(f), final_facets=len(m.faces), source_uncertainty='UNKNOWN')
    dump(R / 'raw' / (key + '_REGULARIZATION.json'), info)
    if si.get('count') != 0 or not m.is_watertight:
        raise ValueError('INVALID_REGULARIZED_SUPPORT ' + str(info))
    path = D / (key + '_R2_support.mesh')
    meshwrite(path, m.vertices, m.faces)
    info.update(path=str(path), sha256=sha(path))
    return (m.vertices.copy(), m.faces.copy(), info)

def supports2(rec):
    key = rec['key']
    tri = load(rec['private_path'])['target']
    (sv, sf) = triangles_mesh(tri)
    try:
        T = close(sv, sf, key + '_native')
    except ValueError:
        (v, f) = meshread(D / (key + '_native_closed.mesh'))
        T = regularize(v, f, key + '_native')
    x = load(rec['outer']['mesh_path'])
    (ev, ef) = (x['vertices'], x['faces'])
    if 'roles' in x:
        (ev, ef) = compact(ev, ef[x['roles'] == 0])
    try:
        E = close(ev, ef, key + '_outer')
    except ValueError:
        (v, f) = meshread(D / (key + '_outer_closed.mesh'))
        E = regularize(v, f, key + '_outer')
    dump(R / 'raw' / (key + '_R2_SUPPORTS.json'), dict(native=T[2], outer=E[2]))
    return (T, E, tri)
if __name__ == '__main__':
    generate.supports = supports2
    generate.run('R2', 3 if '--pilot' in sys.argv else None)
