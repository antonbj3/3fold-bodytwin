from common import *
from loop_crown import clip, loops, annulus
from evaluate import boundary_edges, mf
from scipy.spatial import cKDTree

def cap_for(m):
    if 'prep_cap_vertices' in m:
        return (m['prep_cap_vertices'], m['prep_cap_faces'])
    M = m['M']
    b = float(m['base'])
    (qv, qf) = (m['prep_vertices'], m['prep_faces'])
    active = qf[abs((qv[qf] @ M.T)[:, :, 2] - b).max(1) > 1e-07]
    ll = loops(qv, active, M)
    assert len(ll) == 1
    c = m['native_section']
    v = np.r_[c, qv]
    outer = np.arange(len(c))
    inner = ll[0] + len(c)
    (ring, _) = annulus(v, outer, inner, M)
    f = np.r_[ring, active + len(c)]
    mesh = trimesh.Trimesh(v, f, process=False)
    trimesh.repair.fix_normals(mesh)
    if np.median(mesh.face_normals @ M[2]) < 0:
        mesh.faces = mesh.faces[:, ::-1]
    return (mesh.vertices.copy(), mesh.faces.copy())

def one(row):
    rec = next((x for x in inputs() if x['key'] == row['key']))
    m = load(row['mesh_path'])
    M = m['M']
    b = float(m['base'])
    (pv, pf) = cap_for(m)
    pl = loops(pv, pf, M)
    assert len(pl) == 1
    rim = pl[0]
    (tv, tf, info) = native(rec)
    if info.get('method') == 'winding_isosurface':
        return dict(status='UNKNOWN_SOURCE_REGULARIZED', reason='Actual source finish-line and regularized lower model not silently stitched')
    (lv, lf) = clip(tv, tf, b, M, above=False)
    le = boundary_edges(lf)
    lbd = np.unique(le)
    tree = cKDTree(pv[rim])
    (dd, jj) = tree.query(lv[lbd])
    assert dd.max() < 1e-07, ('UNMATCHED_SOURCE_RIM', dd.max())
    mapping = np.arange(len(lv)) + len(pv)
    mapping[lbd] = rim[jj]
    vv = np.r_[pv, lv]
    fout = []
    bound = {tuple(sorted(e)) for e in le}
    xy = vv @ M.T
    rx = xy[rim, :2]
    for face in lf:
        found = False
        for k in range(3):
            (a, bb) = (int(face[k]), int(face[(k + 1) % 3]))
            opp = int(face[(k + 2) % 3])
            if tuple(sorted((a, bb))) not in bound:
                continue
            aa = xy[mapping[a], :2]
            bbb = xy[mapping[bb], :2]
            vec = bbb - aa
            norm = np.dot(vec, vec)
            u = (rx - aa) @ vec / norm
            perp = np.linalg.norm(rx - (aa + u[:, None] * vec), axis=1)
            ids = np.flatnonzero((perp < 1e-08) & (u > 1e-09) & (u < 1 - 1e-09))
            ids = ids[np.argsort(u[ids])]
            seq = np.r_[mapping[a], rim[ids], mapping[bb]]
            for (x, y) in zip(seq[:-1], seq[1:]):
                if x != y:
                    fout.append([mapping[opp], x, y])
            found = True
            break
        if not found:
            fout.append(mapping[face])
    (v, f) = compact(vv, np.r_[pf, np.array(fout)])
    mm = trimesh.Trimesh(v, f, process=False)
    trimesh.repair.fix_normals(mm, multibody=True)
    if mm.volume < 0:
        mm.invert()
    stem = Path(row['mesh_path']).stem
    si = intersections(mm.vertices, mm.faces, stem + '_FULL_PREP')
    path = D / (stem + '_FULL_PREPARATION_RESEARCH.stl')
    mm.export(path)
    stl = trimesh.load_mesh(path, process=False)
    rounderr = float(np.linalg.norm(stl.triangles - mm.triangles, axis=2).max())
    model = D / (stem + '_FULL_PREPARATION.npz')
    np.savez_compressed(model, vertices=mm.vertices, faces=mm.faces)
    out = dict(status='FROZEN_RESEARCH_APPROXIMATION', path=str(path), sha256=sha(path), mesh_path=str(model), mesh_sha256=sha(model), watertight=bool(mm.is_watertight), self_intersections=si.get('count'), retained_preparation_volume_mm3=float(mm.volume), native_closed_model_volume_mm3=float(trimesh.Trimesh(tv, tf, process=False).volume), source_rim_max_weld_mm=float(dd.max()), STL_max_rounding_mm=rounderr, resolution='PER_TOOTH', exact_subset_scope='Intended P=(T below declared plane) union Q is exact by set identity and stored Q certificate. This emitted P STL is a separately rounded/welded research approximation; no literal exact containment claim for rounded P.', physical_status='No CEJ/pulp validation, no clinical or production release')
    dest = R / 'exports' / (stem + '_FULL_PREPARATION_CERTIFICATE.json')
    dump(dest, out)
    return out

def run(tag):
    rows = []
    for r in read(R / f'RESULTS_{tag}.json')['rows']:
        if r['status'] != 'GENERATED':
            continue
        try:
            x = one(r)
        except Exception as e:
            x = dict(status='EXPORT_FAILED', error=repr(e))
            __import__('traceback').print_exc()
        rows.append(dict(key=r['key'], material=r['material'], round=tag, preparation=x))
        dump(R / f'RESULTS_PREPARATION_EXPORT_{tag}.json', dict(rows=rows))
        print(r['key'], r['material'], x.get('status'), x.get('retained_preparation_volume_mm3'), x.get('error'), flush=True)
if __name__ == '__main__':
    run(sys.argv[1])
