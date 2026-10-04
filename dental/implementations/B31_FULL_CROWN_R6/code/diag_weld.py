from local import *
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

def weld(v, f, eps=1e-07):
    pairs = cKDTree(v).query_pairs(eps, output_type='ndarray')
    parent = np.arange(len(v))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for (a, b) in pairs:
        (a, b) = (root(a), root(b))
        parent[max(a, b)] = min(a, b)
    ids = np.array([root(i) for i in range(len(v))])
    ff = ids[f]
    ok = np.array([len(set(face)) == 3 for face in ff])
    ar = np.linalg.norm(np.cross(v[f][:, 1] - v[f][:, 0], v[f][:, 2] - v[f][:, 0]), axis=1) / 2
    (_, unique) = np.unique(np.sort(ff[ok], axis=1), axis=0, return_index=True)
    ff = ff[ok][np.sort(unique)]
    (w, g) = compact(v, ff)
    return (w, g, dict(merged_vertices=len(v) - len(np.unique(ids)), removed_faces=len(f) - len(g), max_vertex_displacement_mm=float(np.linalg.norm(v - v[ids], axis=1).max()), removed_area_sum_mm2=float(ar[~ok].sum()), max_collapsed_face_area_mm2=float(ar[~ok].max()) if (~ok).any() else 0))

def run():
    d = npz(D6 / '079905ebf9504544_molar_C_inner_raw.npz')
    (v, f) = (d['vertices'], d['faces'])
    (w, g, info) = weld(v, f)
    m = trimesh.Trimesh(w, g, process=False)
    info.update(loops=[len(q) for q in loops(m)], components=len(m.split(only_watertight=False)), raw_vertices=len(v), raw_faces=len(f), vertices=len(w), faces=len(g))
    save(R6 / 'raw/WELD_E.json', info)
    np.savez_compressed(D6 / '079905ebf9504544_molar_E_inner.npz', vertices=w, faces=g)
    print(clean(info))
if __name__ == '__main__':
    run()
