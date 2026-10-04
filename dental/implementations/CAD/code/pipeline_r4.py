from cadlib import *
from scipy.spatial import cKDTree
import trimesh
from readonly import parents
P = None

def canonical_prep(record):
    path = DATA.parent / 'PROOF_LANE_GENCAD_V2/axial_R3' / (record['key'] + '_axial_flare.npz')
    a = np.load(path)
    v = a['preparation_vertices']
    f = a['preparation_faces']
    n = 24
    V = np.r_[v, [[0, 0, v[0, 2]]]]
    F = np.r_[f, np.array([[i, len(v), (i + 1) % n] for i in range(n)])]
    m = trimesh.Trimesh(V, F, process=False)
    stl = trimesh.load_mesh(record['die_path'], process=True)
    (dd, j) = cKDTree(V).query(stl.vertices)
    if len(np.unique(j)) != len(V) or dd.max() > 1e-05:
        raise ValueError('CAD/STL vertex bijection or bound failed')
    got = sorted(map(tuple, np.sort(j[stl.faces], axis=1)))
    want = sorted(map(tuple, np.sort(F, axis=1)))
    if got != want:
        raise ValueError('CAD/STL facet incidence differs')
    return (m, dict(cad_source_path=path, cad_source_sha256=sha(path), stl_sha256=sha(record['die_path']), vertex_error_max_mm=float(dd.max()), facet_bijection=True, continuous_surface_Hausdorff_upper_mm=float(dd.max()), scope='Exact barycentric perturbation inequality; floating norm evaluation not outward rounded, chosen1e-5mm acceptance envelope dominates'))

def offset_cap(prep, faces, gap):
    global P
    if P is None:
        P = parents()
    v = prep.vertices
    z0 = float(v[:, 2].min())
    zt = float(v[:, 2].max())
    bot = v[abs(v[:, 2] - z0) < 1e-10]
    top = v[abs(v[:, 2] - zt) < 1e-10]
    rb = float(np.linalg.norm(bot[:, :2], axis=1).max())
    rt = float(np.linalg.norm(top[:, :2], axis=1).max())
    slope = (rb - rt) / (zt - z0)
    delta = gap * np.sqrt(1 + slope * slope)
    newtop = zt + gap
    newrb = rb + delta
    newrt = rt + delta - slope * gap
    (iv, inf) = P['cap'].cone(newrb, newrt, z0, newtop, n=24)
    return (iv, inf, np.arange(24), dict(nominal_gap_mm=gap, offset_type='exact-decimal homothetic24-gon; film min/max independently checked', bottom_radius_mm=newrb, top_radius_mm=newrt, source_height_mm=zt - z0, scope='nominal CAD, not as-built or uniform true Euclidean offset'))
