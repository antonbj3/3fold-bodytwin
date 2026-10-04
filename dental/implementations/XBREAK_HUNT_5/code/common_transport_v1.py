from dental_release.paths import expand as _release_expand
import os
for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    os.environ[key] = '4'
import numpy as np, json, pathlib, hashlib, time
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = pathlib.Path(_release_expand('@DENTAL_WORK_ROOT@/LANE_XBREAK_HUNT_5'))

def write(path, value):
    p = ROOT / path
    p.parent.mkdir(exist_ok=True, parents=True)
    p.write_text(json.dumps(value, indent=2, allow_nan=False))

def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()

def vtp(path):
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    p = pathlib.Path(path)
    original = p.read_bytes()
    fixed = original.replace(b'\\n', b'\n')
    if fixed != original:
        repaired = DATA / (p.stem + '_transport_fixed.vtp')
        repaired.write_bytes(fixed)
        write('raw/TRANSPORT_' + p.stem + '.json', {'source_sha256': sha(p), 'repaired_sha256': sha(repaired), 'literal_newline_replacements': original.count(b'\\n'), 'operation': 'replace literal backslash-n by newline; base64 alphabet has no backslash; no compressed geometry bytes edited', 'original_file': str(p), 'repaired_file': str(repaired)})
        path = repaired
    rd = vtk.vtkXMLPolyDataReader()
    rd.SetFileName(str(path))
    rd.Update()
    m = rd.GetOutput()
    v = vtk_to_numpy(m.GetPoints().GetData()).astype(np.float64)
    conn = vtk_to_numpy(m.GetPolys().GetData()).reshape(-1, 4)
    if np.any(conn[:, 0] != 3):
        raise ValueError('non-triangular input')
    f = conn[:, 1:].astype(np.int64)
    lab = vtk_to_numpy(m.GetCellData().GetArray('Label')).astype(np.int64)
    if not set(np.unique(lab)).issubset({0, 1}):
        raise ValueError('unexpected labels')
    return (v, f, lab)

def frame(v, f, lab):
    tri = v[f]
    raw = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    area = np.linalg.norm(raw, axis=1) / 2
    centroid = tri.mean(axis=1)
    ct = np.average(centroid[lab == 1], axis=0, weights=area[lab == 1])
    cb = np.average(centroid[lab == 0], axis=0, weights=area[lab == 0])
    a = ct - cb
    a /= np.linalg.norm(a)
    rawtop = raw[lab == 1]
    if rawtop.sum(axis=0) @ a < 0:
        rawtop = -rawtop
    n = rawtop / np.maximum(np.linalg.norm(rawtop, axis=1)[:, None], 1e-300)
    t = np.eye(3)[np.argmin(np.abs(a))]
    e1 = np.cross(a, t)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(a, e1)
    return (a, e1, e2, n, area[lab == 1], ct)

def topology(v, f, lab):
    tf = f[lab == 1]
    edges = np.sort(np.concatenate([tf[:, [0, 1]], tf[:, [1, 2]], tf[:, [2, 0]]]), axis=1)
    (ue, inv, cnt) = np.unique(edges, axis=0, return_inverse=True, return_counts=True)
    b = ue[cnt == 1]
    (verts, deg) = np.unique(b, return_counts=True)
    adj = coo_matrix((np.ones(len(ue) * 2), (np.r_[ue[:, 0], ue[:, 1]], np.r_[ue[:, 1], ue[:, 0]])), shape=(len(v), len(v))).tocsr()
    used = np.unique(tf)
    nc = connected_components(adj[used][:, used], directed=False, return_labels=False)
    ba = coo_matrix((np.ones(len(b) * 2), (np.r_[b[:, 0], b[:, 1]], np.r_[b[:, 1], b[:, 0]])), shape=(len(v), len(v))).tocsr()
    loops = int(connected_components(ba[verts][:, verts], directed=False, return_labels=False)) if len(verts) else 0
    be = np.sort(np.concatenate([f[lab == 0][:, [0, 1]], f[lab == 0][:, [1, 2]], f[lab == 0][:, [2, 0]]]), axis=1)
    sb = set(map(tuple, be))
    missing = sum((tuple(e) not in sb for e in b))
    chi = len(used) - len(ue) + len(tf)
    return ({'top_faces': len(tf), 'used_vertices': len(used), 'boundary_edges': len(b), 'boundary_loops': loops, 'non_degree2_vertices': int(np.sum(deg != 2)), 'nonmanifold_edges': int(np.sum(cnt > 2)), 'top_components': int(nc), 'boundary_edges_not_in_label_interface': int(missing), 'euler_characteristic': int(chi), 'disk_gate': bool(nc == 1 and loops == 1 and np.all(deg == 2) and (np.max(cnt) <= 2) and (missing == 0) and (chi == 1))}, b)

def orientloop(edges):
    neigh = {}
    for (x, y) in edges:
        neigh.setdefault(int(x), []).append(int(y))
        neigh.setdefault(int(y), []).append(int(x))
    if any((len(val) != 2 for val in neigh.values())):
        return None
    start = min(neigh)
    seq = [start]
    prev = None
    cur = start
    while True:
        nxt = next((q for q in neigh[cur] if q != prev), None)
        if nxt == start:
            break
        if nxt is None or nxt in seq:
            return None
        seq.append(nxt)
        (prev, cur) = (cur, nxt)
    return np.array(seq) if len(seq) == len(neigh) else None
