"""Small explicit geometric operators. No reconstruction is called measured truth."""
import numpy as np
from scipy.spatial import cKDTree

def unit(v):
    v = np.asarray(v, float)
    n = np.linalg.norm(v)
    if n < 1e-09:
        raise ValueError('zero direction')
    return v / n

def parse_stl(blob):
    if len(blob) < 84:
        raise ValueError('empty/truncated STL')
    n = int.from_bytes(blob[80:84], 'little')
    if len(blob) == 84 + 50 * n:
        dt = np.dtype([('n', '<f4', (3,)), ('v', '<f4', (3, 3)), ('a', '<u2')])
        tri = np.frombuffer(blob, dtype=dt, count=n, offset=84)['v'].astype(float)
        (v, ix) = np.unique(tri.reshape(-1, 3), axis=0, return_inverse=True)
        return (v, ix.reshape(-1, 3))
    import trimesh, io
    m = trimesh.load_mesh(io.BytesIO(blob), file_type='stl', process=False)
    (v, ix) = np.unique(np.asarray(m.vertices)[m.faces].reshape(-1, 3), axis=0, return_inverse=True)
    return (v, ix.reshape(-1, 3))

def crop(tri, xy, pad=0.0):
    lo = xy.min(0) - pad
    hi = xy.max(0) + pad
    good = np.all(tri[:, :, :2].max(1) >= lo, axis=1) & np.all(tri[:, :, :2].min(1) <= hi, axis=1)
    t = tri[good]
    det = np.cross(t[:, 1, :2] - t[:, 0, :2], t[:, 2, :2] - t[:, 0, :2])
    return t[np.abs(det) > 1e-09]

def height(tri, xy, upper=False):
    """All triangles whose XY projection contains each query; exact affine evaluation in float."""
    out = np.full(len(xy), np.nan)
    if not len(tri):
        return out
    c = tri[:, :, :2].mean(1)
    r = np.linalg.norm(tri[:, :, :2] - c[:, None, :], axis=2).max(1)
    hits = cKDTree(c).query_ball_point(xy, float(r.max()), workers=1)
    for (i, js) in enumerate(hits):
        if not js:
            continue
        t = tri[js]
        (a, b, c) = (t[:, 0, :2], t[:, 1, :2], t[:, 2, :2])
        p = xy[i]
        den = np.cross(b - a, c - a)
        ok = np.abs(den) > 1e-12
        w1 = np.divide(np.cross(b - p, c - p), den, out=np.zeros(len(t)), where=ok)
        w2 = np.divide(np.cross(c - p, a - p), den, out=np.zeros(len(t)), where=ok)
        w3 = 1 - w1 - w2
        ok &= (w1 >= -1e-09) & (w2 >= -1e-09) & (w3 >= -1e-09)
        if ok.any():
            z = w1 * t[:, 0, 2] + w2 * t[:, 1, 2] + w3 * t[:, 2, 2]
            out[i] = z[ok].min() if upper else z[ok].max()
    return out

def grid(half, family, n=13):
    ax = np.linspace(-1, 1, n)
    (u, v) = np.meshgrid(ax, ax, indexing='ij')
    uv = np.c_[u.ravel(), v.ravel()]
    if family == 'bridge3':
        keep = (np.minimum.reduce([(uv[:, 0] - c) ** 2 / 0.42 ** 2 + uv[:, 1] ** 2 for c in [-0.62, 0, 0.62]]) <= 1) | (np.abs(uv[:, 1]) <= 0.38)
    else:
        keep = uv[:, 0] ** 2 + uv[:, 1] ** 2 <= 0.96 ** 2
    if family == 'implant_crown':
        keep &= np.sum(uv ** 2, axis=1) >= 0.22 ** 2
    if family == 'lattice_onlay':
        for x in [-0.42, 0.42]:
            for y in [-0.42, 0.42]:
                keep &= np.sum((uv - [x, y]) ** 2, axis=1) >= 0.2 ** 2
    idx = np.full(n * n, -1, int)
    idx[keep] = np.arange(keep.sum())
    idx = idx.reshape(n, n)
    f = []
    for i in range(n - 1):
        for j in range(n - 1):
            for ij in [[(i, j), (i + 1, j), (i, j + 1)], [(i + 1, j), (i + 1, j + 1), (i, j + 1)]]:
                ids = [idx[a, b] for (a, b) in ij]
                if min(ids) >= 0:
                    f.append(ids)
    f = np.asarray(f, int)
    used = np.unique(f)
    mapping = np.full(keep.sum(), -1, int)
    mapping[used] = np.arange(len(used))
    return (uv[keep][used] * np.asarray(half), mapping[f], uv[keep][used])

def boundary_edges(faces):
    edges = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    (_, ix, ct) = np.unique(np.sort(edges, axis=1), axis=0, return_index=True, return_counts=True)
    return edges[ix[ct == 1]]

def shell(xy, outer, inner, faces):
    n = len(xy)
    v = np.r_[np.c_[xy, outer], np.c_[xy, inner]]
    f = np.r_[faces, faces[:, ::-1] + n]
    side = []
    for (a, b) in boundary_edges(faces):
        side.extend([[a, b + n, b], [a, a + n, b + n]])
    return (v, np.r_[f, np.asarray(side, int)])

def area_weights(xy, f):
    a = np.abs(np.cross(xy[f[:, 1]] - xy[f[:, 0]], xy[f[:, 2]] - xy[f[:, 0]])) / 2
    w = np.zeros(len(xy))
    np.add.at(w, f.ravel(), np.repeat(a / 3, 3))
    return w

def section_weights(xy, faces, x):
    """Integral weights for a vertical y-z section at x; exact piecewise-linear thickness integral."""
    w = np.zeros(len(xy))
    x = float(x) + 1e-08
    for ids in faces:
        t = xy[ids]
        points = []
        barys = []
        for (i, j) in [(0, 1), (1, 2), (2, 0)]:
            if t[i, 0] < x < t[j, 0] or t[j, 0] < x < t[i, 0]:
                a = (x - t[i, 0]) / (t[j, 0] - t[i, 0])
                p = (1 - a) * t[i] + a * t[j]
                b = np.zeros(3)
                b[i] = 1 - a
                b[j] = a
                points.append(p)
                barys.append(b)
        if len(points) == 2:
            length = abs(points[1][1] - points[0][1])
            w[ids] += length * (barys[0] + barys[1]) / 2
    return w

def export_stl(p, v, f):
    from pathlib import Path
    t = v[f]
    dt = np.dtype([('normal', '<f4', (3,)), ('v', '<f4', (3, 3)), ('a', '<u2')])
    a = np.zeros(len(f), dt)
    a['v'] = t
    norm = np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])
    a['normal'] = norm / np.maximum(np.linalg.norm(norm, axis=1)[:, None], 1e-20)
    Path(p).write_bytes(b'GenCAD v2 research surface shell; millimetres'.ljust(80, b' ') + len(f).to_bytes(4, 'little') + a.tobytes())
