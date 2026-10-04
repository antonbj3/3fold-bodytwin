"""Conservative cubical insertion closure of original labelled triangles."""
from common import *
import trimesh

def clip(poly, axis, value, greater):
    if len(poly) == 0:
        return poly
    out = []
    prev = poly[-1]
    pv = prev[axis] - value
    pin = pv >= -1e-12 if greater else pv <= 1e-12
    for cur in poly:
        cv = cur[axis] - value
        cin = cv >= -1e-12 if greater else cv <= 1e-12
        if cin != pin:
            t = pv / (pv - cv)
            out.append(prev + t * (cur - prev))
        if cin:
            out.append(cur)
        (prev, pv, pin) = (cur, cv, cin)
    return np.array(out).reshape(-1, 3)

def coords(v, d):
    helper = np.eye(3)[np.argmin(np.abs(d))]
    e1 = np.cross(d, helper)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(d, e1)
    basis = np.stack([e1, e2, d], axis=1)
    return (v @ basis, basis)

def raster(tri, h):
    lo = tri.reshape(-1, 3).min(axis=0)
    offset = lo - np.array([h, h, 0.5])
    q = tri - offset
    extent = q.reshape(-1, 3).max(axis=0)
    shape = tuple((np.ceil(extent[:2] / h).astype(int) + 2).tolist())
    height = np.zeros(shape, dtype=np.int32)
    empty = 0
    clips = 0
    max_cells = 0
    for face in q:
        first = np.floor(face[:, :2].min(axis=0) / h).astype(int)
        last = np.floor(face[:, :2].max(axis=0) / h).astype(int)
        seen = 0
        for i in range(first[0], last[0] + 1):
            px = clip(clip(face, 0, i * h, True), 0, (i + 1) * h, False)
            if not len(px):
                continue
            for j in range(first[1], last[1] + 1):
                poly = clip(clip(px, 1, j * h, True), 1, (j + 1) * h, False)
                if not len(poly):
                    continue
                seen += 1
                clips += 1
                roof = int(np.ceil(poly[:, 2].max() / h))
                height[i, j] = max(height[i, j], roof)
        max_cells = max(max_cells, seen)
        if seen == 0:
            empty += 1
    pts = q.reshape(-1, 3)
    ix = np.floor(pts[:, :2] / h).astype(int)
    containment_margin = height[ix[:, 0], ix[:, 1]] * h - pts[:, 2]
    return (height, offset, {'source_triangles': len(tri), 'clipped_polygons': clips, 'uncovered_triangles': empty, 'max_cells_per_triangle': max_cells, 'vertex_containment_min': float(containment_margin.min()), 'cover_gate': bool(empty == 0 and containment_margin.min() >= -1e-10)})

def mesh_columns(height, h, offset, basis):
    quads = []
    xy = np.argwhere(height > 0)
    i = xy[:, 0]
    j = xy[:, 1]
    z = height[i, j]

    def horizontal(k):
        return np.stack([np.c_[i, j, k], np.c_[i + 1, j, k], np.c_[i + 1, j + 1, k], np.c_[i, j + 1, k]], axis=1)
    quads.extend([horizontal(z), horizontal(np.zeros_like(z))[:, ::-1]])
    for axis in [0, 1]:
        pad = np.pad(height, ((1, 1), (0, 0))) if axis == 0 else np.pad(height, ((0, 0), (1, 1)))
        left = pad[:-1, :] if axis == 0 else pad[:, :-1]
        right = pad[1:, :] if axis == 0 else pad[:, 1:]
        loc = np.argwhere(left != right)
        counts = np.abs(left[tuple(loc.T)] - right[tuple(loc.T)])
        if not len(loc):
            continue
        expanded = np.repeat(loc, counts, axis=0)
        starts = np.repeat(np.minimum(left[tuple(loc.T)], right[tuple(loc.T)]), counts)
        indices = np.arange(int(counts.sum())) - np.repeat(np.r_[0, np.cumsum(counts)[:-1]], counts)
        zz = starts + indices
        if axis == 0:
            x = expanded[:, 0]
            y = expanded[:, 1]
            q = np.stack([np.c_[x, y, zz], np.c_[x, y + 1, zz], np.c_[x, y + 1, zz + 1], np.c_[x, y, zz + 1]], axis=1)
        else:
            x = expanded[:, 0]
            y = expanded[:, 1]
            q = np.stack([np.c_[x, y, zz], np.c_[x, y, zz + 1], np.c_[x + 1, y, zz + 1], np.c_[x + 1, y, zz]], axis=1)
        flip = left[tuple(expanded.T)] < right[tuple(expanded.T)]
        q[flip] = q[flip, ::-1]
        quads.append(q)
    quad = np.concatenate(quads)
    (lattice, inv) = np.unique(quad.reshape(-1, 3), axis=0, return_inverse=True)
    fq = inv.reshape(-1, 4)
    faces = np.concatenate([fq[:, [0, 1, 2]], fq[:, [0, 2, 3]]])
    vertices = (lattice * h + offset) @ basis.T
    return (trimesh.Trimesh(vertices, faces, process=False), {'quads': len(quad), 'vertices': len(vertices), 'faces': len(faces)})

def distance(mesh, p):
    return trimesh.proximity.closest_point(mesh, p)[1]

def continuous_gap(mesh, tri, tolerance=0.12, maxdepth=8):
    active = tri.copy()
    upper_max = 0.0
    lower_max = 0.0
    witness = None
    iterations = []
    unresolved = 0
    for depth in range(maxdepth + 1):
        cent = active.mean(axis=1)
        r = np.linalg.norm(active - cent[:, None, :], axis=2).max(axis=1)
        dist = distance(mesh, cent)
        low = float(dist.max())
        im = int(np.argmax(dist))
        if low > lower_max:
            lower_max = low
            witness = cent[im].tolist()
        up = dist + r
        clear = up <= tolerance
        if np.any(clear):
            upper_max = max(upper_max, float(up[clear].max()))
        iterations.append({'depth': depth, 'triangles': len(active), 'sampled_max': low, 'upper_max_active': float(up.max()), 'cleared': int(clear.sum())})
        if np.any(dist > tolerance):
            return {'gate': 'FAIL', 'lower_max': lower_max, 'upper_max': max(upper_max, float(up.max())), 'witness': witness, 'iterations': iterations, 'unresolved': int((~clear).sum())}
        bad = active[~clear]
        if len(bad) == 0:
            return {'gate': 'CONDITIONAL_PASS', 'lower_max': lower_max, 'upper_max': upper_max, 'witness': witness, 'iterations': iterations, 'unresolved': 0}
        if depth == maxdepth:
            return {'gate': 'UNKNOWN', 'lower_max': lower_max, 'upper_max': max(upper_max, float(up.max())), 'witness': witness, 'iterations': iterations, 'unresolved': len(bad)}
        (a, b, c) = (bad[:, 0], bad[:, 1], bad[:, 2])
        ab = (a + b) / 2
        bc = (b + c) / 2
        ca = (c + a) / 2
        active = np.concatenate([np.stack([a, ab, ca], axis=1), np.stack([ab, b, bc], axis=1), np.stack([ca, bc, c], axis=1), np.stack([ab, bc, ca], axis=1)])
        if len(active) > 2000000:
            raise MemoryError('frozen lane adaptive array limit reached; do not launch a >4 GB job')
