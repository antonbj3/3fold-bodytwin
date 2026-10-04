"""Triangle operators adapted from DIAG/local_pairs; all distances in mm."""
from common import *
import trimesh
from scipy.spatial import cKDTree
from shapely.geometry import Polygon, Point

def read_obj(p):
    vv = []
    ff = []
    for line in Path(p).open():
        if line.startswith('v '):
            vv.append([float(x) for x in line.split()[1:4]])
        elif line.startswith('f '):
            q = [int(x.split('/')[0]) - 1 for x in line.split()[1:]]
            if min(q) < 0:
                raise ValueError('negative OBJ indices')
            ff.extend([[q[0], q[j], q[j + 1]] for j in range(1, len(q) - 1)])
    return (np.array(vv), np.array(ff, int))

def mesh(t, process=False):
    return trimesh.Trimesh(t.reshape(-1, 3), np.arange(len(t) * 3).reshape(-1, 3), process=process)

def component(t):
    m = mesh(t, True)
    m.update_faces(m.unique_faces())
    m.update_faces(m.nondegenerate_faces())
    m.remove_unreferenced_vertices()
    pp = m.split(only_watertight=False)
    best = max(pp, key=lambda q: q.area)
    return (best, float(1 - best.area / m.area), len(pp))

def sample(t, n):
    a = np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1) / 2
    ids = np.searchsorted(np.cumsum(a), (np.arange(n) + 0.5) * a.sum() / n)
    return t[ids].mean(1)

def closest(t, p):
    m = mesh(t)
    rr = [trimesh.proximity.closest_point(m, q) for q in np.array_split(p, max(1, int(np.ceil(len(p) / 64))))]
    return tuple((np.concatenate([r[k] for r in rr]) for k in range(3)))

def distance(t, p):
    return closest(t, p)[1]

def fit(s, t):
    cs = s.mean(0)
    ct = t.mean(0)
    ss = s - cs
    tt = t - ct
    tree = cKDTree(tt)
    out = []
    for angle in [0, np.pi / 2, np.pi, 3 * np.pi / 2]:
        c = np.cos(angle)
        sn = np.sin(angle)
        R = np.array([[c, -sn, 0], [sn, c, 0], [0, 0, 1.0]])
        u = np.zeros(3)
        for _ in range(60):
            a = ss @ R.T + u
            (d, j) = tree.query(a, workers=1)
            (dr, ir) = cKDTree(a).query(tt, workers=1)
            x = np.r_[ss, ss[ir]]
            y = np.r_[tt[j], tt]
            cx = x.mean(0)
            cy = y.mean(0)
            (U, _, V) = np.linalg.svd((x - cx).T @ (y - cy))
            Rn = V.T @ U.T
            if np.linalg.det(Rn) < 0:
                V[-1] *= -1
                Rn = V.T @ U.T
            R = Rn
            u = cy - cx @ R.T
        a = ss @ R.T + u
        err = float(np.sqrt(np.mean(np.r_[tree.query(a)[0], cKDTree(a).query(tt)[0]] ** 2)))
        out.append((err, R, ct + u - cs @ R.T))
    best = min(out, key=lambda x: x[0])
    return (best[1], best[2], [x[0] for x in out])

def metrics(a, b, n=8192):
    da = distance(b, sample(a, n))
    db = distance(a, sample(b, n))
    return dict(p95_mm=float(max(np.quantile(da, 0.95), np.quantile(db, 0.95))), rms_mm=float(np.sqrt(np.mean(np.r_[da, db] ** 2))), sampled_max_mm=float(max(da.max(), db.max())), directed_p95_mm=[float(np.quantile(d, 0.95)) for d in [da, db]])

def loops(m):
    edges = np.concatenate([m.faces[:, [0, 1]], m.faces[:, [1, 2]], m.faces[:, [2, 0]]])
    (_, ix, ct) = np.unique(np.sort(edges, axis=1), axis=0, return_index=True, return_counts=True)
    if np.any(ct > 2):
        raise ValueError('nonmanifold edge')
    ee = edges[ix[ct == 1]]
    adj = {}
    for (a, b) in ee:
        adj.setdefault(int(a), []).append(int(b))
        adj.setdefault(int(b), []).append(int(a))
    if any((len(x) != 2 for x in adj.values())):
        raise ValueError('branching boundary')
    unseen = set(adj)
    out = []
    while unseen:
        a = min(unseen)
        loop = [a]
        prev = -1
        cur = a
        while True:
            nex = next((j for j in adj[cur] if j != prev))
            if nex == a:
                break
            if nex in loop:
                raise ValueError('boundary repeats')
            loop.append(nex)
            (prev, cur) = (cur, nex)
        unseen.difference_update(loop)
        out.append(np.array(loop, int))
    return out

def cut(m, z):
    out = m.slice_plane([0, 0, z], [0, 0, 1], cap=False)
    out.merge_vertices(digits_vertex=10)
    return out

def boundary(m, z):
    ll = loops(m)
    basal = [l for l in ll if np.max(abs(m.vertices[l, 2] - z)) < 1e-07]
    if len(basal) != 1 or len(ll) != 1:
        raise ValueError('need one planar loop: total=%d basal=%d' % (len(ll), len(basal)))
    return basal[0]

def nearest_curve(v, curve):
    a = curve
    b = np.roll(curve, -1, axis=0)
    d = b - a
    out = []
    for p in v:
        w = np.clip(np.einsum('ij,ij->i', p - a, d) / np.maximum(np.einsum('ij,ij->i', d, d), 1e-30), 0, 1)
        q = a + w[:, None] * d
        out.append(q[np.argmin(np.linalg.norm(q - p, axis=1))])
    return np.array(out)

def curve_error(a, b):

    def dense(q):
        return np.concatenate([q + t * (np.roll(q, -1, axis=0) - q) for t in np.linspace(0, 1, 8, endpoint=False)])
    aa = dense(a)
    bb = dense(b)
    return float(max(np.linalg.norm(aa - nearest_curve(aa, b), axis=1).max(), np.linalg.norm(bb - nearest_curve(bb, a), axis=1).max()))

def grid(t, h=0.25):
    v = t.reshape(-1, 3)
    axes = [np.arange(v[:, k].min(), v[:, k].max() + h / 2, h) for k in range(2)]
    (xx, yy) = np.meshgrid(*axes, indexing='ij')
    i = np.arange(xx.size).reshape(xx.shape)
    f = np.r_[np.stack([i[:-1, :-1], i[1:, :-1], i[:-1, 1:]], -1).reshape(-1, 3), np.stack([i[1:, :-1], i[1:, 1:], i[:-1, 1:]], -1).reshape(-1, 3)]
    return (np.c_[xx.ravel(), yy.ravel()], f)
