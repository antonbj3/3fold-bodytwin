"""Unlabelled PL-mesh feature cycles. Functions never receive a reference line."""
import numpy as np
import trimesh
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

def cycles(edges, vertices):
    edges = np.asarray(edges, int).reshape(-1, 2)
    adj = {}
    for (a, b) in edges:
        adj.setdefault(int(a), []).append(int(b))
        adj.setdefault(int(b), []).append(int(a))
    unseen = set(adj)
    out = []
    rejected = []
    while unseen:
        todo = [min(unseen)]
        cc = set()
        while todo:
            a = todo.pop()
            if a in cc:
                continue
            cc.add(a)
            todo.extend((b for b in adj[a] if b not in cc))
        unseen -= cc
        if len(cc) < 8 or any((len(adj[a]) != 2 for a in cc)):
            rejected.append({'vertices': len(cc), 'reason': 'SMALL_OR_NOT_DEGREE_TWO'})
            continue
        start = min(cc)
        prev = -1
        cur = start
        order = []
        while cur not in order:
            order.append(cur)
            nxt = next((j for j in adj[cur] if j != prev))
            (prev, cur) = (cur, nxt)
        if cur != start or len(order) != len(cc):
            raise ValueError('cycle traversal')
        order = np.asarray(order, int)
        q = np.asarray(vertices)[order]
        out.append(dict(indices=order, points=q, perimeter_mm=float(np.linalg.norm(np.roll(q, -1, axis=0) - q, axis=1).sum()), certificate={'degree_every_vertex': 2, 'components': 1, 'vertices': len(order), 'edges': len(order), 'repeat_count': 0, 'scope': 'Combinatorial simple cycle on supplied PL surface; physical surface embedding unverified'}))
    return (out, rejected)

def feature_cycles(mesh, angle_deg=35.0):
    adj = mesh.face_adjacency
    angles = mesh.face_adjacency_angles
    edges = mesh.face_adjacency_edges[angles >= np.radians(angle_deg)]
    (cs, rej) = cycles(edges, mesh.vertices)
    cs.sort(key=lambda c: -c['perimeter_mm'])
    return (cs, dict(feature_edges=len(edges), closed_components=len(cs), rejected_components=rej))

def detect(mesh):
    (cs, info) = feature_cycles(mesh)
    if not cs:
        return dict(status='UNKNOWN_NO_CLOSED_FEATURE', diagnostics=info)
    return dict(status='PROPOSED', **cs[0], diagnostics=info)

def line_distance(points, curve):
    a = np.asarray(curve)
    d = np.roll(a, -1, axis=0) - a
    dd = np.einsum('ij,ij->i', d, d)
    res = []
    for p in points:
        s = np.clip(np.einsum('ij,ij->i', p - a, d) / np.maximum(dd, 1e-30), 0, 1)
        res.append(np.linalg.norm(a + s[:, None] * d - p, axis=1).min())
    return np.asarray(res)

def dense(q):
    return np.concatenate([q + t * (np.roll(q, -1, axis=0) - q) for t in np.arange(8) / 8])

def errors(a, b):
    d = np.r_[line_distance(dense(a), b), line_distance(dense(b), a)]
    return dict(max_um=float(d.max() * 1000), p95_um=float(np.quantile(d, 0.95) * 1000), mean_um=float(d.mean() * 1000), rms_um=float(np.sqrt(np.mean(d * d)) * 1000), resolution='PER_POINT', sampling='8 equally spaced points per segment per direction; not global Hausdorff certificate', sampling_Lipschitz_remainder_um=float(max(np.linalg.norm(np.roll(a, -1, axis=0) - a, axis=1).max(), np.linalg.norm(np.roll(b, -1, axis=0) - b, axis=1).max()) / 16 * 1000))

def scalar_contour(mesh, values):
    seg = []
    for ids in mesh.faces:
        val = values[ids]
        vv = mesh.vertices[ids]
        pts = []
        for (i, j) in [(0, 1), (1, 2), (2, 0)]:
            if (val[i] < 0) != (val[j] < 0):
                t = val[i] / (val[i] - val[j])
                pts.append(vv[i] + t * (vv[j] - vv[i]))
        if len(pts) == 2:
            seg.append(pts)
    if not seg:
        return []
    p = np.asarray(seg).reshape(-1, 3)
    (unique, inv) = np.unique(np.round(p, 8), axis=0, return_inverse=True)
    (cs, _) = cycles(inv.reshape(-1, 2), unique)
    return sorted(cs, key=lambda c: -c['perimeter_mm'])
