"""Exact rational exposure of a finite union of square frusta (two-core export)."""
from common import *
from carriers import S, H
from fractions import Fraction as Q
import trimesh

def dot(a, b):
    return sum((x * y for (x, y) in zip(a, b)))

def sub(a, b):
    return tuple((x - y for (x, y) in zip(a, b)))

def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])

def unique(poly):
    out = []
    for p in poly:
        if not out or p != out[-1]:
            out.append(p)
    if len(out) > 1 and out[0] == out[-1]:
        out.pop()
    return out

def positive_area(poly):
    return len(poly) >= 3 and any((any((x != 0 for x in cross(sub(poly[j], poly[0]), sub(poly[j + 1], poly[0])))) for j in range(1, len(poly) - 1)))

def clip(poly, n, b, inside=True):
    if not poly:
        return []
    vv = [dot(n, p) - b for p in poly]
    out = []
    for (i, p) in enumerate(poly):
        j = (i + 1) % len(poly)
        q = poly[j]
        a = vv[i]
        c = vv[j]
        kp = a <= 0 if inside else a >= 0
        kq = c <= 0 if inside else c >= 0
        if kp:
            out.append(p)
        if kp != kq:
            f = a / (a - c)
            out.append(tuple((p[k] + f * (q[k] - p[k]) for k in range(3))))
    return unique(out)

def primitive(center, apex, base, r0):
    (x, y) = map(Q, center)
    a = Q(apex)
    z = Q(base)
    r0 = Q(r0)
    rb = r0 + S * (a - z)
    vb = [(x - rb, y - rb, z), (x + rb, y - rb, z), (x + rb, y + rb, z), (x - rb, y + rb, z)]
    vt = [(x - r0, y - r0, a), (x + r0, y - r0, a), (x + r0, y + r0, a), (x - r0, y + r0, a)]
    ns = [(Q(0), Q(-1), S), (Q(1), Q(0), S), (Q(0), Q(1), S), (Q(-1), Q(0), S), (Q(0), Q(0), Q(1)), (Q(0), Q(0), Q(-1))]
    faces = [[vb[j], vb[(j + 1) % 4], vt[(j + 1) % 4], vt[j]] for j in range(4)] + [vt, vb[::-1]]
    return {'faces': faces, 'planes': [(n, dot(n, f[0])) for (n, f) in zip(ns, faces)]}

def exposed(prims):
    polys = []
    tags = []
    for (i, P) in enumerate(prims):
        for (fi, poly) in enumerate(P['faces']):
            pieces = [poly]
            for (j, R) in enumerate(prims):
                if i == j:
                    continue
                (n, b) = P['planes'][fi]
                if i < j and any((n == nn and b == bb for (nn, bb) in R['planes'])):
                    continue
                keep = []
                for piece in pieces:
                    cur = piece
                    for (nn, bb) in R['planes']:
                        vals = [dot(nn, p) - bb for p in cur]
                        if vals and max(vals) > 0:
                            outside = clip(cur, nn, bb, False)
                            if positive_area(outside):
                                keep.append(outside)
                        cur = clip(cur, nn, bb, True)
                        if not positive_area(cur):
                            break
                pieces = keep
            for p in pieces:
                if positive_area(p):
                    polys.append(p)
                    tags.append([i, fi])
    return (polys, tags)

def mesh_exact(prims):
    (polys, tags) = exposed(prims)
    allpts = sorted(set((p for poly in polys for p in poly)))
    seampts = allpts.copy()
    index = {p: i for (i, p) in enumerate(allpts)}
    faces = []
    roles = []
    for (poly, tag) in zip(polys, tags):
        ring = []
        for (i, a) in enumerate(poly):
            b = poly[(i + 1) % len(poly)]
            d = sub(b, a)
            den = dot(d, d)
            points = []
            for p in seampts:
                q = sub(p, a)
                if cross(q, d) == (0, 0, 0) and 0 <= dot(q, d) <= den:
                    points.append((dot(q, d) / den, p))
            ring.extend((p for (t, p) in sorted(points)[:-1]))
        ring = unique(ring)
        center = tuple((sum((p[k] for p in ring)) / len(ring) for k in range(3)))
        ci = len(allpts)
        allpts.append(center)
        for (j, a) in enumerate(ring):
            b = ring[(j + 1) % len(ring)]
            faces.append([ci, index[a], index[b]])
            roles.append(tag[1])
    verts = np.array([[float(q) for q in p] for p in allpts])
    m = trimesh.Trimesh(verts, np.array(faces), process=False)
    volume = Q(0)
    margins = []
    for (f, role) in zip(faces, roles):
        (a, b, c) = [allpts[i] for i in f]
        n = cross(sub(b, a), sub(c, a))
        volume += dot(a, cross(b, c)) / 6
        if role != 5:
            margins.append(n[2] - S * (abs(n[0]) + abs(n[1])))
    return (m, {'exact_volume_mm3': str(volume), 'volume_mm3': float(volume), 'exact_cone_min_unnormalized_margin': str(min(margins)), 'exact_cone_failures': sum((q < 0 for q in margins)), 'vertices': len(verts), 'faces': len(faces), 'watertight': bool(m.is_watertight), 'winding_consistent': bool(m.is_winding_consistent), 'float_serialization_enclosure': 'MISSING'}, allpts, np.array(roles))

def volume1(a, base, r0):
    L = a - base
    return 4 * (r0 * r0 * L + r0 * S * L * L + S * S * L ** 3 / 3)

def volume2(c1, a1, c2, a2, base, r0):
    r1 = r0 + S * (a1 - base)
    r2 = r0 + S * (a2 - base)
    wx = min(c1[0] + r1, c2[0] + r2) - max(c1[0] - r1, c2[0] - r2)
    wy = min(c1[1] + r1, c2[1] + r2) - max(c1[1] - r1, c2[1] - r2)
    if min(wx, wy) <= 0:
        return None
    U = min(a1 - base, a2 - base, wx / (2 * S), wy / (2 * S))
    overlap = wx * wy * U - S * (wx + wy) * U ** 2 + 4 * S * S * U ** 3 / 3
    return volume1(a1, base, r0) + volume1(a2, base, r0) - overlap

def choose_pair(coords, base):
    r0 = H / 2 - S * H / 2
    s = float(S)
    b = float(base)
    rr = float(r0) + s * (coords[:, 2] - b)
    vol = 4 * (float(r0) ** 2 * (coords[:, 2] - b) + float(r0) * s * (coords[:, 2] - b) ** 2 + s * s * (coords[:, 2] - b) ** 3 / 3)
    candidates = []
    for i in range(len(coords) - 1):
        cs = coords[i + 1:]
        rx = rr[i + 1:]
        wx = np.minimum(coords[i, 0] + rr[i], cs[:, 0] + rx) - np.maximum(coords[i, 0] - rr[i], cs[:, 0] - rx)
        wy = np.minimum(coords[i, 1] + rr[i], cs[:, 1] + rx) - np.maximum(coords[i, 1] - rr[i], cs[:, 1] - rx)
        good = (wx > 0) & (wy > 0)
        U = np.maximum(0, np.minimum(np.minimum(coords[i, 2] - b, cs[:, 2] - b), np.minimum(wx, wy) / (2 * s)))
        ov = wx * wy * U - s * (wx + wy) * U ** 2 + 4 * s * s * U ** 3 / 3
        vv = vol[i] + vol[i + 1:] - ov
        for j in np.flatnonzero(good):
            candidates.append((float(vv[j]), i, i + 1 + int(j)))
    if not candidates:
        raise ValueError('NO_CONNECTED_TWO_CORE_PAIR')
    best = max((v[0] for v in candidates))
    near = [q for q in candidates if q[0] >= best - 1e-10]
    out = []
    for (_, i, j) in near:
        c1 = tuple((Q(float(x)) for x in coords[i, :2]))
        c2 = tuple((Q(float(x)) for x in coords[j, :2]))
        a1 = Q(float(coords[i, 2]))
        a2 = Q(float(coords[j, 2]))
        out.append((volume2(c1, a1, c2, a2, Q(float(base)), r0), i, j))
    (v, i, j) = max(out)
    return ([i, j], {'finite_connected_pairs': len(candidates), 'near_best_exact_checked': len(out), 'selected_pair_volume_mm3': float(v), 'selected_pair_exact_volume_mm3': str(v), 'pair_selection_optimality': 'FLOAT_SCREENED; rigorous all-pair ordering enclosure MISSING', 'maximum_single_volume_mm3': float(vol.max())})
