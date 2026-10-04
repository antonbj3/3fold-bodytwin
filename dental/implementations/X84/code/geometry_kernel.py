"""Snapshot of X21 affine triangle kernels; no change to other lane."""
import numpy as np
from numba import njit

@njit(cache=True)
def cross2(ax, ay, bx, by):
    return ax * by - ay * bx

@njit(cache=True)
def affine(T):
    n = len(T)
    out = np.empty((n, 3))
    valid = np.ones(n, np.bool_)
    for i in range(n):
        ax = T[i, 1, 0] - T[i, 0, 0]
        ay = T[i, 1, 1] - T[i, 0, 1]
        az = T[i, 1, 2] - T[i, 0, 2]
        bx = T[i, 2, 0] - T[i, 0, 0]
        by = T[i, 2, 1] - T[i, 0, 1]
        bz = T[i, 2, 2] - T[i, 0, 2]
        det = ax * by - ay * bx
        if abs(det) <= 1e-12:
            valid[i] = False
            out[i, :] = 0.0
            continue
        a = (az * by - bz * ay) / det
        b = (ax * bz - bx * az) / det
        out[i, 0] = a
        out[i, 1] = b
        out[i, 2] = T[i, 0, 2] - a * T[i, 0, 0] - b * T[i, 0, 1]
    return (out, valid)

@njit(cache=True)
def intersection_min(u, l, coeff):
    poly = np.empty((12, 2))
    tmp = np.empty((12, 2))
    for k in range(3):
        poly[k, 0] = u[k, 0]
        poly[k, 1] = u[k, 1]
    n = 3
    det = cross2(l[1, 0] - l[0, 0], l[1, 1] - l[0, 1], l[2, 0] - l[0, 0], l[2, 1] - l[0, 1])
    sgn = 1.0 if det > 0 else -1.0
    for e in range(3):
        ex = l[(e + 1) % 3, 0] - l[e, 0]
        ey = l[(e + 1) % 3, 1] - l[e, 1]
        m = 0
        if n == 0:
            return (np.inf, 0.0, 0.0)
        for k in range(n):
            prev = (k - 1) % n
            px = poly[prev, 0]
            py = poly[prev, 1]
            qx = poly[k, 0]
            qy = poly[k, 1]
            dp = sgn * cross2(ex, ey, px - l[e, 0], py - l[e, 1])
            dq = sgn * cross2(ex, ey, qx - l[e, 0], qy - l[e, 1])
            pin = dp >= -1e-12
            qin = dq >= -1e-12
            if pin != qin:
                den = dp - dq
                if abs(den) > 1e-30:
                    a = dp / den
                    tmp[m, 0] = px + a * (qx - px)
                    tmp[m, 1] = py + a * (qy - py)
                    m += 1
            if qin:
                tmp[m, 0] = qx
                tmp[m, 1] = qy
                m += 1
        n = m
        for k in range(n):
            poly[k, 0] = tmp[k, 0]
            poly[k, 1] = tmp[k, 1]
    best = np.inf
    x = 0.0
    y = 0.0
    for k in range(n):
        g = coeff[0] * poly[k, 0] + coeff[1] * poly[k, 1] + coeff[2]
        if g < best:
            best = g
            x = poly[k, 0]
            y = poly[k, 1]
    return (best, x, y)

@njit(cache=True)
def boxes(T):
    out = np.empty((len(T), 6))
    for i in range(len(T)):
        for a in range(3):
            out[i, a] = min(T[i, 0, a], T[i, 1, a], T[i, 2, a])
            out[i, a + 3] = max(T[i, 0, a], T[i, 1, a], T[i, 2, a])
    return out

@njit(cache=True)
def buckets(box, origin, h, shape):
    count = np.zeros(shape[0] * shape[1], np.int64)
    for i in range(len(box)):
        x0 = max(0, int(np.floor((box[i, 0] - origin[0]) / h)))
        x1 = min(shape[0] - 1, int(np.floor((box[i, 3] - origin[0]) / h)))
        y0 = max(0, int(np.floor((box[i, 1] - origin[1]) / h)))
        y1 = min(shape[1] - 1, int(np.floor((box[i, 4] - origin[1]) / h)))
        for x in range(x0, x1 + 1):
            for y in range(y0, y1 + 1):
                count[x * shape[1] + y] += 1
    ptr = np.empty(len(count) + 1, np.int64)
    ptr[0] = 0
    for i in range(len(count)):
        ptr[i + 1] = ptr[i] + count[i]
    ids = np.empty(ptr[-1], np.int32)
    pos = ptr[:-1].copy()
    for i in range(len(box)):
        x0 = max(0, int(np.floor((box[i, 0] - origin[0]) / h)))
        x1 = min(shape[0] - 1, int(np.floor((box[i, 3] - origin[0]) / h)))
        y0 = max(0, int(np.floor((box[i, 1] - origin[1]) / h)))
        y1 = min(shape[1] - 1, int(np.floor((box[i, 4] - origin[1]) / h)))
        for x in range(x0, x1 + 1):
            for y in range(y0, y1 + 1):
                k = x * shape[1] + y
                ids[pos[k]] = i
                pos[k] += 1
    return (ptr, ids)

@njit(cache=True)
def extrema(U, L, ul, ll, uc, lc, ub, lb, origin, h, shape, ptr, ids, prune=True):
    best = np.full((29, 49), np.inf)
    wi = np.full((29, 49, 2), -1, np.int32)
    xy = np.zeros((29, 49, 2))
    seen = np.full(len(L), -1, np.int32)
    checked = 0
    tested = 0
    for i in range(len(U)):
        a = ul[i]
        x0 = max(0, int(np.floor((ub[i, 0] - origin[0]) / h)))
        x1 = min(shape[0] - 1, int(np.floor((ub[i, 3] - origin[0]) / h)))
        y0 = max(0, int(np.floor((ub[i, 1] - origin[1]) / h)))
        y1 = min(shape[1] - 1, int(np.floor((ub[i, 4] - origin[1]) / h)))
        for x in range(x0, x1 + 1):
            for y in range(y0, y1 + 1):
                bkt = x * shape[1] + y
                for k in range(ptr[bkt], ptr[bkt + 1]):
                    j = ids[k]
                    if seen[j] == i:
                        continue
                    seen[j] = i
                    b = ll[j]
                    if ub[i, 0] > lb[j, 3] or ub[i, 3] < lb[j, 0] or ub[i, 1] > lb[j, 4] or (ub[i, 4] < lb[j, 1]):
                        continue
                    checked += 1
                    if prune and ub[i, 2] - lb[j, 5] >= best[a, b]:
                        continue
                    (g, px, py) = intersection_min(U[i], L[j], uc[i] - lc[j])
                    tested += 1
                    if g < best[a, b]:
                        best[a, b] = g
                        wi[a, b, 0] = i
                        wi[a, b, 1] = j
                        xy[a, b, 0] = px
                        xy[a, b, 1] = py
    return (best, wi, xy, checked, tested)

@njit(cache=True)
def raster(T, lab, coeff, origin, h, shape, upper):
    z = np.full(shape[0] * shape[1], np.inf if upper else -np.inf)
    own = np.full(len(z), -1, np.int32)
    bb = boxes(T)
    for i in range(len(T)):
        det = cross2(T[i, 1, 0] - T[i, 0, 0], T[i, 1, 1] - T[i, 0, 1], T[i, 2, 0] - T[i, 0, 0], T[i, 2, 1] - T[i, 0, 1])
        sign = 1.0 if det > 0 else -1.0
        x0 = max(0, int(np.ceil((bb[i, 0] - origin[0]) / h - 0.5)))
        x1 = min(shape[0] - 1, int(np.floor((bb[i, 3] - origin[0]) / h - 0.5)))
        y0 = max(0, int(np.ceil((bb[i, 1] - origin[1]) / h - 0.5)))
        y1 = min(shape[1] - 1, int(np.floor((bb[i, 4] - origin[1]) / h - 0.5)))
        for x in range(x0, x1 + 1):
            for y in range(y0, y1 + 1):
                px = origin[0] + (x + 0.5) * h
                py = origin[1] + (y + 0.5) * h
                ok = True
                for e in range(3):
                    ex = T[i, (e + 1) % 3, 0] - T[i, e, 0]
                    ey = T[i, (e + 1) % 3, 1] - T[i, e, 1]
                    if sign * cross2(ex, ey, px - T[i, e, 0], py - T[i, e, 1]) < -1e-12:
                        ok = False
                        break
                if not ok:
                    continue
                zz = coeff[i, 0] * px + coeff[i, 1] * py + coeff[i, 2]
                k = x * shape[1] + y
                if upper and zz < z[k] or (not upper and zz > z[k]):
                    z[k] = zz
                    own[k] = i
    return (z, own)
