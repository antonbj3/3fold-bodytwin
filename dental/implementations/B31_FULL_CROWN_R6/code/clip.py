"""Reused unchanged R5 construct_e continuous shared-edge clipping; parent hash locked."""
from local import *

def continuous_cut(v, f, values):
    verts = v.tolist()
    out = []
    cache = {}

    def intersect(a, b):
        key = tuple(sorted((int(a), int(b))))
        if key not in cache:
            (aa, bb) = key
            t = values[aa] / (values[aa] - values[bb])
            cache[key] = len(verts)
            verts.append(((1 - t) * v[aa] + t * v[bb]).tolist())
        return cache[key]
    for face in f:
        poly = []
        for (a, b) in zip(face, np.roll(face, -1)):
            if values[a] >= 0:
                poly.append(int(a))
            if (values[a] >= 0) != (values[b] >= 0):
                poly.append(intersect(a, b))
        if len(poly) >= 3:
            out.extend([[poly[0], poly[k], poly[k + 1]] for k in range(1, len(poly) - 1)])
    return compact(np.array(verts), np.array(out, int))
