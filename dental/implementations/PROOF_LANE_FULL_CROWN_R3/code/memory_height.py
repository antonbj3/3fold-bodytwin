"""Bound candidate lists; preserve the frozen affine triangle calculation/order."""
import numpy as np
from scipy.spatial import cKDTree

def height(tri, xy, upper=False):
    out = np.full(len(xy), np.nan)
    if not len(tri):
        return out
    centers = tri[:, :, :2].mean(1)
    radii = np.linalg.norm(tri[:, :, :2] - centers[:, None, :], axis=2).max(1)
    tree = cKDTree(centers)
    for start in range(0, len(xy), 16):
        hits = tree.query_ball_point(xy[start:start + 16], float(radii.max()), workers=1)
        for (off, js) in enumerate(hits):
            if not js:
                continue
            i = start + off
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
