"""Classical ball-envelope control. Exact gates evaluate its final tessellation.

The continuous construction is downward_closure(Q-r*ez) + B_(r+spacer),
using a solid of revolution enclosing the regular-polygon prep. The sampled
profile is not assumed to inherit the continuous construction's guarantees.
"""
import math
import numpy as np
from ..geometry import rounded, SECTORS
from ..checks.exact import q

def profile(preparation, tool_radius=0.5, spacer=0.08):
    zc = np.r_[preparation['z_mm'], preparation['apex_mm']]
    rc = np.r_[np.max(preparation['radii_mm'], axis=1), 0.0]
    R = tool_radius + spacer
    top = zc[-1] + spacer
    zs = np.linspace(0, top, 18)[:-1]
    rr = []
    for z in zs:
        best = 0.0
        for j in range(len(zc) - 1):
            m = (rc[j + 1] - rc[j]) / (zc[j + 1] - zc[j])
            low = max(zc[j], z + tool_radius - R)
            high = min(zc[j + 1], z + tool_radius + R)
            if low > high:
                continue
            optimum = z + tool_radius + m * R / math.sqrt(1 + m * m)
            c = min(high, max(low, optimum))
            radius = rc[j] + m * (c - zc[j]) + math.sqrt(max(0, R * R - (z - c + tool_radius) ** 2))
            best = max(best, radius)
        rr.append(max(best, 0.0001))
    points = [(q(round(float(z), 6)), q(round(float(r), 6))) for (z, r) in zip(zs, rr)] + [(q(round(top, 6)), q(0))]
    hull = []
    for p in points:
        hull.append(p)
        while len(hull) >= 3:
            (a, b, c) = hull[-3:]
            if (b[1] - a[1]) * (c[0] - b[0]) >= (c[1] - b[1]) * (b[0] - a[0]):
                break
            hull.pop(-2)
    return dict(z_mm=[float(p[0]) for p in hull[:-1]], radii_mm=[[float(p[1])] * SECTORS for p in hull[:-1]], apex_mm=round(top, 6), coordinate_model='homothetic_decimal12')
