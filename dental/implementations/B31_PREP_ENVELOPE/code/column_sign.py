from common import *
from geometry import distances

def column_sign(closed, axes):
    (nx, ny, nz) = map(len, axes)
    hits = {}
    projected_vertical = 0
    for tri in closed.triangles:
        (a, b, c) = tri
        u = b[:2] - a[:2]
        v = c[:2] - a[:2]
        den = u[0] * v[1] - u[1] * v[0]
        if abs(den) < 1e-12:
            projected_vertical += 1
            continue
        ix = np.flatnonzero((axes[0] >= tri[:, 0].min() - 1e-12) & (axes[0] <= tri[:, 0].max() + 1e-12))
        iy = np.flatnonzero((axes[1] >= tri[:, 1].min() - 1e-12) & (axes[1] <= tri[:, 1].max() + 1e-12))
        if not len(ix) or not len(iy):
            continue
        (xx, yy) = np.meshgrid(axes[0][ix], axes[1][iy], indexing='ij')
        dx = xx - a[0]
        dy = yy - a[1]
        be = (dx * v[1] - dy * v[0]) / den
        ga = (u[0] * dy - u[1] * dx) / den
        mask = (be >= -1e-12) & (ga >= -1e-12) & (be + ga <= 1 + 1e-12)
        zz = a[2] + be * (b[2] - a[2]) + ga * (c[2] - a[2])
        for (i, j) in np.argwhere(mask):
            hits.setdefault((int(ix[i]), int(iy[j])), []).append(float(zz[i, j]))
    inside = np.zeros((nx, ny, nz), bool)
    dedup = 0
    for ((i, j), vals) in hits.items():
        vals = np.sort(vals)
        unique = vals[np.r_[True, np.diff(vals) > 1e-10]]
        dedup += len(vals) - len(unique)
        inside[i, j] = (len(unique) - np.searchsorted(unique, axes[2], side='right')) % 2 == 1
    pts = np.stack(np.meshgrid(*axes, indexing='ij'), -1).reshape(-1, 3)
    ids = np.linspace(0, len(pts) - 1, min(256, len(pts))).astype(int)
    control = np.r_[closed.contains(pts[ids[:128]]), closed.contains(pts[ids[128:]])]
    ours = inside.reshape(-1)[ids]
    away = distances(closed, pts[ids]) > 1e-07
    mis = (ours != control) & away
    return (inside.reshape(-1), {'control_points': len(ids), 'away_surface_control_points': int(away.sum()), 'interior_control_mismatches': int(mis.sum()), 'mismatch_points_mm': pts[ids[mis]], 'projected_vertical_facets': projected_vertical, 'crossing_duplicates_removed': dedup, 'arithmetic_enclosure': 'MISSING; floating parity under declared tolerance', 'physical_inside': 'UNKNOWN_SOURCE_CLOSURE'})
