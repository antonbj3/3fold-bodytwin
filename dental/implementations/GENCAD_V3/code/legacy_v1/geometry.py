import numpy as np
import trimesh
from .checks.exact import vec, sub, cross, dot, norm2, strings
SECTORS = 16
THETA = np.array([0.0, 0.22, 0.44, 0.66, 0.88, 1.1, 1.32])

def rounded(x):
    return np.round(np.asarray(x, dtype=float), 6).tolist()

def dome_profile(radius, height):
    return dict(z_mm=rounded(height * np.sin(THETA)), radii_mm=rounded(np.repeat((radius * np.cos(THETA))[:, None], SECTORS, axis=1)), apex_mm=round(float(height), 6))

def cap(profile):
    z = np.asarray(profile['z_mm'], float)
    rr = np.asarray(profile['radii_mm'], float)
    if rr.shape != (len(z), SECTORS) or len(z) < 2 or z[0] != 0 or np.any(np.diff(z) <= 0):
        raise ValueError('Invalid ring topology')
    if not np.isfinite(rr).all() or not np.isfinite(z).all() or np.min(rr) <= 0:
        raise ValueError('Nonfinite or nonpositive geometry')
    H = float(profile['apex_mm'])
    if not np.isfinite(H) or H <= z[-1]:
        raise ValueError('Invalid apex')
    a = np.arange(SECTORS) * 2 * np.pi / SECTORS
    precision = profile.get('coordinate_model', 'decimal6')
    if precision == 'homothetic_decimal12':
        cx = np.round(np.cos(a), 6)
        cy = np.round(np.sin(a), 6)
        v = np.c_[(rr * cx).ravel(), (rr * cy).ravel(), np.repeat(z, SECTORS)]
        digits = 12
    elif precision == 'decimal6':
        v = np.c_[(rr * np.cos(a)).ravel(), (rr * np.sin(a)).ravel(), np.repeat(z, SECTORS)]
        digits = 6
    else:
        raise ValueError('Unsupported coordinate model')
    v = np.round(np.vstack([v, [0, 0, H]]), digits)
    P = list(map(vec, v))
    f = []
    for j in range(len(z) - 1):
        for k in range(SECTORS):
            a0 = j * SECTORS + k
            b = j * SECTORS + (k + 1) % SECTORS
            c = b + SECTORS
            d = a0 + SECTORS
            n = cross(sub(P[b], P[a0]), sub(P[c], P[a0]))
            if dot(n, sub(P[d], P[a0])) > 0:
                f.extend([[a0, b, d], [b, c, d]])
            else:
                f.extend([[a0, b, c], [a0, c, d]])
    for k in range(SECTORS):
        f.append([(len(z) - 1) * SECTORS + k, (len(z) - 1) * SECTORS + (k + 1) % SECTORS, len(v) - 1])
    return (v, np.array(f, dtype=int))

def shell(outer, inner):
    (ov, of) = cap(outer)
    (iv, inf) = cap(inner)
    n = len(ov)
    f = np.vstack([of, inf[:, ::-1] + n])
    rim = []
    for k in range(SECTORS):
        j = (k + 1) % SECTORS
        rim.extend([[k, n + j, j], [k, n + k, n + j]])
    return trimesh.Trimesh(np.vstack([ov, iv]), np.vstack([f, rim]), process=False)

def planes_from_cap(profile, require_convex=True):
    """Exact supporting planes of the submitted triangulation, open toward -z.
    All faces must be supporting: no float-convexification of the submitted shape.
    """
    (v, f) = cap(profile)
    P = list(map(vec, v))
    planes = []
    seen = set()
    for tri in f:
        (a, b, c) = [P[k] for k in tri]
        normal = cross(sub(b, a), sub(c, a))
        bb = dot(normal, a)
        if not norm2(normal):
            raise ValueError('Degenerate face')
        scale = max((abs(x) for x in normal))
        normal = tuple((x / scale for x in normal))
        bb /= scale
        if normal[2] < 0:
            raise ValueError('Cavity has an overhang')
        if require_convex and any((dot(normal, p) > bb for p in P)):
            raise ValueError('Cap is not exactly convex')
        row = tuple(strings(normal + (bb,)))
        if row not in seen:
            seen.add(row)
            planes.append(list(row))
    return planes

def profile_from_points(points):
    """Fixed coarse anatomical representation; no target-dependent registration."""
    p = np.asarray(points, float)
    H = float(np.quantile(p[:, 2], 0.995))
    if H <= 0.5:
        raise ValueError('Insufficient crown height')
    angle = np.arctan2(p[:, 1], p[:, 0])
    rad = np.linalg.norm(p[:, :2], axis=1)
    zs = np.sin(THETA) * H
    out = []
    for zz in zs:
        ring = []
        for aa in np.arange(SECTORS) * 2 * np.pi / SECTORS:
            da = np.angle(np.exp(1j * (angle - aa)))
            dist = (da / 0.5) ** 2 + ((p[:, 2] - zz) / (H * 0.18)) ** 2
            ii = np.argsort(dist, kind='stable')[:max(5, min(30, len(p) // 30))]
            ring.append(max(0.1, float(np.median(rad[ii]))))
        out.append(ring)
    return dict(z_mm=rounded(zs), radii_mm=rounded(out), apex_mm=round(H + 0.01, 6))

def surface_samples(v, f, count=2048, seed=6101):
    v = np.asarray(v, float)
    f = np.asarray(f, int)
    t = v[f]
    area = np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1)
    if not np.isfinite(area).all() or area.sum() <= 0:
        raise ValueError('Empty surface')
    rng = np.random.default_rng(seed)
    ii = rng.choice(len(f), count, p=area / area.sum())
    uv = rng.random((count, 2))
    flip = uv.sum(1) > 1
    uv[flip] = 1 - uv[flip]
    return t[ii, 0] + uv[:, 0, None] * (t[ii, 1] - t[ii, 0]) + uv[:, 1, None] * (t[ii, 2] - t[ii, 0])
