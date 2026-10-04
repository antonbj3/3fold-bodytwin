"""D R5 square-frustum dictionary, with its rational box containment rule."""
from common import *
from fractions import Fraction as Q
from scipy import ndimage
from skimage.measure import marching_cubes
from shapely.geometry import box
from geometry import source_mesh, boundary
import trimesh
H = Q(1, 5)
S = Q(float(np.tan(np.deg2rad(3)))).limit_denominator(10 ** 12)

def ceilq(x):
    return -(-x.numerator // x.denominator)

def floorq(x):
    return x.numerator // x.denominator

def erode(a, m):
    return a if m == 0 else ndimage.minimum_filter(a.astype(np.uint8), size=2 * m + 1, mode='constant', cval=0) > 0

def dilate(a, m):
    return a if m == 0 else ndimage.maximum_filter(a.astype(np.uint8), size=2 * m + 1, mode='constant', cval=0) > 0

def dictionary(T):
    (nx, ny, nz) = T.shape
    maxm = ceilq(1 + S * (nz - 1)) - 1
    cache = {(j, m): erode(T[:, :, j], m) for j in range(nz) for m in range(maxm + 1)}
    F = np.zeros_like(T)
    for k in range(nz):
        a = T[:, :, k].copy()
        for j in range(k + 1):
            a &= cache[j, ceilq(1 + S * (k - j)) - 1]
        F[:, :, k] = a
    upper = np.zeros_like(T)
    lower = np.zeros_like(T)
    for k in range(nz):
        for j in range(k + 1):
            upper[:, :, j] |= dilate(F[:, :, k], ceilq(1 + S * (k - j)) - 1)
            if j < k:
                lower[:, :, j] |= dilate(F[:, :, k], floorq(S * (k - j - 1)))
    (ix, iy) = np.nonzero(F.any(2))
    iz = np.array([np.flatnonzero(F[x, y])[-1] for (x, y) in zip(ix, iy)])
    C = np.c_[ix, iy, iz].astype(int)
    return (F, C, lower, upper)

def explicit_check(T, F, C):
    checks = []
    indices = np.argwhere(F)
    for idx in np.linspace(0, len(indices) - 1, min(128, len(indices))).astype(int) if len(indices) else []:
        (x, y, k) = map(int, indices[idx])
        bad = []
        for j in range(k + 1):
            m = ceilq(1 + S * (k - j)) - 1
            for xx in range(x - m, x + m + 1):
                for yy in range(y - m, y + m + 1):
                    if not (0 <= xx < T.shape[0] and 0 <= yy < T.shape[1]) or not T[xx, yy, j]:
                        bad.append([xx, yy, j])
        checks.append({'apex_index': [x, y, k], 'outside_boxes': bad})
    return checks

def union_surface(C, origin, base, pitch=0.1):
    coords = origin + float(H) * C
    coords[:, 2] = origin[2] + float(H) * C[:, 2]
    r0 = float(H / 2 - S * H / 2)
    s = float(S)
    maxr = r0 + s * (coords[:, 2].max() - base)
    lo = np.r_[coords[:, :2].min(0) - maxr - pitch * 2, base - pitch * 2]
    hi = np.r_[coords[:, :2].max(0) + maxr + pitch * 2, coords[:, 2].max() + pitch * 2]
    axes = [np.arange(np.floor(lo[i] / pitch) * pitch, np.ceil(hi[i] / pitch) * pitch + pitch / 2, pitch) for i in range(3)]
    phi = np.full(tuple((len(a) for a in axes)), np.inf)
    for (x, y, a) in coords:
        rad = r0 + s * (a - base)
        bounds = [(x - rad - pitch, x + rad + pitch), (y - rad - pitch, y + rad + pitch), (base - pitch, a + pitch)]
        ss = [slice(max(0, np.searchsorted(ax, b[0]) - 1), min(len(ax), np.searchsorted(ax, b[1]) + 1)) for (ax, b) in zip(axes, bounds)]
        xx = axes[0][ss[0]][:, None, None]
        yy = axes[1][ss[1]][None, :, None]
        zz = axes[2][ss[2]][None, None, :]
        q = np.maximum(np.maximum(abs(xx - x) - r0 - s * (a - zz), abs(yy - y) - r0 - s * (a - zz)), np.maximum(zz - a, base - zz))
        phi[tuple(ss)] = np.minimum(phi[tuple(ss)], q)
    phi = np.minimum(phi, pitch * 4).astype(np.float32)
    (v, f, _, _) = marching_cubes(phi, 0, spacing=(pitch, pitch, pitch), allow_degenerate=False)
    v += np.array([a[0] for a in axes])
    m = trimesh.Trimesh(v, f, process=False)
    trimesh.repair.fix_normals(m)
    if m.volume < 0:
        m.invert()
    return (m, coords, {'pitch_mm': pitch, 'field_shape': phi.shape, 'r0_mm': r0, 'slope_exact': str(S), 'base_mm': base, 'field_is_SDF': False, 'redistance': 'NOT_USED_NO_METRIC_SDF_OPERATION; analytic halfspace union field'})

def parameter(v):
    l = np.linalg.norm(np.roll(v, -1, axis=0) - v, axis=1)
    return np.r_[0, np.cumsum(l[:-1])] / l.sum()

def crown_shell(original, closed_inner, base):
    o = source_mesh(original)
    io = boundary(o)
    rim = o.vertices[io]
    m = closed_inner.copy()
    keep = ~np.all(abs(m.triangles[:, :, 2] - base) < 2e-06, axis=1)
    m.update_faces(keep)
    m.remove_unreferenced_vertices()
    ii = boundary(m)
    inner = m.vertices[ii]
    if np.max(abs(inner[:, 2] - base)) > 2e-06:
        raise ValueError('UNRESOLVED_INNER_BASE_RING')
    if PolygonSign(rim[:, :2]) * PolygonSign(inner[:, :2]) < 0:
        ii = ii[::-1]
        inner = m.vertices[ii]
    j = int(np.argmin(np.linalg.norm(inner - rim[0], axis=1)))
    ii = np.roll(ii, -j)
    inner = m.vertices[ii]
    u = parameter(rim)
    w = parameter(inner)
    N = len(o.vertices)
    v = np.r_[o.vertices, m.vertices]
    f = o.faces.tolist() + (m.faces + N).tolist()
    roles = [0] * len(o.faces) + [1] * len(m.faces)
    a = b = 0
    while a < len(io) or b < len(ii):
        ua = u[a + 1] if a + 1 < len(io) else 1.0
        wb = w[b + 1] if b + 1 < len(ii) else 1.0
        if a < len(io) and (b == len(ii) or ua <= wb):
            f.append([int(io[a % len(io)]), int(io[(a + 1) % len(io)]), int(N + ii[b % len(ii)])])
            a += 1
        else:
            f.append([int(io[a % len(io)]), int(N + ii[(b + 1) % len(ii)]), int(N + ii[b % len(ii)])])
            b += 1
        roles.append(2)
    out = trimesh.Trimesh(v, np.array(f), process=False)
    trimesh.repair.fix_normals(out)
    if out.volume < 0:
        out.invert()
    return (out, np.array(roles), m)

def PolygonSign(p):
    return np.sum(p[:, 0] * np.roll(p[:, 1], -1) - p[:, 1] * np.roll(p[:, 0], -1))
