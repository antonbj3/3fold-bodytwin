"""Continuous finite-cylinder distance to the union of annotated voxel boxes."""
import numpy as np
from scipy.optimize import minimize

def project_cylinder(points, entry, axis, length, radius):
    points = np.asarray(points, float)
    rel = points - entry
    z = rel @ axis
    lat = rel - z[..., None] * axis
    rr = np.linalg.norm(lat, axis=-1)
    factor = np.minimum(1, np.divide(radius, rr, out=np.ones_like(rr), where=rr > 0))
    return entry + np.clip(z, 0, length)[..., None] * axis + lat * factor[..., None]

def voxel_cylinder_bracket(centers, spacing, entry, axis, length, radius, tol=1e-08):
    half = np.asarray(spacing, float) / 2
    dist = np.linalg.norm(centers - project_cylinder(centers, entry, axis, length, radius), axis=1)
    initial_upper = float(dist.min())
    halfdiag = float(np.linalg.norm(half))
    active = np.flatnonzero(dist - halfdiag <= initial_upper + tol)
    records = []
    for i in active:
        center = centers[i]
        lo = center - half
        hi = center + half

        def value_grad(p):
            q = project_cylinder(p, entry, axis, length, radius)
            v = p - q
            d = float(np.linalg.norm(v))
            return (d, v / d if d > 0 else np.zeros(3))
        opt = minimize(value_grad, center, jac=True, bounds=list(zip(lo, hi)), method='L-BFGS-B', options={'ftol': 1e-15, 'gtol': 1e-11, 'maxiter': 200, 'maxls': 50})
        pp = np.clip(opt.x, lo, hi)
        (up, normal) = value_grad(pp)
        corner = np.where(normal >= 0, lo, hi)
        lower = max(0.0, up + float(normal @ (corner - pp)))
        records.append({'voxel_index_in_center_array': int(i), 'center_zyx_mm': center.tolist(), 'witness_zyx_mm': pp.tolist(), 'normal_zyx': normal.tolist(), 'lower_mm': lower, 'upper_mm': up, 'optimizer_success': bool(opt.success), 'optimizer_message': str(opt.message)})
    lower = min((r['lower_mm'] for r in records))
    upper = min((r['upper_mm'] for r in records))
    return {'lower_mm': lower, 'upper_mm': upper, 'gap_mm': upper - lower, 'active_voxel_boxes': len(active), 'all_voxel_boxes': len(centers), 'initial_center_upper_mm': initial_upper, 'halfdiag_mm': halfdiag, 'witness': min(records, key=lambda r: r['upper_mm']), 'records': records}

def write_stl(path, triangles):
    triangles = np.asarray(triangles, np.float32)
    normal = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    nn = np.linalg.norm(normal, axis=1)
    normal /= np.maximum(nn, 1e-30)[:, None]
    records = np.zeros(len(triangles), dtype=[('normal', '<f4', (3,)), ('vertices', '<f4', (3, 3)), ('attribute', '<u2')])
    records['normal'] = normal
    records['vertices'] = triangles
    with open(path, 'wb') as f:
        f.write(b'X8 RESEARCH PHANTOM; mm; xyz; annotation-derived'.ljust(80, b' '))
        f.write(np.array([len(triangles)], dtype='<u4').tobytes())
        f.write(records.tobytes())
    return len(triangles)

def voxel_triangles(mask, spacing):
    """Exact exposed faces of an occupancy grid, local xyz millimetres.

    Voxel centers are (index+.5)*spacing. Exports may retain nonmanifold diagonal
    contacts; STL is a fabrication candidate until the slicer/metrology reviews it.
    """
    chunks = []
    for axis in range(3):
        for sign in [-1, 1]:
            neighbor = np.zeros_like(mask)
            a = [slice(None)] * 3
            b = [slice(None)] * 3
            if sign == 1:
                a[axis] = slice(None, -1)
                b[axis] = slice(1, None)
            else:
                a[axis] = slice(1, None)
                b[axis] = slice(None, -1)
            neighbor[tuple(a)] = mask[tuple(b)]
            idx = np.argwhere(mask & ~neighbor)
            if len(idx) == 0:
                continue
            xyz = idx[:, ::-1].astype(float)
            ax = 2 - axis
            other = [i for i in range(3) if i != ax]
            base = xyz.copy()
            base[:, ax] += 1 if sign == 1 else 0
            u = np.zeros(3)
            v = np.zeros(3)
            u[other[0]] = 1
            v[other[1]] = 1
            corners = np.stack([base, base + u, base + u + v, base + v], 1) * np.asarray(spacing)[::-1]
            tri = np.concatenate([corners[:, [0, 1, 2]], corners[:, [0, 2, 3]]])
            n = np.cross(tri[0, 1] - tri[0, 0], tri[0, 2] - tri[0, 0])
            if np.sign(n[ax]) != sign:
                tri = tri[:, [0, 2, 1]]
            chunks.append(tri)
    return np.concatenate(chunks)

def cylinder_triangles(entry_xyz, axis_xyz, length, radius, facets=64):
    a = np.asarray(axis_xyz)
    e = np.asarray(entry_xyz)
    u = np.cross(a, [1, 0, 0] if abs(a[0]) < 0.9 else [0, 1, 0])
    u /= np.linalg.norm(u)
    v = np.cross(a, u)
    phi = np.linspace(0, 2 * np.pi, facets, endpoint=False)
    ring = radius * (np.cos(phi)[:, None] * u + np.sin(phi)[:, None] * v)
    p = e + ring
    q = p + length * a
    tris = []
    for i in range(facets):
        j = (i + 1) % facets
        tris.extend([[p[i], p[j], q[j]], [p[i], q[j], q[i]], [e, p[j], p[i]], [e + length * a, q[i], q[j]]])
    return np.array(tris)
