"""Native signed-field mesh and independent topology/distance audit."""
from collections import defaultdict
import numpy as np
from scipy import ndimage
from scipy.optimize import minimize
from full_geometry import project_cylinder
OFFSETS = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1], [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]])
TETS = [[0, 4, 6, 7], [0, 4, 5, 7], [0, 2, 6, 7], [0, 2, 3, 7], [0, 1, 5, 7], [0, 1, 3, 7]]

def mesh_from_mask(mask, spacing):
    mask = np.pad(mask, 1, constant_values=False)
    sp = np.array(spacing)
    f = ndimage.distance_transform_edt(~mask, sampling=sp) - ndimage.distance_transform_edt(mask, sampling=sp)
    values = np.stack([f[tuple((slice(int(o), int(o) + n - 1) for (o, n) in zip(off, f.shape)))] for off in OFFSETS], axis=-1)
    index = np.argwhere((values.min(axis=-1) < 0) & (values.max(axis=-1) > 0))
    vv = values[tuple(index.T)]
    parts = []
    for tet in TETS:
        val = vv[:, tet]
        codes = ((val < 0) * np.array([1, 2, 4, 8])).sum(axis=1)
        for code in range(1, 15):
            selected = np.flatnonzero(codes == code)
            if not len(selected):
                continue
            pos = (index[selected, None, :] + OFFSETS[tet][None, :, :] - 0.5) * sp
            ft = val[selected]
            neg = [i for i in range(4) if code & 1 << i]
            plus = [i for i in range(4) if i not in neg]

            def zero(i, j):
                return (ft[:, j, None] * pos[:, i, :] - ft[:, i, None] * pos[:, j, :]) / (ft[:, j, None] - ft[:, i, None])
            if len(neg) == 1:
                poly = np.stack([zero(neg[0], j) for j in plus], axis=1)
                tri = poly[:, None, :, :]
            elif len(plus) == 1:
                poly = np.stack([zero(i, plus[0]) for i in neg], axis=1)
                tri = poly[:, None, :, :]
            else:
                (i, k) = neg
                (j, l) = plus
                poly = np.stack([zero(i, j), zero(i, l), zero(k, l), zero(k, j)], axis=1)
                tri = np.stack([poly[:, [0, 1, 2]], poly[:, [0, 2, 3]]], axis=1)
            direction = pos[:, plus, :].mean(axis=1) - pos[:, neg, :].mean(axis=1)
            flat = tri.reshape(-1, 3, 3)
            direction = np.repeat(direction, tri.shape[1], axis=0)
            n = np.cross(flat[:, 1] - flat[:, 0], flat[:, 2] - flat[:, 0])
            reverse = (n * direction).sum(axis=1) < 0
            flat[reverse] = flat[reverse][:, [0, 2, 1]]
            parts.append(flat[:, :, ::-1])
    return np.concatenate(parts).astype(np.float32)

def topology(tri):
    (vertices, idx) = np.unique(tri.reshape(-1, 3), axis=0, return_inverse=True)
    faces = idx.reshape(-1, 3)
    edges = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    orient = np.where(edges[:, 0] < edges[:, 1], 1, -1)
    ordered = np.sort(edges, axis=1)
    (_, inverse, counts) = np.unique(ordered, axis=0, return_inverse=True, return_counts=True)
    oriented_sum = np.bincount(inverse, weights=orient)
    link = [defaultdict(set) for _ in vertices]
    for (a, b, c) in faces:
        for (center, u, v) in [(a, b, c), (b, c, a), (c, a, b)]:
            link[center][u].add(v)
            link[center][v].add(u)
    bad_links = 0
    for graph in link:
        if any((len(v) != 2 for v in graph.values())):
            bad_links += 1
            continue
        todo = [next(iter(graph))]
        seen = set()
        while todo:
            u = todo.pop()
            if u in seen:
                continue
            seen.add(u)
            todo.extend(graph[u] - seen)
        if len(seen) != len(graph):
            bad_links += 1
    return {'triangles': len(tri), 'unique_vertices': len(vertices), 'boundary_edges': int((counts == 1).sum()), 'edges_incidence_not2': int((counts != 2).sum()), 'orientation_inconsistent_edges': int((oriented_sum != 0).sum()), 'bad_vertex_links': bad_links, 'duplicate_triangles': len(faces) - len(np.unique(np.sort(faces, axis=1), axis=0)), 'zero_area_triangles': int((np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1) == 0).sum())}

def triangle_cylinder_bracket(tri, entry, axis, length, radius):
    tri = tri.astype(float)
    centers = tri.mean(axis=1)
    rho = np.linalg.norm(tri - centers[:, None, :], axis=2).max(axis=1)
    dc = np.linalg.norm(centers - project_cylinder(centers, entry, axis, length, radius), axis=1)
    initial = float(dc.min())
    active = np.flatnonzero(dc - rho <= initial + 1e-08)
    records = []
    for i in active:
        verts = tri[i]
        v0 = verts[0]
        basis = (verts[1:] - v0).T

        def fun(uv):
            p = v0 + basis @ uv
            v = p - project_cylinder(p, entry, axis, length, radius)
            d = float(np.linalg.norm(v))
            return (d, basis.T @ (v / d if d else np.zeros(3)))
        ans = minimize(fun, [1 / 3, 1 / 3], jac=True, bounds=[(0, 1), (0, 1)], constraints=[{'type': 'ineq', 'fun': lambda uv: 1 - uv.sum(), 'jac': lambda uv: -np.ones(2)}], method='SLSQP', options={'ftol': 1e-12, 'maxiter': 200})
        uv = np.maximum(ans.x, 0)
        if uv.sum() > 1:
            uv /= uv.sum()
        point = v0 + basis @ uv
        delta = point - project_cylinder(point, entry, axis, length, radius)
        up = float(np.linalg.norm(delta))
        normal = delta / up if up else np.zeros(3)
        low = max(0.0, up + float(((verts - point) @ normal).min()))
        records.append({'triangle_index': int(i), 'vertices_xyz_mm': verts.tolist(), 'barycentric_uv': uv.tolist(), 'point_xyz_mm': point.tolist(), 'normal_xyz': normal.tolist(), 'lower_mm': low, 'upper_mm': up, 'optimizer_success': bool(ans.success)})
    lower = min((r['lower_mm'] for r in records))
    upper = min((r['upper_mm'] for r in records))
    return {'lower_mm': lower, 'upper_mm': upper, 'gap_mm': upper - lower, 'active_triangles': len(active), 'all_triangles': len(tri), 'records': records}

def inside_mesh(point, tri):
    """Generic rays; ambiguity is rejected rather than selecting a convenient ray."""
    tri = tri.astype(float)
    answer = []
    e1 = tri[:, 1] - tri[:, 0]
    e2 = tri[:, 2] - tri[:, 0]
    for vec in [[1, 0.371, 0.127], [0.239, 1, 0.413], [0.311, 0.193, 1]]:
        direction = np.array(vec)
        direction /= np.linalg.norm(direction)
        h = np.cross(direction, e2)
        a = (e1 * h).sum(axis=1)
        ok = np.abs(a) > 1e-12
        inv = np.divide(1, a, out=np.zeros_like(a), where=ok)
        s = point - tri[:, 0]
        u = inv * (s * h).sum(axis=1)
        q = np.cross(s, e1)
        v = inv * (q * direction).sum(axis=1)
        t = inv * (e2 * q).sum(axis=1)
        hits = ok & (u >= 0) & (v >= 0) & (u + v <= 1) & (t > 1e-10)
        answer.append(bool(hits.sum() % 2))
    if len(set(answer)) != 1:
        raise ValueError('Point containment ray ambiguity')
    return answer[0]
