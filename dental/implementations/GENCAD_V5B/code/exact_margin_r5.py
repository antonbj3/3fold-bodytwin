import numpy as np
import trimesh
import shapely
from scipy.spatial import cKDTree

def clip_cap(m, margin):
    q = trimesh.intersections.slice_mesh_plane(m, [0, 0, 1], [0, 0, margin], cap=False)
    q.merge_vertices(digits_vertex=10)
    edges = np.sort(np.vstack([q.faces[:, [0, 1]], q.faces[:, [1, 2]], q.faces[:, [2, 0]]]), axis=1)
    (unique, count) = np.unique(edges, axis=0, return_counts=True)
    boundary = unique[count == 1]
    if not len(boundary):
        raise ValueError('Missing exact cut boundary')
    nodes = np.unique(boundary)
    if np.max(abs(q.vertices[nodes, 2] - margin)) > 1e-06:
        raise ValueError('Nonplanar cut boundary')
    q.vertices[nodes, 2] = margin
    polys = trimesh.path.polygons.edges_to_polygons(boundary, q.vertices[:, :2])
    tree = cKDTree(q.vertices[nodes, :2])
    cap = []
    for poly in polys:
        triangles = shapely.constrained_delaunay_triangles(poly)
        for tri in triangles.geoms:
            xy = np.asarray(tri.exterior.coords)[:3]
            (dist, ix) = tree.query(xy)
            if np.max(dist) > 1e-07:
                raise ValueError('Cap invented a boundary vertex')
            ids = nodes[ix]
            if np.cross(xy[1] - xy[0], xy[2] - xy[0]) > 0:
                ids = ids[::-1]
            cap.append(ids)
    if not cap:
        raise ValueError('Empty annular cap')
    adj = {int(i): [] for i in nodes}
    for (a, b) in boundary:
        adj[int(a)].append(int(b))
        adj[int(b)].append(int(a))
    original_edges = {tuple(sorted(map(int, e))) for e in boundary}
    cf = np.asarray(cap, int)
    ee = np.sort(np.vstack([cf[:, [0, 1]], cf[:, [1, 2]], cf[:, [2, 0]]]), axis=1)
    (uu, cc) = np.unique(ee, axis=0, return_counts=True)
    paths = {}
    for (a, b) in uu[cc == 1]:
        (a, b) = (int(a), int(b))
        key = tuple(sorted((a, b)))
        if key in original_edges:
            continue
        aa = q.vertices[a, :2]
        d = q.vertices[b, :2] - aa
        l2 = float(d @ d)
        found = None
        for first in adj[a]:
            path = [a, first]
            previous = a
            cur = first
            for _ in range(len(nodes)):
                if cur == b:
                    found = path
                    break
                xy = q.vertices[cur, :2]
                u = float((xy - aa) @ d / max(l2, 1e-300))
                distance = float(np.linalg.norm(xy - aa - u * d))
                if not 0 < u < 1 or distance > 1e-08:
                    break
                nxt = [x for x in adj[cur] if x != previous]
                if len(nxt) != 1 or nxt[0] in path:
                    break
                (previous, cur) = (cur, nxt[0])
                path.append(cur)
            if found is not None:
                break
        if found is None:
            raise ValueError('Cap boundary is not a collinear subdivision of original boundary')
        paths[key] = found
    pending = list(cap)
    cap = []
    while pending:
        ids = np.asarray(pending.pop(), int)
        split = False
        for j in range(3):
            (a, b, c) = map(int, [ids[j], ids[(j + 1) % 3], ids[(j + 2) % 3]])
            key = tuple(sorted((a, b)))
            if key in paths:
                edge = paths.pop(key)
                if edge[0] != a:
                    edge = edge[::-1]
                pending.extend([[edge[k], edge[k + 1], c] for k in range(len(edge) - 1)])
                split = True
                break
        if not split:
            cap.append(ids)
    out = trimesh.Trimesh(q.vertices, np.vstack([q.faces, cap]), process=False)
    if not out.is_watertight or not out.is_winding_consistent:
        raise ValueError('Exact cap did not preserve closed oriented boundary')
    return out
