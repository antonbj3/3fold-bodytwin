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
    out = trimesh.Trimesh(q.vertices, np.vstack([q.faces, cap]), process=False)
    if not out.is_watertight or not out.is_winding_consistent:
        raise ValueError('Exact cap did not preserve closed oriented boundary')
    return out
