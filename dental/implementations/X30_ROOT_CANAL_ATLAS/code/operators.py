"""Frozen R1 measurement operators; annotation geometry, not clinical landmarks."""
import numpy as np
from scipy import ndimage as ndi
from scipy.signal import savgol_filter
H = 0.3
ROOT_RANGE = (0.1, 0.75)

def root_slices(p):
    z = np.argwhere(p)[:, 0]
    (lo, hi) = (int(z.min()), int(z.max()))
    a = int(np.ceil(lo + ROOT_RANGE[0] * (hi - lo)))
    b = int(np.floor(lo + ROOT_RANGE[1] * (hi - lo)))
    return (a, b)

def persistent_count(counts, h=H):
    result = 0
    for k in range(1, max(counts, default=0) + 1):
        runs = ndi.label(np.asarray(counts) >= k)[0]
        longest = max(np.bincount(runs)[1:], default=0)
        if (longest - 1) * h >= 1.5 - 1e-09:
            result = k
    return int(result)

def circumcircle(a, b, c):
    (u, v) = (b - a, c - a)
    cross = float(np.linalg.norm(np.cross(u, v)))
    sides = [float(np.linalg.norm(u)), float(np.linalg.norm(c - b)), float(np.linalg.norm(v))]
    radius = np.prod(sides) / (2 * cross) if cross > 1e-12 else float('inf')
    angle = float(np.degrees(np.arccos(np.clip(np.dot(b - a, c - b) / max(1e-12, sides[0] * sides[1]), -1, 1))))
    return (float(radius), angle)

def path_geometry(points):
    p = np.asarray(points, float)
    length = float(np.linalg.norm(np.diff(p, axis=0), axis=1).sum())
    w = min(7, len(p) if len(p) % 2 else len(p) - 1)
    smooth = savgol_filter(p, w, 2, axis=0) if w >= 3 else p
    triplets = [(i, *circumcircle(smooth[i - 5], smooth[i], smooth[i + 5])) for i in range(5, len(p) - 5)]
    finite = [(i, r, a) for (i, r, a) in triplets if np.isfinite(r)]
    if finite:
        (i, r, a) = min(finite, key=lambda t: t[1])
        selected = smooth[[i - 5, i, i + 5]].tolist()
    else:
        (i, r, a, selected) = (None, None, None, None)
    return dict(annotation_path_length_mm=length, smoothed_path_length_mm=float(np.linalg.norm(np.diff(smooth, axis=0), axis=1).sum()), window_radius_mm=r, window_turn_deg=a, curvature_triplet_mm=selected, end_to_end_mm=float(np.linalg.norm(p[-1] - p[0])), points_mm=p.tolist(), smooth_points_mm=smooth.tolist())

def section_graph(p, origin=(0, 0, 0), h=H):
    (lo, hi) = root_slices(p)
    nodes = []
    edges = []
    counts = []
    prev = None
    ids_prev = {}
    for z in range(lo, hi + 1):
        (labels, n) = ndi.label(p[z], structure=np.ones((3, 3)))
        sizes = np.bincount(labels.ravel())
        keep = [i for i in range(1, n + 1) if sizes[i] >= 2]
        ids = {}
        counts.append(len(keep))
        for i in keep:
            xy = np.argwhere(labels == i).astype(float)
            eig = np.linalg.eigvalsh(np.cov(xy.T, bias=True) + np.eye(2) / 12)
            point = (np.array([z, *xy.mean(0)]) + np.asarray(origin)) * h
            idx = len(nodes)
            ids[i] = idx
            nodes.append(dict(id=idx, z=z, point_mm=point.tolist(), area_mm2=float(sizes[i] * h * h), aspect_ratio=float(np.sqrt(eig[-1] / eig[0])), voxel_count=int(sizes[i])))
        if prev is not None:
            for (old, oid) in ids_prev.items():
                adjacent = ndi.binary_dilation(prev == old, structure=np.ones((3, 3)))
                for i in np.unique(labels[adjacent]):
                    if int(i) in ids:
                        edges.append([oid, ids[int(i)]])
        (prev, ids_prev) = (labels, ids)
    parents = [[] for _ in nodes]
    children = [[] for _ in nodes]
    for (a, b) in edges:
        children[a].append(b)
        parents[b].append(a)
    best = {}
    for n in reversed(nodes):
        idx = n['id']
        options = children[idx]
        if options:
            nxt = max(options, key=lambda j: (len(best[j]), nodes[j]['area_mm2'], -j))
            best[idx] = [idx] + best[nxt]
        else:
            best[idx] = [idx]
    paths = []
    for n in nodes:
        if parents[n['id']]:
            continue
        ids = best[n['id']]
        if (nodes[ids[-1]]['z'] - n['z']) * h < 3 - 1e-09:
            continue
        if nodes[ids[-1]]['z'] < hi - 1:
            continue
        rec = path_geometry([nodes[j]['point_mm'] for j in ids])
        rec['node_ids'] = ids
        rec['apical_annotation_endpoint_mm'] = nodes[ids[0]]['point_mm']
        rec['coronal_cut_endpoint_mm'] = nodes[ids[-1]]['point_mm']
        paths.append(rec)
    return dict(z_range=[lo, hi], nodes=nodes, edges=edges, counts=counts, persistent_count=persistent_count(counts, h), paths=paths, split_nodes=sum((len(c) > 1 for c in children)), merge_nodes=sum((len(a) > 1 for a in parents)))

def measure(p, t, origin=(0, 0, 0), h=H):
    graph = section_graph(p, origin, h)
    (lo, hi) = graph['z_range']
    boundary = t & ~ndi.binary_erosion(t, structure=ndi.generate_binary_structure(3, 1), border_value=0)
    pb = p & ~ndi.binary_erosion(p, structure=ndi.generate_binary_structure(3, 1), border_value=0)
    pb[:lo] = False
    pb[hi + 1:] = False
    pts = np.argwhere(pb)
    dist = ndi.distance_transform_edt(~boundary, sampling=h)[pb]
    root = p.copy()
    root[:lo] = False
    root[hi + 1:] = False
    (_, ncc) = ndi.label(root, structure=np.ones((3, 3, 3)))
    aspects = [n['aspect_ratio'] for n in graph['nodes']]
    areas = [n['area_mm2'] for n in graph['nodes']]
    lengths = [r['annotation_path_length_mm'] for r in graph['paths']]
    radii = [r['window_radius_mm'] for r in graph['paths'] if r['window_radius_mm'] is not None]
    rec = dict(persistent_section_count=graph['persistent_count'], whole_rootward_connected_components=int(ncc), path_count=len(lengths), max_annotation_path_mm=max(lengths, default=None), min_window_radius_mm=min(radii, default=None), section_aspect_median=float(np.median(aspects)) if aspects else None, section_area_median_mm2=float(np.median(areas)) if areas else None, rootward_clearance_min_mm=float(dist.min()) if len(dist) else None, rootward_clearance_p05_mm=float(np.quantile(dist, 0.05)) if len(dist) else None, split_nodes=graph['split_nodes'], merge_nodes=graph['merge_nodes'], clinical_working_length_mm=None, clinical_canal_count=None, Schneider_angle_deg=None, Pruett_radius_mm=None, apex_sinus_mm=None, segmentation_error_mm=None)
    return (rec, graph, pts, np.argwhere(boundary), dist)
