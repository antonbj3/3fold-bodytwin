"""Continuous sufficient conditions in the represented PL roof. No clinical CAM claim."""
from common import *
from scipy.sparse import csr_matrix

def attach(meta, a):
    t = dict(meta)
    t.update({k: a[k] for k in ['xy', 'faces', 'uv', 'prior', 'ceiling', 'weights']})
    t['A'] = csr_matrix((a['a_data'], a['a_indices'], a['a_indptr']), shape=(int(a['a_rows']), len(a['xy'])))
    t['obstacle_b'] = a['b'] - t['requirements']['clearance_mm']
    t['preparation_z'] = np.full(len(a['xy']), t['preparation_height_mm'])
    return t

def basis(xy):
    t = (xy[:, 0] - xy[:, 0].min()) / np.ptp(xy[:, 0])
    return np.c_[1 - t, t]

def roof_checks(t, z, inner):
    f = t['faces']
    xy = t['xy']
    tri = np.c_[xy, inner][f]
    E = tri[:, 1:] - tri[:, :1]
    det = E[:, 0, 0] * E[:, 1, 1] - E[:, 0, 1] * E[:, 1, 0]
    gx = (E[:, 0, 2] * E[:, 1, 1] - E[:, 1, 2] * E[:, 0, 1]) / det
    gy = (E[:, 0, 0] * E[:, 1, 2] - E[:, 1, 0] * E[:, 0, 2]) / det
    L = float(np.hypot(gx, gy).max())
    lower = float(np.min(z - inner) / np.sqrt(1 + L * L))
    slack = t['obstacle_b'] - t['A'] @ z
    return dict(wall_lower_mm=lower, wall_upper_mm=float(np.min(z - inner)), wall_pass=lower >= t['requirements']['wall_mm'] - 1e-06, intaglio_milling_pass=L < 1e-10, intaglio_Lipschitz=L, antagonist_min_slack_mm=float(slack.min()) if len(slack) else None, antagonist_pass=bool(len(slack) and slack.min() >= -1e-06), scope='Continuous Euclidean separation of PL roof graphs; vertical rim/complete axial crown excluded. Intaglio plane accessible from lower unbounded halfspace with ideal ball.', roundoff_enclosure='MISSING')

def ball_access(xy, faces, z, r=0.3, shaft=0.2):
    if r <= 0 or shaft < 0 or shaft > r:
        return dict(status='UNKNOWN_TOOL_CONTRACT')
    tri = np.c_[xy, z][faces]
    norm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    norm /= np.linalg.norm(norm, axis=1)[:, None]
    norm *= np.where(norm[:, 2] < 0, -1, 1)[:, None]
    if np.min(norm[:, 2]) <= 0:
        return dict(status='UNKNOWN_NON_GRAPH')
    qmin = tri.min(1)
    qmax = tri.max(1)
    bad = []
    ray_bad = []
    minimum = 1000000000.0
    for (i, (p, n)) in enumerate(zip(tri, norm)):
        centers = p + r * n
        cmin = centers.min(0)
        cmax = centers.max(0)
        dbox = np.linalg.norm(np.maximum(0, np.maximum(qmin - cmax, cmin - qmax)), axis=1)
        dplane = float((centers @ n).min()) - (tri @ n).max(1)
        lb = np.maximum(dbox, np.maximum(0, dplane))
        minimum = min(minimum, float(lb.min()))
        if lb.min() < r - 1e-09:
            bad.append(dict(face=i, obstacle_face=int(lb.argmin()), lower_distance_mm=float(lb.min())))
        dx = np.maximum(0, np.maximum(qmin[:, :2] - cmax[:2], cmin[:2] - qmax[:, :2]))
        nearby = np.linalg.norm(dx, axis=1) <= r + 1e-09
        if np.any(qmax[nearby, 2] > cmin[2] + 1e-09):
            ray_bad.append(i)
    return dict(status='CONDITIONAL_IDEAL_TOOL_PASS' if not bad and (not ray_bad) else 'UNKNOWN_NOT_CERTIFIED', tool_radius_mm=r, shaft_radius_mm=shaft, faces=len(faces), uncertified_ball_faces=len(bad), uncertified_vertical_ray_faces=len(ray_bad), minimum_separating_distance_lower_mm=minimum, first_obstructions=bad[:3], scope='Every point of every upper PL face: tangent sphere center triangle, separating plane/AABB distance, vertical approach and semi-infinite coaxial shaft<=ball radius. Rim/intaglio use independent plane certificate. Holder/fixture and full axial crown UNKNOWN.', arithmetic='Analytic sufficient inequalities evaluated in float64 with1e-9mm comparison tolerance; rigorous rounding enclosure MISSING')

def shell(xy, f, z, inner):
    faces = []
    edges = {}
    n = len(xy)
    for (a, b, c) in f:
        if np.cross(xy[b] - xy[a], xy[c] - xy[a]) < 0:
            (b, c) = (c, b)
        faces.extend([[a, b, c], [n + c, n + b, n + a]])
        for (i, j) in [(a, b), (b, c), (c, a)]:
            k = tuple(sorted((int(i), int(j))))
            if k in edges:
                edges.pop(k)
            else:
                edges[k] = (i, j)
    for (i, j) in edges.values():
        faces.extend([[j, i, n + i], [j, n + i, n + j]])
    return (np.vstack([np.c_[xy, z], np.c_[xy, inner]]), np.array(faces, int))

def export_stl(path, vertices, faces):
    import struct
    tris = vertices[faces]
    n = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    n /= np.maximum(np.linalg.norm(n, axis=1)[:, None], 1e-30)
    with Path(path).open('wb') as f:
        f.write(b'X85 millimetres, experimental roof, not clinical crown'.ljust(80, b' '))
        f.write(struct.pack('<I', len(tris)))
        for (no, tt) in zip(n, tris):
            f.write(struct.pack('<12fH', *no, *tt.ravel(), 0))
