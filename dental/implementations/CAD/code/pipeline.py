"""Finite convex-cap CAD tool. Geometric construction and physical qualification differ."""
from cadlib import *
from margin import detect, cycles
import trimesh
from scipy.spatial import HalfspaceIntersection, ConvexHull, cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from readonly import cone_module

def cap_from_cycle(m, ids):
    banned = {tuple(sorted((int(a), int(b)))) for (a, b) in zip(ids, np.roll(ids, -1))}
    edges = m.face_adjacency_edges
    ok = np.array([tuple(sorted(map(int, e))) not in banned for e in edges])
    adj = m.face_adjacency[ok]
    (n, labels) = connected_components(coo_matrix((np.ones(len(adj)), (adj[:, 0], adj[:, 1])), shape=(len(m.faces), len(m.faces))).tocsr(), directed=False)
    if n != 2:
        raise ValueError('Margin must separate closed mesh into two connected face components')
    target = labels[int(np.argmax(m.triangles_center[:, 2]))]
    mask = labels == target
    return (m.faces[mask], dict(components=2, cap_faces=int(mask.sum()), removed_support_faces=int((~mask).sum()), face_indices=np.flatnonzero(mask)))

def offset_cap(prep, faces, gap):
    tri = prep.vertices[faces]
    norm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    norm /= np.linalg.norm(norm, axis=1)[:, None]
    b = np.einsum('ij,ij->i', norm, tri[:, 0])
    P = prep.vertices
    if (norm @ P.T - b[:, None]).max() > 1e-07:
        raise ValueError('UNKNOWN_NONCONVEX_PREPARATION: halfspace offset would change preparation')
    planes = np.unique(np.round(np.c_[norm, b], 10), axis=0)
    nn = planes[:, :3]
    bb = planes[:, 3]
    nnlen = np.linalg.norm(nn, axis=1)
    bb += gap * nnlen
    z0 = float(P[:, 2].min())
    hs = np.r_[np.c_[nn, -bb], [[0, 0, -1, z0]]]
    inside = P.mean(0)
    inside[2] = (P[:, 2].min() + P[:, 2].max()) / 2
    if np.max(hs[:, :3] @ inside + hs[:, 3]) >= 0:
        raise ValueError('No interior seed for convex offset')
    v = HalfspaceIntersection(hs, inside).intersections
    hull = ConvexHull(v)
    m = trimesh.Trimesh(v, hull.simplices, process=True)
    trimesh.repair.fix_normals(m)
    if m.volume < 0:
        m.invert()
    basal = np.all(abs(m.vertices[m.faces, 2] - z0) < 1e-07, axis=1)
    f = m.faces[~basal]
    edges = np.concatenate([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]])
    (uu, counts) = np.unique(np.sort(edges, axis=1), axis=0, return_counts=True)
    (cs, _) = cycles(uu[counts == 1], m.vertices)
    if len(cs) != 1:
        raise ValueError('Offset mouth not one cycle')
    residual = float(np.max(nn @ P.T - bb[:, None]))
    return (m.vertices, f, cs[0]['indices'], dict(nominal_gap_mm=gap, max_preparation_halfspace_residual_mm=residual, offset_type='convex plane offset; corner gap can exceed nominal', scope='nominal CAD geometry, physical fit UNKNOWN'))

def stitch(ov, of, iv, inf, ib, ob):
    z0 = iv[ib, 2].mean()
    if np.ptp(iv[ib, 2]) > 1e-07 or np.max(abs(ov[ob, 2] - z0)) > 1e-07:
        raise ValueError('UNKNOWN_NONPLANAR_MARGIN')
    ob = np.array(ob)[np.argsort(np.mod(np.arctan2(ov[ob, 1], ov[ob, 0]), 2 * np.pi))]
    ib = np.array(ib)[np.argsort(np.mod(np.arctan2(iv[ib, 1], iv[ib, 0]), 2 * np.pi))]
    oa = np.mod(np.arctan2(ov[ob, 1], ov[ob, 0]), 2 * np.pi)
    ia = np.mod(np.arctan2(iv[ib, 1], iv[ib, 0]), 2 * np.pi)
    i = j = 0
    rim = []
    n = len(ob)
    k = len(ib)
    N = len(ov)
    while i < n or j < k:
        no = oa[(i + 1) % n] + (2 * np.pi if i + 1 >= n else 0) if i < n else np.inf
        ni = ia[(j + 1) % k] + (2 * np.pi if j + 1 >= k else 0) if j < k else np.inf
        if no <= ni:
            rim.append([ob[i % n], N + ib[j % k], ob[(i + 1) % n]])
            i += 1
        else:
            rim.append([ob[i % n], N + ib[j % k], N + ib[(j + 1) % k]])
            j += 1
    v = np.r_[ov, iv]
    f = np.r_[of, inf[:, ::-1] + N, np.asarray(rim, int)]
    roles = np.r_[np.zeros(len(of), int), np.ones(len(inf), int), np.full(len(rim), 2, int)]
    m = trimesh.Trimesh(v, f, process=False)
    trimesh.repair.fix_normals(m)
    return (m, roles)

def create(prep, t, profile, parents, donor_tri):
    d = detect(prep)
    if d['status'] != 'PROPOSED':
        return dict(status='ABSTAIN_MARGIN', margin=d)
    (pf, pc) = cap_from_cycle(prep, np.asarray(d['indices']))
    capm = trimesh.Trimesh(prep.vertices, pf, process=False)
    ax = cone_module().classify(capm.face_normals)
    direction = np.asarray(ax.get('unit_direction', [0, 0, 1]), float)
    dots = capm.face_normals @ direction
    if ax['status'] != 'YES_STRICT':
        return dict(status='ABSTAIN_INSERTION', margin=d, insertion=ax, undercut_dot=dots)
    if direction[2] < 1 - 1e-06:
        raise ValueError('UNKNOWN_AXIS_FRAME_ADAPTER: nonaxial insertion needs full-3D outer solver')
    gap = profile['nominal_internal_gap_mm']
    (iv, inf, ib, offset) = offset_cap(prep, pf, gap)
    top = float(iv[:, 2].max())
    tt = dict(t)
    tt['requirements'] = dict(t['requirements'], wall_mm=profile['wall_mm'])
    tt['preparation_z'] = np.full(len(t['xy']), top - tt['requirements']['film_min_mm'])
    tri = np.array(donor_tri).copy()
    tri[:, :, 0] *= -1
    tri = tri[:, ::-1]
    center = np.mean(tri.reshape(-1, 3), axis=0)
    tri[:, :, :2] -= center[:2]
    donor = trimesh.Trimesh(tri.reshape(-1, 3), np.arange(tri.size // 3).reshape(-1, 3), process=True)
    pv = donor.vertices
    lo = np.quantile(pv[:, 2], 0.3)
    collar = pv[pv[:, 2] <= lo]
    edges = parents['v2g'].boundary_edges(t['faces'])
    bid = np.unique(edges)
    target = np.c_[t['xy'][bid], np.full(len(bid), np.median(collar[:, 2]))]
    scale = float(np.median(np.linalg.norm(target[:, :2], axis=1)) / max(np.median(np.linalg.norm(collar[:, :2], axis=1)), 1e-12))
    donor.vertices[:, :2] *= scale
    collar = donor.vertices[donor.vertices[:, 2] <= lo]
    (R, u, err) = parents['r4'].fit(collar[np.linspace(0, len(collar) - 1, min(512, len(collar)), dtype=int)], target)
    aligned = donor.vertices @ R.T + u
    dtri = aligned[donor.faces]
    prior = parents['v2g'].height(dtri, t['xy'])
    coverage = float(np.isfinite(prior).mean())
    if coverage < 0.5:
        return dict(status='ABSTAIN_DONOR_COVERAGE', coverage=coverage, margin=d, insertion=ax)
    missing = ~np.isfinite(prior)
    if missing.any():
        prior[missing] = prior[~missing][cKDTree(t['xy'][~missing]).query(t['xy'][missing])[1]]
    tt['prior'] = prior
    proposal = parents['gen'].optimize(tt)
    flat = parents['gen'].optimize(dict(tt, prior=np.full(len(prior), top + profile['wall_mm'] + 0.05)))
    if proposal['status'] != 'DESIGN':
        return dict(status='ABSTAIN_ROOF_CONSTRAINTS', margin=d, insertion=ax, reason=proposal.get('reason'), control=flat)
    des = dict(proposal)
    des['inner_vertices'] = np.c_[t['xy'], np.full(len(prior), prep.vertices[:, 2].min() + 3.0)].tolist()
    cap = parents['cap'].complete(tt, des, max(0, profile['wall_mm'] - 0.35) + 0.15)
    (ov, of) = (cap['outer_vertices'], cap['outer_faces'])
    bottom = np.flatnonzero(abs(ov[:, 2] - prep.vertices[:, 2].min()) < 1e-07)
    (m, roles) = stitch(ov, of, iv, inf, ib, bottom)
    return dict(status='GENERATED', mesh=m, roles=roles, margin=d, insertion=ax, cap=pc, undercut_dot=dots, undercut_faces=pf, offset=offset, ov=ov, of=of, iv=iv, inf=inf, outer_roof=np.asarray(proposal['outer_vertices']), prior=prior, donor_fit=dict(R=R, translation_mm=u, scale_xy=scale, start_rms_mm=err, raw_ray_coverage=coverage, missing_nearest_filled=int(missing.sum()), pose='registered source arch, reflected tooth and approximate collar fit'), control=flat, roof_solver=proposal['diagnostics'])
