from common import *
from scipy.sparse import csr_matrix, coo_matrix, diags
from scipy.spatial import cKDTree
import struct

def stencil(v, f, xy):
    """Every exterior hit and barycentric owner, immutable XY footprint."""
    tri = v[f]
    c = tri[:, :, :2].mean(1)
    radius = np.linalg.norm(tri[:, :, :2] - c[:, None, :], axis=2).max(1)
    hits = cKDTree(c).query_ball_point(xy, float(radius.max()), workers=1)
    rows = []
    cols = []
    vals = []
    z = []
    qr = []
    owner = np.full(len(xy), -1, int)
    height = np.full(len(xy), np.nan)
    for (i, js) in enumerate(hits):
        if not js:
            continue
        t = tri[js]
        (a, b, c) = (t[:, 0, :2], t[:, 1, :2], t[:, 2, :2])
        den = np.cross(b - a, c - a)
        ok = np.abs(den) > 1e-12
        p = xy[i]
        w1 = np.divide(np.cross(b - p, c - p), den, out=np.zeros(len(t)), where=ok)
        w2 = np.divide(np.cross(c - p, a - p), den, out=np.zeros(len(t)), where=ok)
        w3 = 1 - w1 - w2
        ok &= (w1 >= -1e-09) & (w2 >= -1e-09) & (w3 >= -1e-09)
        local = []
        for j in np.flatnonzero(ok):
            weights = np.maximum(0.0, [w1[j], w2[j], w3[j]])
            weights /= weights.sum()
            h = float(weights @ t[j, :, 2])
            rid = len(z)
            rows.extend([rid] * 3)
            cols.extend(f[js[j]])
            vals.extend(weights)
            z.append(h)
            qr.append(i)
            local.append(rid)
        if local:
            best = max(local, key=lambda k: z[k])
            owner[i] = best
            height[i] = z[best]
    H = csr_matrix((vals, (rows, cols)), shape=(len(z), len(v)))
    return (H, np.array(qr, int), np.array(z), owner, height)

def prepare(r, contact, reuse):
    with np.load(r['mesh_path'], allow_pickle=False) as a:
        v = a['vertices'].copy()
        f = a['faces']
        roles = a['face_roles']
    with np.load(r['preparation_path'], allow_pickle=False) as a:
        margin = float(a['margin_z'])
    with np.load(r['reference_path'], allow_pickle=False) as a:
        source = dict(a)
    tri = source['source_triangles']
    tri = tri[tri.mean(1)[:, 2] >= margin]
    (xy, ff) = reuse.grid(tri)
    ceil = reuse.VF.height(source['antagonist_triangles'], xy, True)
    rg = ceil - reuse.VF.height(tri, xy)
    extf = f[roles == 0]
    (H, qr, z, owner, top) = stencil(v, extf, xy)
    g0 = ceil - top
    legacy_g = ceil - reuse.VF.height(v[extf], xy)
    finite = np.isfinite(g0) & np.isfinite(legacy_g)
    parity = float(abs(g0[finite] - legacy_g[finite]).max()) if finite.any() else 0.0
    protected = np.unique(f[roles != 0])
    nz = np.cross(v[extf[:, 1], :2] - v[extf[:, 0], :2], v[extf[:, 2], :2] - v[extf[:, 0], :2])
    nonupper = np.unique(extf[nz <= 1e-12])
    editable = np.setdiff1d(np.unique(extf), np.union1d(protected, nonupper))
    taper = np.clip((v[:, 2] - margin) / 0.6, 0, 1)
    taper[np.setdiff1d(np.arange(len(v)), editable)] = 0.0
    area = np.linalg.norm(np.cross(v[extf[:, 1]] - v[extf[:, 0]], v[extf[:, 2]] - v[extf[:, 0]]), axis=1) / 2
    vw = np.bincount(extf.ravel(), weights=np.repeat(area / 3, 3), minlength=len(v))
    commonfaces = ff[np.isfinite(g0[ff]).all(1) & np.isfinite(rg[ff]).all(1)]
    qgood = np.zeros(len(xy), bool)
    qgood[np.unique(commonfaces)] = True
    weights = np.bincount(commonfaces.ravel(), weights=np.repeat(abs(np.cross(xy[commonfaces[:, 1]] - xy[commonfaces[:, 0]], xy[commonfaces[:, 2]] - xy[commonfaces[:, 0]])) / 6, 3), minlength=len(xy))
    ndown = nz < -1e-12
    max_down = float(v[extf[ndown], 2].max()) if ndown.any() else None
    dest = DATA / 'inputs' / (r['key'] + '.npz')
    dest.parent.mkdir(exist_ok=True, parents=True)
    np.savez_compressed(dest, vertices=v, faces=f, face_roles=roles, xy=xy, grid_faces=ff, ceiling=ceil, reference_gap=rg, original_gap=g0, margin_z=margin, protected=protected, editable=editable, taper=taper, vertex_area_weights=vw, query_weights=weights, H_data=H.data, H_indices=H.indices, H_indptr=H.indptr, H_shape=H.shape, hit_queries=qr, hit_z=z, owner=owner, query_good=qgood, max_downward_exterior_z=max_down if max_down is not None else -1e+30)
    return dict(key=r['key'], path=str(dest), sha256=sha(dest), vertices=len(v), faces=len(f), editable_vertices=int(np.count_nonzero(taper)), grid_points=len(xy), common_faces=len(commonfaces), hit_rows=len(z), legacy_gap_max_difference_mm=parity, reference_positive_band_area_mm2=contact.compare(xy, ff, g0, rg).get('reference', {}).get('area_mm2'))

def load_case(r):
    p = DATA / 'inputs' / (r['key'] + '.npz')
    with np.load(p, allow_pickle=False) as a:
        t = {k: a[k] for k in a.files}
    t['H'] = csr_matrix((t.pop('H_data'), t.pop('H_indices'), t.pop('H_indptr')), shape=tuple(t.pop('H_shape')))
    return t

def gap(t, d):
    g = np.full(len(t['xy']), np.inf)
    qr = t['hit_queries']
    keep = np.isfinite(t['ceiling'][qr])
    val = t['ceiling'][qr] - t['hit_z'] + t['H'] @ d
    np.minimum.at(g, qr[keep], val[keep])
    g[np.isinf(g)] = np.nan
    return g

def field_basis(t, regional=True):
    n = len(t['vertices'])
    used = np.flatnonzero(t['taper'] > 0)
    if not regional:
        return (csr_matrix((t['taper'][used], (used, np.arange(len(used)))), shape=(n, len(used))), dict(kind='ordinary_mesh_vertex', coefficients=len(used)))
    h = 0.75
    lo = t['vertices'][:, :2].min(0)
    p = (t['vertices'][used, :2] - lo) / h
    ij = np.floor(p).astype(int)
    uv = p - ij
    shape = np.ceil((t['vertices'][:, :2].max(0) - lo) / h).astype(int) + 2
    rows = []
    cols = []
    vals = []
    for (dx, dy) in [(0, 0), (1, 0), (0, 1), (1, 1)]:
        w = (uv[:, 0] if dx else 1 - uv[:, 0]) * (uv[:, 1] if dy else 1 - uv[:, 1]) * t['taper'][used]
        rows.extend(used)
        cols.extend((ij[:, 0] + dx) * shape[1] + ij[:, 1] + dy)
        vals.extend(w)
    B = csr_matrix((vals, (rows, cols)), shape=(n, int(np.prod(shape))))
    active = np.flatnonzero(np.asarray(B.sum(0)).ravel() > 0)
    B = B[:, active]
    return (B, dict(kind='regional_bilinear', coefficients=B.shape[1], spacing_mm=h))

def wall_and_removal(t, d, r):
    v = t['vertices']
    f = t['faces']
    extf = f[t['face_roles'] == 0]
    new = v.copy()
    new[:, 2] -= d
    protected = t['protected']
    identity = float(np.max(abs(new[protected] - v[protected])))
    moving = np.any(d[extf] > 0, axis=1)
    nz = np.cross(v[extf[:, 1], :2] - v[extf[:, 0], :2], v[extf[:, 2], :2] - v[extf[:, 0], :2])
    up = bool(np.all(nz[moving] > 0))
    lowest = float(new[extf[moving], 2].min()) if moving.any() else None
    down = float(t['max_downward_exterior_z'])
    separator = not moving.any() or lowest > down
    cap = float(d.max(initial=0))
    lb = r['wall_lower_mm'] - cap
    return (new, dict(protected_identity_error_mm=identity, wall_lower_mm=lb, max_relief_mm=cap, area_integral_abs_displacement_mm3=float(t['vertex_area_weights'] @ d), wall_model_pass=lb >= 0.5 and identity == 0.0, downward_vertex_only=bool(np.all(d >= -1e-12)), moving_faces_upper=up, external_lower_boundary_separator_pass=separator, lowest_moving_face_z_mm=lowest, max_fixed_downward_exterior_z_mm=down, removal_status='CONDITIONAL_SUBSET' if up and separator and (lb > 0) else 'UNKNOWN_FULL_SOLID_SUBSET', scope='Represented exterior, basal annulus excluded; fixed intaglio. Sufficient vertical-boundary separator; original shell validity inherited.', rigorous_machine_enclosure='MISSING', anatomical_hausdorff_upper_mm=cap))

def export_stl(p, v, f):
    dt = np.dtype([('normal', '<f4', (3,)), ('vertices', '<f4', (3, 3)), ('attr', '<u2')])
    out = np.zeros(len(f), dt)
    tri = v[f]
    normal = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    normal /= np.maximum(np.linalg.norm(normal, axis=1)[:, None], 1e-30)
    out['normal'] = normal
    out['vertices'] = tri
    p = Path(p)
    p.parent.mkdir(exist_ok=True, parents=True)
    p.write_bytes(b'X95 private experimental geometry; mm; no clinical validation'.ljust(80, b' ') + struct.pack('<I', len(f)) + out.tobytes())
