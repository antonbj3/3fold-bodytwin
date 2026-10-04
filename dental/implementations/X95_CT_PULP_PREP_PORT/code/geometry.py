"""Discrete complete solids and exact rectangular voxel-face clearances, in mm.

Clinical dentin and product applicability are never inferred from tooth occupancy.
"""
import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
SIX = ndi.generate_binary_structure(3, 1)

def cap(t, n, base):
    zz = np.arange(t.shape[0])[:, None, None]
    roof = np.where(t, zz, -1).max(axis=0)
    return t & ~((zz > roof[None] - n) & (zz >= base))

def regional(t, base, prefer_low):
    p0 = cap(t, 4, base)
    target = cap(t, 5, base)
    pmax = cap(t, 7, base)
    need = int(p0.sum() - target.sum())
    candidates = np.argwhere(p0 & ~pmax)
    (yy, xx) = np.nonzero(t.any(axis=0))
    mid = float(np.median(xx))
    zz = np.arange(t.shape[0])[:, None, None]
    roof = np.where(t, zz, -1).max(axis=0)
    depth = roof[candidates[:, 1], candidates[:, 2]] - candidates[:, 0]
    half = candidates[:, 2] <= mid
    priority = (~half if prefer_low else half).astype(int)
    order = np.lexsort((candidates[:, 2], candidates[:, 1], depth, priority))
    selected = candidates[order[:need]]
    if len(selected) != need:
        raise ValueError('Insufficient regional capacity')
    p = p0.copy()
    p[tuple(selected.T)] = False
    assert int((t & ~p).sum()) == int((t & ~target).sum())
    return (p, {'target_removed_voxels': int((t & ~target).sum()), 'extra_voxels': need, 'geometric_halfplane_x_voxel': mid, 'prefer_low_x': bool(prefer_low)})

def s13(t, h, anterior):
    """Source-specific geometric scenario using tooth only, no pulp input.

    Anatomical finish line, crown third and clinical applicability remain unknown.
    Radius/taper are explicit closures, never a fitted manufacturer's parameter.
    """
    occupied = np.flatnonzero(t.any(axis=(1, 2)))
    base = int(occupied[int(np.floor(0.6 * (len(occupied) - 1)))])
    p = cap(t, int(np.ceil(1 / h - 1e-12)), base)
    r = h
    for z in range(base, t.shape[0]):
        dz = (z - base) * h
        a = 1 + r
        depth = a - r + np.sqrt(max(0, r * r - (dz - r) ** 2)) if dz < r else a
        if anterior:
            depth += max(0, dz - r) * np.tan(np.deg2rad(10.0))
        inward = ndi.distance_transform_edt(np.pad(t[z], 1)) * h - h / 2
        p[z] &= inward[1:-1, 1:-1] >= depth
    return (p, {'virtual_finish_z_voxel': base, 'tooth_only_plane': True, 'shoulder_width_mm': 1.0, 'fillet_radius_mm': r, 'axial_above_fillet_mm': 1 + r, 'anterior_total_taper_deg': 20.0 if anterior else None, 'product_applicability': 'UNKNOWN', 'source': 'S13 p1 dated2016-11-14', 'limitations': 'Voxelized geometric realization; true CEJ, cusp/crown-third and sharpness bound UNKNOWN'})

def faces(p, neighbor=None):
    """All exterior faces, or only faces shared with specified removed cells."""
    centers = []
    extents = []
    owners = []
    normals = []
    for axis in range(3):
        for sign in [-1, 1]:
            pad = np.pad(p if neighbor is None else neighbor, 1)
            slices = [slice(1, 1 + n) for n in p.shape]
            slices[axis] = slice(1 + sign, 1 + sign + p.shape[axis])
            near = pad[tuple(slices)]
            selected = p & (~near if neighbor is None else near)
            c = np.argwhere(selected)
            f = c.astype(float)
            f[:, axis] += sign * 0.5
            e = np.full((len(c), 3), 0.5)
            e[:, axis] = 0.0
            normal = np.zeros((len(c), 3), np.int8)
            normal[:, axis] = sign
            centers.append(f)
            extents.append(e)
            owners.append(c)
            normals.append(normal)
    return tuple((np.concatenate(a) for a in [centers, extents, owners, normals]))

def rectangle_distance(q, centers, extents, h, cube=False):
    """Exhaustive point/cube to union-of-rectangles distance, not nearest-centre approximation."""
    if not len(centers):
        return np.full(len(q), np.inf)
    result = np.full(len(q), np.inf)
    for start in range(0, len(q), 128):
        v = q[start:start + 128]
        best = np.full(len(v), np.inf)
        for j in range(0, len(centers), 1024):
            d = np.maximum(np.abs(v[:, None, :] - centers[None, j:j + 1024, :]) - extents[None, j:j + 1024, :] - (0.5 if cube else 0.0), 0.0)
            best = np.minimum(best, np.min(np.sum(d * d, axis=2), axis=1))
        result[start:start + len(v)] = np.sqrt(best) * h
    return result

def measure(t, u, p, h):
    removed = t & ~p
    (f, e, owners, normals) = faces(p, removed)
    if not len(f):
        raise ValueError('Empty cut boundary')
    b = np.zeros(t.shape, bool)
    b[tuple(owners.T)] = True
    q = np.argwhere(u)
    center = ndi.distance_transform_edt(~b, sampling=h)[tuple(q.T)]
    kd = cKDTree(np.argwhere(b) * h).query(q * h, workers=1)[0]
    pt = rectangle_distance(q, f, e, h)
    cube = rectangle_distance(q, f, e, h, cube=True)
    allb = p & ~ndi.binary_erosion(p, structure=SIX, border_value=0)
    whole = ndi.distance_transform_edt(~allb, sampling=h)[tuple(q.T)]
    eps = np.sqrt(3) * h
    low = np.maximum(0, center - eps)
    high = center + eps
    enclosure_error = max(float(np.max(np.maximum(low - pt, 0))), float(np.max(np.maximum(pt - high, 0))), float(np.max(np.maximum(low - cube, 0))), float(np.max(np.maximum(cube - high, 0))))
    exposure = int((u & ~p).sum())
    counts = []
    c = np.argwhere(removed)
    mid = np.median(np.argwhere(t), axis=0)
    for ry in [0, 1]:
        for rx in [0, 1]:
            counts.append(int(np.sum(((c[:, 1] > mid[1]) == bool(ry)) & ((c[:, 2] > mid[2]) == bool(rx)))))
    arg = int(np.argmin(cube))
    summary = {'removed_voxels': int(removed.sum()), 'removed_volume_mm3': float(removed.sum() * h ** 3), 'removed_region_voxels_native_YX_quadrants': counts, 'retained_voxels': int(p.sum()), 'pulp_voxels': int(u.sum()), 'pulp_removed_voxels': exposure, 'cut_boundary_faces': len(f), 'cut_center_min_mm': float(center.min()), 'cut_exact_point_min_mm': float(pt.min()), 'cut_exact_pulp_cube_min_mm': float(cube.min()), 'cut_min_arg_pulp_local_zyx_voxel': q[arg].tolist(), 'whole_prep_center_min_mm': float(whole.min()), 'whole_prep_lower_mm': float(max(0, whole.min() - eps)), 'whole_prep_upper_mm': float(whole.min() + eps), 'cut_lower_mm': float(low.min()), 'cut_upper_mm': float(high.min()), 'epsilon_mm': float(eps), 'EDT_KDTree_max_error_mm': float(np.max(np.abs(center - kd))), 'exact_face_enclosure_error_mm': enclosure_error, 'geometry_gate': 'REJECT_PULP_REMOVED' if exposure else 'DIGITAL_POSITIVE_CUT_CLEARANCE' if cube.min() > 0 else 'DIGITAL_CONTACT_OR_ZERO_CLEARANCE', 'physical_geometry_gate': 'UNKNOWN_UNMEASURED_BOUNDARY_ERRORS', 'remaining_dentin_mm': 'UNKNOWN_NO_DEJ', 'resolution': 'PER_POINT before min aggregation PER_TOOTH', 'timescale': 'HANDOVER'}
    fields = {'pulp_local_zyx_voxel': q, 'center_cut_distance_mm': center, 'exact_point_cut_distance_mm': pt, 'exact_pulp_cube_cut_distance_mm': cube, 'whole_prep_center_distance_mm': whole, 'lower_cut_mm': low, 'upper_cut_mm': high, 'cut_face_centers_zyx_voxel': f, 'cut_face_extents_voxel': e, 'cut_face_normals_zyx': normals}
    return (summary, fields)

def digital_geometry_accepts(u, p):
    return bool(np.all(p[u]))

def validate_port(case, source_sha, canonical_fdi, source_fdi, spacing, expected):
    if case != expected['case']:
        raise ValueError('Foreign case pulp provenance')
    if source_sha != expected['sha256']:
        raise ValueError('Pulp source object hash drift')
    if canonical_fdi != expected['canonical_fdi'] or source_fdi != expected['source_pulpy_fdi']:
        raise ValueError('Patientwise FDI permutation drift')
    if not np.allclose(spacing, [0.3] * 3, rtol=0, atol=1e-12):
        raise ValueError('Invalid physical scale')

def validate_pulp_payload(pulp, expected, crown_sign=1):
    """Bind actual mask content, so merely spoofing case/hash metadata cannot pass."""
    from identity import source_path, sha
    path = source_path(expected['source_path'])
    if sha(path) != expected['sha256']:
        raise ValueError('Source pulp object hash drift')
    with np.load(path) as source:
        reference = source['pulp']
        if crown_sign < 0:
            reference = reference[::-1]
        if pulp.dtype != reference.dtype or pulp.shape != reference.shape or (not np.array_equal(pulp, reference)):
            raise ValueError('Pulp payload differs from exact same-case source mask')

def global_mm(local_zyx, shape, origin, sign, header):
    g = np.asarray(local_zyx, dtype=float).copy()
    if sign < 0:
        g[:, 0] = shape[0] - 1 - g[:, 0]
    g += np.asarray(origin)
    matrix = np.array([float(v) for v in header['TransformMatrix'].split()]).reshape(3, 3)
    offset = np.array([float(v) for v in header['Offset'].split()])
    spacing = np.array([float(v) for v in header['ElementSpacing'].split()])
    return g[:, ::-1] * spacing @ matrix.T + offset
