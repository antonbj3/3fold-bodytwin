"""Source voxel geometry/intensity ports. No HU, mineral density or clinical inference."""
import numpy as np
from scipy import ndimage as ndi
ENVELOPE = [1, 3, 4, 10] + list(range(31, 39)) + list(range(41, 49))
TEETH = list(range(11, 19)) + list(range(21, 29)) + list(range(31, 39)) + list(range(41, 49))

def ray_grid(shape, spacing, origin, direction, lower, upper):
    """Partition a ray at every voxel-box face, including volume limits."""
    ts = [lower, upper]
    for k in range(3):
        if abs(direction[k]) > 1e-14:
            q = ((np.arange(shape[k] + 1) - 0.5) * spacing[k] - origin[k]) / direction[k]
            ts.extend(q[(q > lower) & (q < upper)].tolist())
    t = np.unique(ts)
    mid = (t[:-1] + t[1:]) / 2
    ix = np.floor((origin + mid[:, None] * direction) / spacing + 0.5).astype(int)
    good = np.all((ix >= 0) & (ix < np.asarray(shape)), axis=1)
    return (t[:-1][good], t[1:][good], ix[good])

def connected_run(lo, hi, labels, anchor=0.0):
    mask = np.isin(labels, ENVELOPE)
    hit = np.flatnonzero(mask & (lo <= anchor) & (hi >= anchor))
    if not len(hit):
        return None
    j = int(hit[0])
    a = j
    b = j
    while a and mask[a - 1] and (abs(hi[a - 1] - lo[a]) < 1e-09):
        a -= 1
    while b + 1 < len(mask) and mask[b + 1] and (abs(hi[b] - lo[b + 1]) < 1e-09):
        b += 1
    return (float(lo[a]), float(hi[b]))

def slab_control(lab, sp, o, u, lower, upper, anchor=0.0):
    """Independent enumeration of ALL voxel boxes in the ray bounding rectangle."""
    ends = np.array([o + lower * u, o + upper * u])
    il = np.maximum(0, np.floor(ends.min(0) / sp - 0.5).astype(int))
    ih = np.minimum(lab.shape, np.ceil(ends.max(0) / sp + 0.5).astype(int) + 1)
    if np.any(ih <= il):
        return None
    ix = np.indices(tuple(ih - il)).reshape(3, -1).T + il
    centers = ix * sp
    left = np.full(len(ix), lower)
    right = np.full(len(ix), upper)
    for k in range(3):
        if abs(u[k]) < 1e-14:
            bad = abs(centers[:, k] - o[k]) > sp[k] / 2
            left[bad] = 1
            right[bad] = 0
        else:
            a = (centers[:, k] - sp[k] / 2 - o[k]) / u[k]
            b = (centers[:, k] + sp[k] / 2 - o[k]) / u[k]
            left = np.maximum(left, np.minimum(a, b))
            right = np.minimum(right, np.maximum(a, b))
    good = left < right - 1e-12
    left = left[good]
    right = right[good]
    ix = ix[good]
    order = np.argsort(left)
    return connected_run(left[order], right[order], lab[tuple(ix[order].T)], anchor)

def lateral_axis(lab, sp, e, a):
    b = np.array([0.0, 1.0, 0.0])
    b -= a * np.dot(a, b)
    if np.linalg.norm(b) < 0.1:
        b = np.array([0.0, 0.0, 1.0])
        b -= a * np.dot(a, b)
    b /= np.linalg.norm(b)
    sl = tuple((slice(max(0, int(v / sp[k]) - 50), min(lab.shape[k], int(v / sp[k]) + 51)) for (k, v) in enumerate(e)))
    teeth = np.argwhere(np.isin(lab[sl], TEETH))
    off = np.array([s.start for s in sl])
    if len(teeth):
        c = (teeth + off).mean(0) * sp
        if np.dot(e - c, b) < 0:
            b = -b
    return b

def arch_centers(lab, sp):
    objects = ndi.find_objects(lab, max_label=48)
    out = []
    for k in list(range(31, 39)) + list(range(41, 49)):
        sl = objects[k - 1]
        if sl is None:
            continue
        pts = np.argwhere(lab[sl] == k)
        if len(pts):
            out.append((k, (pts.mean(0) + np.array([s.start for s in sl])) * sp))
    return out

def arch_normal(centers, e, a, lab, sp):
    if len(centers) < 3:
        return (lateral_axis(lab, sp, e, a), {'status': 'FALLBACK_TOO_FEW_TEETH', 'n_teeth': len(centers)})
    near = sorted(centers, key=lambda v: float(np.linalg.norm(v[1] - e)))[:4]
    xyz = np.array([v[1] for v in near])
    x = xyz - xyz.mean(0)
    (eig, V) = np.linalg.eigh(x.T @ x)
    t = V[:, -1]
    t -= a * np.dot(a, t)
    if np.linalg.norm(t) < 0.1:
        return (lateral_axis(lab, sp, e, a), {'status': 'FALLBACK_DEGENERATE_TANGENT', 'n_teeth': len(centers)})
    t /= np.linalg.norm(t)
    b = np.cross(a, t)
    b /= np.linalg.norm(b)
    c = np.mean([v[1] for v in centers], axis=0)
    if np.dot(e - c, b) < 0:
        b = -b
    return (b, {'status': 'LOCAL_ARCH_PROXY', 'source_fdis': [v[0] for v in near], 'tangent_zyx': t.tolist(), 'normal_axis_dot': float(np.dot(b, a)), 'normal_tangent_dot': float(np.dot(b, t)), 'tangent_residual_rms_mm': float(np.sqrt(np.mean(np.sum((x - (x @ t)[:, None] * t) ** 2, axis=1)))), 'true_anatomical_buccolingual_plane': 'UNKNOWN'})

def sample_nearest(arr, sp, pts):
    ix = np.floor(pts / sp + 0.5).astype(int)
    good = np.all((ix >= 0) & (ix < np.asarray(arr.shape)), axis=1)
    vals = np.zeros(len(pts), float)
    vals[good] = arr[tuple(ix[good].T)]
    return (vals, good, ix)

def shell_profile(img, lab, sp, o, inward, cap):
    s = np.arange(0, 5.0001, 0.15)
    pts = o + s[:, None] * inward
    (labels, ok, ix) = sample_nearest(lab, sp, pts)
    (nearest, okg, _) = sample_nearest(img, sp, pts)
    gray = ndi.map_coordinates(img, (pts / sp).T, order=1, mode='constant', cval=0, prefilter=False)
    reason = None
    if not np.all(ok):
        reason = 'PROFILE_OUTSIDE_FOV'
    elif np.any(np.isin(labels, TEETH + [8, 9, 10])):
        reason = 'PROFILE_CROSSES_TOOTH_OR_METAL'
    elif np.any(gray >= cap - 1):
        reason = 'SATURATED_PROFILE'
    core = float(np.median(gray[s >= 3]))
    bright = float(np.quantile(gray[s <= 3], 0.9))
    contrast = bright - core
    if contrast <= 0:
        reason = reason or 'NO_POSITIVE_SHELL_CONTRAST'
    thickness = []
    for tau in [0.4, 0.5, 0.6]:
        high = gray >= core + tau * contrast
        first = np.flatnonzero(high & (s <= sp.max()))
        if not len(first):
            thickness.append(None)
            continue
        j = int(first[0])
        k = j
        while k + 1 < len(s) and high[k + 1]:
            k += 1
        thickness.append(None if k + 1 == len(s) else float(s[k] + 0.15))
    if any((v is None for v in thickness)):
        reason = reason or 'NO_RESOLVED_INNER_CROSSING'
    return {'gray_nearest': nearest, 'gray_interp': gray, 'depth_mm': s, 'labels': labels.astype(np.uint8), 'ix': ix, 'apparent_layer_mm_tau_04_05_06': thickness, 'rejection_reason': reason, 'contrast_gray': contrast, 'anatomical_cortex_mm': None, 'anatomical_accuracy': 'UNKNOWN_NO_MATCHED_INNER_BOUNDARY'}

def core_samples(img, lab, sp, e, a):
    pts = np.array([e, e + 8 * a])
    low = np.maximum(0, np.floor((pts.min(0) - 5) / sp).astype(int))
    high = np.minimum(lab.shape, np.ceil((pts.max(0) + 5) / sp).astype(int) + 1)
    sl = tuple((slice(int(x), int(y)) for (x, y) in zip(low, high)))
    L = lab[sl]
    ix = np.indices(L.shape).reshape(3, -1).T
    xyz = (ix + low) * sp
    delta = xyz - e
    z = delta @ a
    rad = np.linalg.norm(delta - z[:, None] * a, axis=1)
    envelope = np.isin(L, ENVELOPE)
    dist = ndi.distance_transform_edt(envelope, sampling=sp).ravel()
    toothdist = ndi.distance_transform_edt(~np.isin(L, TEETH + [8, 9, 10]), sampling=sp).ravel()
    mask = (L.ravel() == 1) & (z >= 1) & (z <= 8) & (rad >= 2.2) & (rad <= 4) & (dist >= 1) & (toothdist >= 1)
    ids = ix[mask] + low
    gray = img[tuple(ids.T)]
    return (ids, gray)
