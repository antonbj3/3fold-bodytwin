from dental_release.paths import expand as _release_expand
import hashlib
import json
import os
import time
from pathlib import Path
import numpy as np
from scipy import ndimage as ndi
ROOT = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get('X22_DATA', _release_expand('@DENTAL_WORK_ROOT@/X22-caries-synthetic')))
CLASSES = ['H0', 'E1', 'E2', 'D1', 'D2', 'D3']

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda : f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def write(p, x):

    def convert(v):
        if isinstance(v, np.ndarray):
            return v.tolist()
        if isinstance(v, np.generic):
            return v.item()
        if isinstance(v, Path):
            return str(v)
        raise TypeError(type(v).__name__)
    Path(p).parent.mkdir(parents=True, exist_ok=True)
    Path(p).write_text(json.dumps(x, indent=2, default=convert, allow_nan=False) + '\n')

def state(stage, last, next_op):
    write(ROOT / 'CURRENT_WORK_STATE.json', {'lane': 'X22-caries-synthetic', 'claim_type': 'capability', 'status': stage, 'updated_at_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), 'latest_gate': last, 'next_operation': next_op, 'data_directory': str(DATA), 'review_state': 'PENDING_INDEPENDENT_REVIEW'})

def tissue(mask, h, shell):
    dt = ndi.distance_transform_edt(mask, sampling=h)
    labels = np.zeros(mask.shape, np.uint8)
    labels[mask] = 2
    labels[mask & (dt <= shell + 1e-09)] = 1
    labels[mask & (dt >= 3.0)] = 3
    return (labels, dt)

def coefficient(labels):
    mu = np.array([0.0, 2.75 * 0.99, 1.44 * 0.99 + 0.5 * 0.24, 0.2683]) / 10
    return mu[labels]

def project(mu, h, angle):
    theta = np.deg2rad(angle)
    mid = (mu.shape[1] - 1) / 2
    (x, z) = np.meshgrid(np.arange(mu.shape[0]), np.arange(mu.shape[2]), indexing='ij')
    out = np.zeros(x.shape, dtype=np.float64)
    for y in range(mu.shape[1]):
        xx = x + np.tan(theta) * (y - mid)
        out += ndi.map_coordinates(mu, [xx, np.full_like(xx, y), z], order=1, mode='constant', cval=0, prefilter=False)
    return out * h / np.cos(theta)

def transmission(tau, h, psf_mm=0.12, scatter=0.1):
    return ndi.gaussian_filter(np.exp(-tau), psf_mm / h, mode='nearest') + scatter

def noisy_image(expected, photons, rng):
    return (rng.poisson(expected * photons) + rng.normal(0, 3, expected.shape)) / photons

def ray_columns(labels, axis):
    a = np.moveaxis(labels, axis, 0)[::-1]
    hard = (a == 1) | (a == 2)
    first = np.argmax(a > 0, axis=0)
    n = a.shape[0]
    idx = np.arange(n)[:, None, None]
    interior = idx >= first[None]
    dej = np.argmax((a != 1) & interior, axis=0)
    end = np.argmax(((a == 0) | (a == 3)) & (idx >= dej[None]), axis=0)
    end = np.where(end > dej, end, n // 2)
    e = dej - first
    d = end - dej
    valid = (a[first, np.arange(a.shape[1])[:, None], np.arange(a.shape[2])[None, :]] == 1) & (e >= 4) & (d >= 8)
    return (a, first, dej, end, e, d, valid)

def choose_site(labels, site, h, origin, shell):
    axis = 0 if site == 'proximal' else 2
    (a, first, dej, end, e, d, valid) = ray_columns(labels, axis)
    if site == 'proximal':
        zs = np.arange(labels.shape[2]) * h + origin[2]
        ztop = np.argwhere(labels > 0)[:, 2].max() * h + origin[2]
        valid &= (zs[None, :] > ztop - 5) & (zs[None, :] < ztop - 2)
        target = np.array([(labels.shape[1] - 1) / 2, (ztop - 3.5 - origin[2]) / h])
    else:
        target = (np.array(a.shape[1:]) - 1) / 2
    q = np.argwhere(valid)
    if not len(q):
        raise ValueError('No eligible path with >=4 enamel and >=8 dentin voxels')
    center = q[np.argmin(np.sum((q - target) ** 2, axis=1))]
    return (axis, center)

def lesion(labels, h, site, center, name):
    axis = 0 if site == 'proximal' else 2
    (a, first, dej, end, e, d, valid) = ray_columns(labels, axis)
    (p, q) = np.indices(a.shape[1:])
    radius = 0.9 if site == 'proximal' else 0.7
    rr = np.sqrt((p - center[0]) ** 2 + (q - center[1]) ** 2) * h / radius
    footprint = (rr <= 1) & valid
    fractions = {'E1': 0.3, 'E2': 0.75, 'D1': 1 / 6, 'D2': 0.5, 'D3': 5 / 6}
    tip = e * fractions[name] if name.startswith('E') else e + d * fractions[name]
    k = np.arange(a.shape[0])[:, None, None]
    depth = k - first[None] + 0.5
    profile = np.clip(1 - rr, 0, 1)[None, :] ** 0.5
    taper = np.maximum(0.35, 1 - 0.65 * np.minimum(depth / np.maximum(e[None], 1), 1))
    if name.startswith('D'):
        taper = np.where(depth > e[None], 0.35 + 0.4 * (depth - e[None]) / np.maximum(d[None], 1), taper)
    support = footprint[None] & (depth > 0) & (depth <= tip[None]) & ((a == 1) | (a == 2)) & (rr[None] <= taper)
    loss = support * profile * (0.3 + 0.7 * np.clip(1 - depth / np.maximum(tip[None], 0.1), 0, 1) ** 0.35)
    loss = np.moveaxis(loss[::-1], 0, axis).astype(np.float32)
    return loss

def verify_truth(labels, loss, h, site):
    axis = 0 if site == 'proximal' else 2
    (a, first, dej, end, e, d, valid) = ray_columns(labels, axis)
    b = np.moveaxis(loss, axis, 0)[::-1]
    supported = b > 0
    k = np.arange(a.shape[0])[:, None, None]
    depth = k - first[None] + 0.5
    ratio = np.where(depth <= e[None], depth / np.maximum(e[None], 1), 1 + (depth - e[None]) / np.maximum(d[None], 1))
    maxratio = float(np.max(ratio[supported])) if supported.any() else 0.0
    bounds = [0.5, 1, 1 + 1 / 3, 1 + 2 / 3, 2 + 1e-06]
    achieved = CLASSES[1 + next((i for (i, bound) in enumerate(bounds) if maxratio <= bound), 4)] if maxratio else 'H0'
    return {'actual_class': achieved, 'max_normalized_depth': maxratio, 'voxels': int(supported.sum()), 'pulp_overlap': int(np.sum(supported & (a == 3))), 'outside_overlap': int(np.sum(supported & (a == 0))), 'max_depth_mm': float(depth[supported].max() * h) if supported.any() else 0.0, 'mineral_fraction_max': float(loss.max())}

def roi_center(center, site, labels, angle):
    axis = 0 if site == 'proximal' else 2
    (a, first, dej, end, e, d, valid) = ray_columns(labels, axis)
    (p, q) = map(int, center)
    i = a.shape[0] - 1 - (first[p, q] + 0.5 * e[p, q])
    if axis == 0:
        (x, y, z) = (i, p, q)
    else:
        (x, y, z) = (p, q, i)
    u = x - np.tan(np.deg2rad(angle)) * (y - (labels.shape[1] - 1) / 2)
    return (float(u), float(z))

def crop(image, center, size=32):
    xy = np.array(center) - size / 2 + 0.5
    (x, z) = np.meshgrid(np.arange(size) + xy[0], np.arange(size) + xy[1], indexing='ij')
    return ndi.map_coordinates(image, [x, z], order=1, mode='nearest')

def features(patch):
    a = -np.log(np.clip(patch - 0.1, 0.001, 2))
    b = ndi.zoom(a, 0.25, order=1, prefilter=False).ravel()
    high = a - ndi.gaussian_filter(a, 2)
    return np.r_[b, a.mean(), a.std(), np.quantile(high, 0.25), np.quantile(high, 0.75)]

def wilson(k, n):
    if not n:
        return [None, None]
    z = 1.96
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    r = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return [float(c - r), float(c + r)]
