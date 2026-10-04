"""Conditional airway geometry, in mm. Expert labels are measured input.
No tissue parameter in this module is an externally measured patient parameter.
"""
import numpy as np
from scipy import ndimage

def profile(mask, spacing, z_offset=0, margin=5.0):
    (dz, dy, dx) = spacing
    area = mask.sum(axis=(1, 2), dtype=np.int64) * dx * dy
    lateral = mask.any(axis=1).sum(axis=1) * dx
    zs = np.arange(len(area)) * dz + z_offset * dz
    nonzero = np.where(area > 0)[0]
    eligible = (area >= 1) & (zs >= zs[nonzero[0]] + margin) & (zs <= zs[nonzero[-1]] - margin)
    if not eligible.any():
        raise ValueError('No interior slices at fixed margin')
    return dict(z_mm=zs, area_mm2=area, lateral_support_mm=lateral, eligible=eligible, volume_mm3=float(area.sum() * dz))

def pose(points_zyx, pivot_zyx, incisor_zyx, advancement, angle_deg, ap_sign, si_sign):
    """Positive angle rotates anterior point inferiorly. Solve incisor AP command."""
    p = np.asarray(points_zyx, float)
    c = np.asarray(pivot_zyx, float)
    inc = np.asarray(incisor_zyx, float)
    a = (p[..., 1] - c[1]) * ap_sign
    z = (p[..., 0] - c[0]) * si_sign
    ia = (inc[1] - c[1]) * ap_sign
    iz = (inc[0] - c[0]) * si_sign
    th = np.deg2rad(angle_deg)
    co = np.cos(th)
    sn = np.sin(th)
    translation = advancement - (co * ia + sn * iz - ia)
    delta_a = co * a + sn * z + translation - a
    delta_z = -sn * a + co * z - z
    out = p.copy()
    out[..., 1] += delta_a * ap_sign
    out[..., 0] += delta_z * si_sign
    return (out, float(translation))

def landmarks(label, spacing):
    sp = np.array(spacing)

    def cen(ids):
        pts = np.argwhere(np.isin(label, ids))
        if len(pts) == 0:
            raise ValueError('required tooth landmark missing: ' + str(ids))
        return pts.mean(axis=0) * sp
    inc = cen([31, 41])
    uinc = cen([11, 21])
    mol = cen([36, 37, 46, 47])
    ap = 1 if inc[1] > mol[1] else -1
    si = 1 if uinc[0] > inc[0] else -1
    jaw = np.argwhere(label == 1).astype(np.float32) * sp
    if len(jaw) == 0:
        raise ValueError('mandible missing')
    superior = jaw[:, 0] * si
    posterior = jaw[:, 1] * ap
    subset = jaw[(superior >= np.quantile(superior, 0.85)) & (posterior <= np.quantile(posterior, 0.4))]
    split = (float(jaw[:, 2].min()) + float(jaw[:, 2].max())) / 2
    sides = [subset[subset[:, 2] < split], subset[subset[:, 2] >= split]]
    if any((len(s) == 0 for s in sides)):
        raise ValueError('bilateral condyle proxy incomplete')
    condyles = np.array([s.mean(axis=0) for s in sides])
    if abs(condyles[1, 2] - condyles[0, 2]) < 30:
        raise ValueError('bilateral proxy separation below 30 mm')
    return dict(incisor_zyx_mm=inc.tolist(), upper_incisor_zyx_mm=uinc.tolist(), condyle_proxies_zyx_mm=condyles.tolist(), pivot_zyx_mm=condyles.mean(axis=0).tolist(), ap_sign=ap, si_sign=si, proxy_status='CONSTITUTIVE_CLOSURE', added_incisor_vertical_opening_at_6deg_mm=float(-(pose(inc, condyles.mean(axis=0), inc, 6, 6, ap, si)[0] - inc)[0] * si))

def crop_pharynx(label):
    mask = label == 7
    (lab, n) = ndimage.label(mask, structure=ndimage.generate_binary_structure(3, 1))
    count = np.bincount(lab.ravel())
    count[0] = 0
    if not n:
        raise ValueError('pharynx absent')
    chosen = int(count.argmax())
    objects = ndimage.find_objects(lab)
    sl = objects[chosen - 1]
    crop = lab[sl] == chosen
    return (crop, tuple((s.start for s in sl)), float(count[chosen] / mask.sum()), n)

def scenarios(pr, lm, spacing, crop_offset, mask, advancements, rotations, transfers):
    pts = []
    for i in range(len(mask)):
        yx = np.argwhere(mask[i])
        yz = yx.mean(axis=0) if len(yx) else np.zeros(2)
        pts.append([pr['z_mm'][i], (yz[0] + crop_offset[1]) * spacing[1], (yz[1] + crop_offset[2]) * spacing[2]])
    pts = np.array(pts)
    rows = []
    for m in advancements:
        for theta in rotations:
            (moved, t) = pose(pts, lm['pivot_zyx_mm'], lm['incisor_zyx_mm'], m, theta, lm['ap_sign'], lm['si_sign'])
            d = (moved[:, 1] - pts[:, 1]) * lm['ap_sign']
            for k in transfers:
                a = pr['area_mm2'] + pr['lateral_support_mm'] * k * d
                eligible = pr['eligible']
                ids = np.where(eligible)[0]
                ix = int(ids[np.argmin(a[eligible])])
                vol = float(a.sum() * spacing[0])
                factor = vol / pr['volume_mm3']
                rows.append(dict(advancement_mm=m, rotation_deg=theta, transfer=k, minimum_mm2=float(a[ix]), minimum_z_mm=float(pr['z_mm'][ix]), volume_mm3=vol, volume_control_minimum_mm2=float(pr['area_mm2'][eligible].min() * factor), anterior_condyle_translation_mm=t, negative_interior_area=bool((a[eligible] < 0).any())))
    return rows

def gate(value, tolerance):
    return bool(np.isfinite(value) and value <= tolerance)
