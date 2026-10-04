"""Reference-owned quality operators. No participant objective enters these scores."""
import numpy as np

def clip(poly, values):
    """Clip a barycentric polygon by sum(bary*values)>=0."""
    if len(poly) == 0:
        return poly
    out = []
    for (a, b) in zip(poly, np.roll(poly, -1, axis=0)):
        (va, vb) = (a @ values, b @ values)
        if va >= 0:
            out.append(a)
        if (va >= 0) != (vb >= 0):
            out.append(a + (b - a) * (va / (va - vb)))
    return np.array(out).reshape(-1, 3)

def area(poly, xy):
    if len(poly) < 3:
        return 0.0
    q = poly @ xy
    return float(abs(np.sum(q[:, 0] * np.roll(q[:, 1], -1) - q[:, 1] * np.roll(q[:, 0], -1))) / 2)

def band_polygon(g, band):
    return clip(clip(np.eye(3), g), band - g)

def contact_map(xy, faces, pred_gap, ref_gap, band=0.1):
    """Exact polygon integration on specified PL gap fields, float arithmetic."""
    predicted = reference = intersection = negative = support = 0.0
    for ids in faces:
        (pg, rg) = (pred_gap[ids], ref_gap[ids])
        if not np.isfinite(pg).all() or not np.isfinite(rg).all():
            continue
        q = xy[ids]
        p = band_polygon(pg, band)
        r = band_polygon(rg, band)
        support += area(np.eye(3), q)
        predicted += area(p, q)
        reference += area(r, q)
        intersection += area(clip(clip(p, rg), band - rg), q)
        negative += area(clip(np.eye(3), -pg), q) if np.min(pg) < 0 else 0.0
    sym = max(0.0, predicted + reference - 2 * intersection)
    union = predicted + reference - intersection
    if support == 0:
        return dict(contact_symdiff_mm2=None, contact_area_error_mm2=None, reference_contact_mm2=None, predicted_contact_mm2=None, contact_iou=None, negative_gap_area_mm2=None, contact_support_mm2=0.0, empty_contact_union=None)
    return dict(contact_symdiff_mm2=sym, contact_area_error_mm2=abs(predicted - reference), reference_contact_mm2=reference, predicted_contact_mm2=predicted, contact_iou=intersection / union if union > 1e-12 else None, negative_gap_area_mm2=negative, contact_support_mm2=support, empty_contact_union=union <= 1e-12)

def roof_metrics(t, z, inner, ref):
    w = np.asarray(t['weights'])
    ok = np.isfinite(ref)
    den = w[ok].sum()
    with np.errstate(over='raise', invalid='raise'):
        if not np.isfinite(z).all() or not np.isfinite(inner).all():
            raise ValueError('nonfinite prediction')
        e = z[ok] - ref[ok]
        rmse = float(np.sqrt(np.sum(w[ok] * e * e) / den)) if den > 0 else None
        gap = inner - t['preparation_z']
        out = dict(anatomy_rmse_mm=rmse, anatomy_max_mm=float(np.max(abs(e))) if len(e) else None, reference_coverage=float(den / w.sum()), nominal_gap_mean_mm=float(w @ gap / w.sum()), nominal_gap_min_mm=float(gap.min()), nominal_gap_max_mm=float(gap.max()), nominal_gap_range_mm=float(np.ptp(gap)))
        out.update(contact_map(t['xy'], t['faces'], t['ceiling'] - z, t['ceiling'] - ref))
    return out

def nondominated(values, tolerances):
    x = np.asarray(values, float)
    tol = np.asarray(tolerances, float)
    if not np.isfinite(x).all():
        raise ValueError('Unknown coordinates cannot be ranked')
    return np.array([not any((np.all(y <= v + tol) and np.any(y < v - tol) for (j, y) in enumerate(x) if i != j)) for (i, v) in enumerate(x)])
