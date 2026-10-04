"""Independent full triangle enumeration, no KD-tree and no inherited height code."""
import numpy as np

def direct_height(tri, xy):
    (a, b, c) = (tri[:, 0], tri[:, 1], tri[:, 2])
    B = b[:, :2] - a[:, :2]
    C = c[:, :2] - a[:, :2]
    det = B[:, 0] * C[:, 1] - B[:, 1] * C[:, 0]
    nondeg = np.abs(det) > 1e-09
    out = []
    for p in xy:
        d = p - a[:, :2]
        v = np.divide(d[:, 0] * C[:, 1] - d[:, 1] * C[:, 0], det, out=np.zeros(len(det)), where=nondeg)
        w = np.divide(B[:, 0] * d[:, 1] - B[:, 1] * d[:, 0], det, out=np.zeros(len(det)), where=nondeg)
        hit = nondeg & (v >= -1e-09) & (w >= -1e-09) & (v + w <= 1 + 1e-09)
        z = a[:, 2] + v * (b[:, 2] - a[:, 2]) + w * (c[:, 2] - a[:, 2])
        out.append(float(z[hit].max()) if hit.any() else np.nan)
    return np.array(out)

def run(bundle):
    rows = []
    policy = bundle.json('PREREG_R1.json')['metrics']
    for p in sorted(bundle.files):
        if not p.startswith('payload/private/audit/') or not p.endswith('.npz'):
            continue
        a = bundle.npz(p)
        tri = (a['original_world_triangles'] - a['base']) @ a['R']
        truth = direct_height(tri, a['xy'])
        ref = a['reference']
        ok = np.isfinite(truth) & np.isfinite(ref)
        error = float(np.max(np.abs(truth[ok] - ref[ok]))) if ok.any() else None
        mask_match = np.array_equal(np.isfinite(truth), np.isfinite(ref))
        gate = mask_match and error is not None and (error <= policy['reference_recompute_tolerance_mm'])
        shift = direct_height(tri + [0, 0, 0.2], a['xy'])
        bad = float(np.max(np.abs(shift[ok] - ref[ok]))) > policy['reference_recompute_tolerance_mm']
        case = p.rsplit('/', 1)[-1][:-4]
        rname = 'payload/private/references/' + case + '.npz'
        live = bundle.npz(rname)['molar_crown']
        binding = np.array_equal(ref, live, equal_nan=True)
        rows.append(dict(case_key=case, points=int(ok.sum()), max_error_mm=error, finite_mask_match=mask_match, reference_binding=binding, gate=gate and binding, injected_surface_shift_rejected=bad, resolution='PER_POINT'))
    return dict(rows=rows, source_cases=len(rows), all_pass=bool(rows) and all((r['gate'] and r['injected_surface_shift_rejected'] for r in rows)), external_referent=dict(kind='published_dataset', locator='https://ditto.ing.unimore.it/bite2text/', compared_quantity='Original STL pointwise occlusal height at locked site coordinates, recomputed by independent full triangle enumeration', refutes_us=True))
