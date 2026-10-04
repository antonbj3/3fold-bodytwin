import numpy as np
from .bounds import support

def polygon_integral(poly, xy, values):
    if len(poly) < 3:
        return (0.0, 0.0)
    total = 0.0
    integ = 0.0
    for j in range(1, len(poly) - 1):
        b = poly[[0, j, j + 1]]
        v = b @ xy
        area = abs((v[1, 0] - v[0, 0]) * (v[2, 1] - v[0, 1]) - (v[1, 1] - v[0, 1]) * (v[2, 0] - v[0, 0])) / 2
        total += area
        integ += area * float(np.mean(b @ values))
    return (total, integ)

def reachable(q, t, lo, hi, r, with_polygons=False):
    xy = t['xy']
    ff = support(t)
    (R, ar, _, _) = q.polygons(xy, ff, r)
    P = {}
    D = {}
    unreachable = {}
    ad = 0.0
    ap = 0.0
    over_possible = 0.0
    over_definite = 0.0
    addition_int = 0.0
    addition_max = 0.0
    for (fi, ids) in enumerate(ff):
        vl = lo[ids]
        vh = hi[ids]
        possible = q.Q.clip(q.Q.clip(np.eye(3), 0.1 - vl), vh)
        definite = q.Q.clip(q.Q.clip(np.eye(3), vl), 0.1 - vh)
        P[fi] = possible
        D[fi] = definite
        (pa, _) = q.measure(possible, xy[ids])
        (da, _) = q.measure(definite, xy[ids])
        ap += pa
        ad += da
        if fi in R:
            rp = R[fi]
            op = q.Q.clip(q.Q.clip(rp, 0.1 - vl), vh)
            (oparea, _) = q.measure(op, xy[ids])
            over_possible += oparea
            od = q.Q.clip(q.Q.clip(rp, vl), 0.1 - vh)
            (odarea, _) = q.measure(od, xy[ids])
            over_definite += odarea
            addpoly = q.Q.clip(rp, vl - 0.1)
            (aa, ii) = polygon_integral(addpoly, xy[ids], vl - 0.1)
            addition_int += ii
            if aa > 1e-12:
                addition_max = max(addition_max, float((addpoly @ (vl - 0.1)).max()))
            if with_polygons and ar[list(R).index(fi)] - oparea > 1e-12:
                unreachable[str(fi)] = dict(grid_vertex_ids=ids, source_polygon_barycentric=rp, source_polygon_xy_mm=rp @ xy[ids], gap_lower_mm=vl, gap_upper_mm=vh, minimum_positive_addition_polygon_barycentric=addpoly)
    ra = float(ar.sum())
    missing = max(0.0, ra - over_possible)
    spurious = max(0.0, ad - over_definite)
    return dict(unavoidable_missing_contact_mm2=missing, unavoidable_spurious_contact_mm2=spurious, unavoidable_symdiff_floor_mm2=missing + spurious, reference_area_mm2=ra, possible_predicted_area_mm2=ap, definite_predicted_area_mm2=ad, minimum_added_height_integral_mm3=addition_int, minimum_added_height_max_mm=addition_max, polygons=unreachable if with_polygons else None)

def direct_partition_floor(q, t, lo, hi, r):
    """Independent Boolean arrangement evaluation, rather than the set-area difference formula."""
    ff = support(t)
    total = 0.0
    cells = 0
    for ids in ff:
        vl = lo[ids]
        vh = hi[ids]
        vr = r[ids]
        polys = [np.eye(3)]
        for values in [vl, vl - 0.1, vh, vh - 0.1, vr, vr - 0.1]:
            out = []
            for p in polys:
                vv = p @ values
                if vv.min() >= 0 or vv.max() <= 0:
                    out.append(p)
                else:
                    for pp in [q.Q.clip(p, values), q.Q.clip(p, -values)]:
                        if len(pp) >= 3:
                            out.append(pp)
            polys = out
        for p in polys:
            (ar, _) = q.measure(p, t['xy'][ids])
            b = p.mean(0)
            l = float(b @ vl)
            h = float(b @ vh)
            rr = float(b @ vr)
            target = 0 <= rr <= 0.1
            match_possible = l <= 0.1 and h >= 0 if target else l < 0 or h > 0.1
            if not match_possible:
                total += ar
            cells += 1
    return (float(total), cells)

def box_lower(q, t, lo, hi, r, delta):
    ff = support(t)
    xy = t['xy']
    (Rd, ar, _, _) = q.polygons(xy, ff, r, delta, 0.1 - delta)
    (Ru, _, _, _) = q.polygons(xy, ff, r, -delta, 0.1 + delta)
    possible = {}
    definite = {}
    da = 0.0
    overR = 0.0
    overD = 0.0
    for (fi, ids) in enumerate(ff):
        p = q.Q.clip(q.Q.clip(np.eye(3), 0.1 + delta - lo[ids]), hi[ids] + delta)
        d = q.Q.clip(q.Q.clip(np.eye(3), lo[ids] - delta), 0.1 - delta - hi[ids])
        possible[fi] = p
        definite[fi] = d
        (a, _) = q.measure(d, xy[ids])
        da += a
        if fi in Rd:
            z = q.Q.clip(q.Q.clip(Rd[fi], 0.1 + delta - lo[ids]), hi[ids] + delta)
            (a, _) = q.measure(z, xy[ids])
            overR += a
        z = q.Q.clip(q.Q.clip(d, r[ids] + delta), 0.1 + delta - r[ids])
        (a, _) = q.measure(z, xy[ids])
        overD += a
    return max(0.0, float(ar.sum() - overR)) + max(0.0, da - overD)

def sufficiency(q):
    left = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    xy = np.r_[left, left + [1.0, 0.0]]
    ff = np.array([[0, 1, 2], [3, 4, 5]])
    a = np.array([0.05, 0.05, 0.05, 0.2, 0.2, 0.2])
    b = a[[3, 4, 5, 0, 1, 2]]
    aa = q.compare(xy, ff, a, a)
    bb = q.compare(xy, ff, b, a)
    da = a - 0.05
    db = b - 0.05
    Ja = 0.5 * (da[:3].mean() + da[3:].mean())
    Jb = 0.5 * (db[:3].mean() + db[3:].mean())
    sa = np.array([aa['predicted']['area_mm2'], aa['predicted']['count'], np.quantile(a, 0.95), np.mean(np.sort(a)), Ja, 0.8 - da.max()])
    sb = np.array([bb['predicted']['area_mm2'], bb['predicted']['count'], np.quantile(b, 0.95), np.mean(np.sort(b)), Jb, 0.8 - db.max()])
    identity = float(abs(sa - sb).max())
    return dict(summary_names=['contact_area_mm2', 'component_count', 'gap_p95_mm', 'sorted_gap_mean_mm', 'area_integral_abs_displacement_mm3', 'wall_Lipschitz_lower_mm'], state_A_summary=sa, state_B_summary=sb, summary_identity_error=identity, bit_identical=bool(np.array_equal(sa, sb)), mask_translation_mm=1.0, contact_area_change_mm2=bb['predicted']['area_mm2'] - aa['predicted']['area_mm2'], overlap_A_mm2=aa['intersection_mm2'], overlap_B_mm2=bb['intersection_mm2'], downstream_symdiff_A_mm2=aa['symdiff_mm2'], downstream_symdiff_B_mm2=bb['symdiff_mm2'], downstream_difference_mm2=bb['symdiff_mm2'] - aa['symdiff_mm2'], translated_mask_rejected=bb['intersection_mm2'] < aa['reference']['area_mm2'] - 1e-07, untranslated_mask_accepted=aa['intersection_mm2'] >= aa['reference']['area_mm2'] - 1e-07, minimal_extension='For fixed symmetric-difference query: reference-overlap area. For a repair action or changed reference/pose: spatial mask support and case/frame are required.', resolution=dict(areas='PER_SURFACE_REGION', gap_and_wall='PER_POINT', displacement_integral='PER_TOOTH'), timescale='SIMULTANEOUS', external_referent=dict(kind='our_own_fixture', locator='code/reachability.py:sufficiency', compared_quantity='Exact-summary sufficiency counterexample and1mm mask translation operator', refutes_us=True), force_summary='UNKNOWN: no invented Newton value')
