from common import *
from fractions import Fraction as F
os.environ['NUMBA_CACHE_DIR'] = str(DATA / 'numba')
saved = sys.modules.get('contact')
path = list(sys.path)
try:
    X = module('x21_contact_readonly', ROOT.parent / 'LANE_X21_CONTACT_MAP/code/contact.py')
    sys.modules['contact'] = X
    C = module('x21_certificates_readonly', ROOT.parent / 'LANE_X21_CONTACT_MAP/code/certificates.py')
finally:
    if saved is not None:
        sys.modules['contact'] = saved
    else:
        sys.modules.pop('contact', None)
    sys.path[:] = path
M = module('tandlast_readonly', TAND / 'code/model.py')

def pair_intervals(U, L, ul, ll, control=False):
    (uc, ok) = X.affine(U)
    U = np.ascontiguousarray(U[ok])
    ul = np.ascontiguousarray(ul[ok])
    uc = uc[ok]
    (lc, ok2) = X.affine(L)
    L = np.ascontiguousarray(L[ok2])
    ll = np.ascontiguousarray(ll[ok2])
    lc = lc[ok2]
    ub = X.boxes(U)
    lb = X.boxes(L)
    origin = np.minimum(ub[:, :2].min(0), lb[:, :2].min(0)) - 0.8
    hi = np.maximum(ub[:, 3:5].max(0), lb[:, 3:5].max(0))
    shape = np.ceil((hi - origin) / 0.8).astype(np.int64) + 2
    (ptr, ids) = X.buckets(lb, origin, 0.8, shape)
    (best, wi, xy, checked, tested) = X.extrema(U, L, ul, ll, uc, lc, ub, lb, origin, 0.8, shape, ptr, ids, True)
    (lower, _) = C.lower_bounds(ub, lb, ul, ll, origin, 0.8, shape, ptr, ids)
    rows = []
    for (a, b) in np.argwhere(np.isfinite(best)):
        (i, j) = wi[a, b]
        (hiQ, point) = C.rational_min(U[i], L[j])
        if hiQ is None:
            raise ArithmeticError('Floating witness not in exact projected intersection')
        lo = float(np.nextafter(lower[a, b], -np.inf))
        high = float(np.nextafter(float(hiQ), np.inf))
        if not lo <= best[a, b] + 1e-07 or not best[a, b] <= high + 1e-07 or lo > high:
            raise ArithmeticError('Gap enclosure rejected')
        rows.append(dict(upper=int(a), lower=int(b), gap_bounds_mm=[lo, high], nominal_gap_mm=float(best[a, b]), witness_xy_mm=[float(x) for x in point], witness_upper_triangle=U[i], witness_lower_triangle=L[j], exact_upper=str(hiQ)))
    check = None
    if control and rows:
        (a, b) = (rows[0]['upper'], rows[0]['lower'])
        uu = U[ul == a]
        lt = L[ll == b]
        uub = X.boxes(uu)
        llb = X.boxes(lt)
        (cu, _) = X.affine(uu)
        (cl, _) = X.affine(lt)
        (p0, k0) = X.buckets(llb, origin, 0.8, shape)
        (truth, _, _, _, _) = X.extrema(uu, lt, np.full(len(uu), a, np.int32), np.full(len(lt), b, np.int32), cu, cl, uub, llb, origin, 0.8, shape, p0, k0, False)
        err = abs(float(truth[a, b]) - best[a, b])
        lo2 = C.independent_pair_lb(uub, llb)
        check = dict(pair=[a, b], unpruned_error_mm=float(err), lower_bound_error_mm=abs(lo2 - rows[0]['gap_bounds_mm'][0]), injected_plus1_rejected=abs(float(truth[a, b]) - (best[a, b] + 1)) > 1e-07)
    return (rows, dict(aabb_pairs=int(checked), intersections=int(tested), upper_degenerate_excluded=int((~ok).sum()), lower_degenerate_excluded=int((~ok2).sum()), control=check))

def edges(rows, delta=0.0):
    return [(r['upper'], r['lower'], F(r['gap_bounds_mm'][0]) - F(str(delta)), F(r['gap_bounds_mm'][1]) + F(str(delta))) for r in rows]

def classes(rows, delta=0.0):
    e = edges(rows, delta)
    out = {}
    for t in M.FDI:
        a = M.classify(e, [t])
        out[str(t)] = dict(classification=a['classification'], share_hull_pp=list(map(float, a['share_hull_pp'])), force_hull_N=list(map(float, a['force_hull_N'])), infimum_attained=a['infimum_attained'], necessity_margin_mm=float(a['necessity_margin_mm']) if a['necessity_margin_mm'] is not None else None)
    return out
