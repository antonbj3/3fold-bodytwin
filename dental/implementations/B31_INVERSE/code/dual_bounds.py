from common import *
from parameter_bounds import Response, ratio
from scipy.special import logsumexp
import heapq, time

def dual_cell(ref, opt, a, b, l, h):
    angle = (a + b) / 2
    m0 = (l + h) / 2
    t = np.deg2rad(angle)
    S = ref.S[0] * np.cos(t) + ref.S[1] * np.sin(t)
    (vals, vec) = np.linalg.eigh(S)
    sp = np.maximum(vals[:, -1], 0)
    v = vec[:, :, -1]
    logq = np.full(len(sp), -np.inf)
    ok = sp > 0
    logq[ok] = (m0 - 1) * np.log(sp[ok] / sp.max())
    q = np.exp(logq)
    g = np.array([np.sum(ref.v * q * np.einsum('ei,eij,ej->e', v, s, v)) for s in ref.S])
    candidates = []
    for aa in [a, b]:
        theta = np.deg2rad(aa)
        numer = float(g @ np.array([np.cos(theta), np.sin(theta)]))
        if numer <= 0:
            return -float('inf')
        for m in [l, h]:
            conjugate = m / (m - 1)
            logD = logsumexp(np.log(ref.v[ok]) + conjugate * logq[ok]) / conjugate
            candidates.append(np.log(numer) - logD - opt.logH(opt.principal(aa), m) / m)
    return float(min(candidates))

def bound(ref, opt, maxcells=2048, tol=0.0002):
    queue = []
    counter = 0
    upper = float('inf')
    arg = None
    history = []

    def push(a, b, l, h):
        nonlocal counter, upper, arg
        lower = dual_cell(ref, opt, a, b, l, h)
        for (aa, mm) in [((a + b) / 2, (l + h) / 2), (a, l), (b, h)]:
            val = ratio(ref, opt, aa, mm)
            if val < upper:
                upper = val
                arg = [aa, mm]
        heapq.heappush(queue, (lower, counter, a, b, l, h))
        counter += 1
    push(0.0, 30.0, 2.0, 15.6)
    while counter + 2 <= maxcells and upper - queue[0][0] > tol:
        (lo, _, a, b, l, h) = heapq.heappop(queue)
        if (b - a) / 30 >= (h - l) / 13.6:
            c = (a + b) / 2
            push(a, c, l, h)
            push(c, b, l, h)
        else:
            c = (l + h) / 2
            push(a, b, l, c)
            push(a, b, c, h)
        if counter % 101 == 0:
            history.append([counter, queue[0][0], upper])
    return dict(minimum_force_ratio_lower=float(np.exp(queue[0][0])), minimum_force_ratio_sample_upper=float(np.exp(upper)), log_bound_gap=float(upper - queue[0][0]), cells_evaluated=counter, covering_cells=len(queue), argmin_sample_angle_m=arg, converged=bool(upper - queue[0][0] <= tol), history=history, cover=[list(x) for x in queue])

def midpoint_bound_ok(log_lower, log_ratio):
    return bool(np.isfinite(log_lower) and log_lower <= log_ratio + 1e-10)

def run():
    start = time.perf_counter()
    z = np.load(DATA / 'B_FIELDS.npz')
    old = read(ROOT / 'raw/B_BOUNDS.json')
    rows = {}
    for s in ['8600', '18000', 'rigid']:
        ref = Response(z['reference__' + s + '__stress_basis'], z['reference__' + s + '__vol'])
        opt = Response(z['optimized__' + s + '__stress_basis'], z['optimized__' + s + '__vol'])
        r = bound(ref, opt)
        cover = r.pop('cover')
        save(ROOT / f'raw/D_COVER_{s}.json', cover)
        errors = []
        record_checks = []
        for (lo, _, a, b, l, h) in cover:
            actual = ratio(ref, opt, (a + b) / 2, (l + h) / 2)
            errors.append(max(0.0, lo - actual))
            record_checks.append(midpoint_bound_ok(lo, actual))
        r.update(max_cover_midpoint_violation=float(max(errors)), midpoint_records_pass=all(record_checks), wrong_bound_rejected=not midpoint_bound_ok(actual + 0.01, actual), old_log_gap=old['contrasts']['optimized'][s]['continuous']['log_bound_gap'])
        assert r['midpoint_records_pass'] and r['wrong_bound_rejected']
        rows[s] = r
        print(s, json.dumps(r), flush=True)
    out = dict(claim_type='capability', round='D', supports=rows, seconds=time.perf_counter() - start, formal_floating_certificate='MISSING', FE_continuum_enclosure='MISSING', scope='Analytic real-arithmetic continuous m/angle enclosure for stored discrete tensors; finite support set', external_referent=read(ROOT / 'PREREG_D.json')['external_referent'])
    save(ROOT / 'raw/D.json', out)
    freeze(ROOT / 'FROZEN_PREDICTIONS_D.json', dict(prereg_sha256=sha(ROOT / 'PREREG_D.json'), field_sha256=sha(DATA / 'B_FIELDS.npz'), supports=rows, physical_measurement='NOT_RUN', absolute_force_N=None))
    state('D_COMPLETE', dict(bounds={k: v['minimum_force_ratio_lower'] for (k, v) in rows.items()}), 'Assemble independently checkable demo and physical no-go/next measurement contract')
if __name__ == '__main__':
    run()
