"""Analytic continuous m/angle enclosure of a discrete FE response.

Weyl bounds tensile eigenvalues for all angles in a cell. Convexity of
log(sum(v*s**m)) supplies lower tangents and upper chords in m. No claim
of directed-rounding eigenvalue/FEM continuum error certification.
"""
from common import *
from elastic_primitives import stress_tensor
from scipy.special import logsumexp
import heapq, time

class Response:

    def __init__(self, stress, volume):
        self.S = [stress_tensor(stress[:, :, k]) for k in range(2)]
        self.v = volume
        self.norms = [np.linalg.norm(s, axis=(1, 2)) for s in self.S]
        self.cache = {}

    def principal(self, angle):
        if angle not in self.cache:
            t = np.deg2rad(angle)
            self.cache[angle] = np.maximum(np.linalg.eigvalsh(self.S[0] * np.cos(t) + self.S[1] * np.sin(t))[:, -1], 0)
        return self.cache[angle]

    def span(self, a, b):
        c = (a + b) / 2
        t = np.deg2rad(c)
        (ta, tb) = np.deg2rad([a, b])
        dc = max(abs(np.cos(ta) - np.cos(t)), abs(np.cos(tb) - np.cos(t)))
        ds = max(abs(np.sin(ta) - np.sin(t)), abs(np.sin(tb) - np.sin(t)))
        radius = dc * self.norms[0] + ds * self.norms[1]
        p = self.principal(c)
        return (np.maximum(p - radius, 0), p + radius)

    def logH(self, p, m, derivative=False):
        ok = p > 0
        if not ok.any():
            return (-float('inf'), 0) if derivative else -float('inf')
        ls = np.log(p[ok])
        z = np.log(self.v[ok]) + m * ls
        h = logsumexp(z)
        return (h, float(np.exp(z - h) @ ls)) if derivative else float(h)

def ratio(ref, opt, angle, m):
    return (ref.logH(ref.principal(angle), m) - opt.logH(opt.principal(angle), m)) / m

def cell_bounds(ref, opt, a, b, l, h):
    (rlo, rhi) = ref.span(a, b)
    (olo, ohi) = opt.span(a, b)
    mid = (l + h) / 2

    def lower(A, ap, B, bp):
        (v, d) = A.logH(ap, mid, True)
        if not np.isfinite(v):
            return -float('inf')
        return min(((v + d * (m - mid) - B.logH(bp, m)) / m for m in [l, h]))
    lo = lower(ref, rlo, opt, ohi)
    hi = -lower(opt, olo, ref, rhi)
    return (lo, hi)

def bound_minimum(ref, opt, maxboxes=2048, tol=0.002):
    count = 0
    queue = []
    best = float('inf')
    arg = None

    def push(a, b, l, h):
        nonlocal count, best, arg
        (lo, hi) = cell_bounds(ref, opt, a, b, l, h)
        for (aa, mm) in [((a + b) / 2, (l + h) / 2), (a, l), (b, h)]:
            v = ratio(ref, opt, aa, mm)
            if v < best:
                best = v
                arg = [aa, mm]
        heapq.heappush(queue, (lo, count, a, b, l, h))
        count += 1
    push(0.0, 30.0, 2.0, 15.6)
    while count + 2 <= maxboxes and best - queue[0][0] > tol:
        (lo, _, a, b, l, h) = heapq.heappop(queue)
        if (b - a) / 30 >= (h - l) / 13.6:
            c = (a + b) / 2
            push(a, c, l, h)
            push(c, b, l, h)
        else:
            c = (l + h) / 2
            push(a, b, l, c)
            push(a, b, c, h)
    return dict(minimum_force_ratio_lower=float(np.exp(queue[0][0])), minimum_force_ratio_sample_upper=float(np.exp(best)), log_bound_gap=float(best - queue[0][0]), argmin_sample_angle_m=arg, cells_evaluated=count, active_covering_cells=len(queue), converged=bool(best - queue[0][0] <= tol), scope='Exact-arithmetic analytic parameter enclosure of stored FE fields; formal floating-point and FE continuum enclosures MISSING')

def run():
    start = time.perf_counter()
    z = np.load(DATA / 'B_FIELDS.npz')
    out = {}
    grids = {}
    angs = np.linspace(0, 30, 31)
    ms = np.linspace(2, 15.6, 35)
    for name in ['R3', 'R4', 'R5', 'optimized']:
        out[name] = {}
        grids[name] = {}
        for support in ['8600', '18000', 'rigid']:
            ref = Response(z['reference__' + support + '__stress_basis'], z['reference__' + support + '__vol'])
            opt = Response(z[name + '__' + support + '__stress_basis'], z[name + '__' + support + '__vol'])
            grid = np.array([[np.exp(ratio(ref, opt, a, m)) for m in ms] for a in angs])
            rr = bound_minimum(ref, opt) if name == 'optimized' else {}
            out[name][support] = dict(sampled_min=float(grid.min()), sampled_max=float(grid.max()), sampled_gate=bool(grid.min() > 1.02), continuous=rr)
            grids[name][support] = grid
            print(name, support, out[name][support], flush=True)
    save(ROOT / 'raw/B_BOUNDS.json', dict(claim_type='capability', contrasts=out, seconds=time.perf_counter() - start, angle_grid_deg=angs, m_grid=ms, relative_force_definition='(volume integral reference tensile^m / candidate tensile^m)^(1/m)', resolution='PER_TOOTH; underlying tensors PER_POINT', numeric_certificate='Analytical real-arithmetic only; no formal floating-point enclosure'))
    np.savez_compressed(DATA / 'B_GRIDS.npz', angles=angs, ms=ms, **{n + '__' + s: g for (n, ss) in grids.items() for (s, g) in ss.items()})
    state('B_CONTINUOUS_BOUNDS', out['optimized'], 'Preserve outcome and select a changed physical/design operation for next round')
if __name__ == '__main__':
    run()
