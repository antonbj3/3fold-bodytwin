from common import *
from geometry import gap
import heapq, time

def support(t):
    ff = t['grid_faces']
    p = t['original_gap']
    r = t['reference_gap']
    return ff[np.isfinite(p[ff]).all(1) & np.isfinite(r[ff]).all(1)]

def score(q, t, p, r=None):
    if r is None:
        r = t['reference_gap']
    ff = support(t)
    (P, pa, _, _) = q.polygons(t['xy'], ff, p)
    (_, ra, _, _) = q.polygons(t['xy'], ff, r)
    overlap = q.intersection_area(t['xy'], ff, P, r)
    return max(0.0, float(pa.sum() + ra.sum() - 2 * overlap))

def set_bounds(q, t, p, r, dp, dr):
    """Real-arithmetic enclosure for arbitrary independent point errors; common offset is included."""
    ff = support(t)
    xy = t['xy']
    band = 0.1
    (Dp, adp, _, _) = q.polygons(xy, ff, p, dp, band - dp) if dp <= band / 2 else ({}, np.array([]), [], [])
    (Dr, adr, _, _) = q.polygons(xy, ff, r, dr, band - dr) if dr <= band / 2 else ({}, np.array([]), [], [])
    (Up, aup, _, _) = q.polygons(xy, ff, p, -dp, band + dp)
    (Ur, aur, _, _) = q.polygons(xy, ff, r, -dr, band + dr)
    low = adp.sum() + adr.sum() - q.intersection_area(xy, ff, Dp, r, -dr, band + dr) - q.intersection_area(xy, ff, Dr, p, -dp, band + dp)
    upper = aup.sum() + aur.sum() - q.intersection_area(xy, ff, Up, r, -dr, band + dr) - q.intersection_area(xy, ff, Dp, r, dr, band - dr)
    return (max(0.0, float(low)), max(0.0, float(upper)))

def uniform(q, t, r, pr):
    start = time.perf_counter()
    cap = r['relief_cap_mm']
    cache = {}

    def at(c):
        key = float(c)
        if key not in cache:
            cache[key] = score(q, t, gap(t, t['taper'] * key))
        return cache[key]
    pts = np.linspace(0, cap, 17)
    best = min(pts, key=lambda c: (at(c), c))
    value = at(best)
    queue = []
    counter = 0

    def push(a, b):
        nonlocal counter, best, value
        mid = (a + b) / 2
        v = at(mid)
        if v < value:
            best = mid
            value = v
        (low, _) = set_bounds(q, t, gap(t, t['taper'] * mid), t['reference_gap'], (b - a) / 2, 0.0)
        heapq.heappush(queue, (low, counter, a, b))
        counter += 1
    push(0.0, cap)
    while queue and value - queue[0][0] > pr['uniform_area_bracket_tolerance_mm2'] and (counter < pr['uniform_max_intervals']):
        (low, _, a, b) = heapq.heappop(queue)
        if low >= value:
            continue
        mid = (a + b) / 2
        push(a, mid)
        push(mid, b)
    lower = min(value, queue[0][0]) if queue else value
    return (t['taper'] * best, dict(offset_mm=float(best), symdiff_optimum_bracket_mm2=[lower, value], bracket_width_mm2=value - lower, complete_parameter_cover=True, status='CONDITIONAL_REAL_ARITHMETIC_BRACKET' if value - lower <= pr['uniform_area_bracket_tolerance_mm2'] else 'UNKNOWN_TOLERANCE_NOT_REACHED', evaluations=len(cache), intervals=counter, seconds=time.perf_counter() - start, rigorous_machine_enclosure='MISSING', max_cap_mm=cap))

def box(q, t, p, r, delta, pr):
    start = time.perf_counter()
    queue = []
    count = 0
    best = -1.0
    argmax = None
    low_best = 0.0
    cache = {}

    def at(s):
        nonlocal best, argmax
        if s not in cache:
            cache[s] = score(q, t, p + s, r + s)
        if cache[s] > best:
            best = cache[s]
            argmax = s
        return cache[s]

    def push(a, b):
        nonlocal count, low_best
        mid = (a + b) / 2
        at(a)
        at(mid)
        at(b)
        (low, upper) = set_bounds(q, t, p + mid, r + mid, (b - a) / 2, (b - a) / 2)
        low_best = max(low_best, low)
        heapq.heappush(queue, (-upper, count, a, b))
        count += 1
    push(-delta, delta)
    while queue and -queue[0][0] - best > pr['box_area_bracket_tolerance_mm2'] and (count < pr['box_max_intervals']):
        (neg, _, a, b) = heapq.heappop(queue)
        if -neg <= best:
            continue
        mid = (a + b) / 2
        push(a, mid)
        push(mid, b)
    upper = max(best, -queue[0][0]) if queue else best
    (all_low, all_up) = set_bounds(q, t, p, r, delta, delta)
    return dict(delta_mm=delta, worst_symdiff_bracket_mm2=[best, upper], all_offsets_symdiff_lower_mm2=all_low, observed_argmax_offset_mm=argmax, bracket_width_mm2=upper - best, complete_interval_cover=True, status='CONDITIONAL_REAL_ARITHMETIC_ENCLOSURE' if upper - best <= pr['box_area_bracket_tolerance_mm2'] else 'UNKNOWN_TOLERANCE_NOT_REACHED', rigorous_machine_enclosure='MISSING', source_geometry_enclosure='MISSING', intervals=count, evaluations=len(cache), seconds=time.perf_counter() - start)
