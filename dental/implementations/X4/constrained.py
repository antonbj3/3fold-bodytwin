"""Move finite graft thickness inside the straight segment search."""
import itertools
import numpy as np

def cut_feasible(nodes, radius=6.0):
    u = np.diff(nodes, axis=0)
    lengths = np.linalg.norm(u, axis=1)
    u = u / lengths[:, None]
    if len(u) == 1:
        return (bool(lengths[0] >= 10), [0.0, 0.0])
    cos = np.sum(u[:-1] * u[1:], axis=1)
    extent = radius * np.sqrt(np.maximum(1 - cos, 0) / np.maximum(1 + cos, 1e-15))
    trim = np.r_[0.0, extent, 0.0]
    good = np.all(lengths >= 10) and np.all(lengths > trim[:-1] + trim[1:] + 1e-06)
    return (bool(good), trim.tolist())

def constrained_plan(curve, E, radius=6.0, force_k=None):
    n = len(curve)
    best = {}
    details = {}
    examined = 0
    rejected = 0
    for k in [force_k] if force_k else range(1, 4):
        value = np.inf
        bestidx = None
        besttrim = None
        for inside in itertools.combinations(range(1, n - 1), k - 1):
            ids = (0,) + inside + (n - 1,)
            err = max((E[ids[i], ids[i + 1]] for i in range(k)))
            examined += 1
            if not np.isfinite(err):
                continue
            (ok, trim) = cut_feasible(curve[list(ids)], radius)
            if not ok:
                rejected += 1
                continue
            if err < value:
                value = err
                bestidx = ids
                besttrim = trim
        if bestidx:
            best[k] = bestidx
            details[k] = dict(max_chord_error_mm=float(value), miter_axis_extents_mm=besttrim)
    if not best:
        raise ValueError('UNKNOWN_NO_THICKNESS_FEASIBLE_PARTITION')
    reaches = [k for k in best if details[k]['max_chord_error_mm'] <= 2.0]
    k = min(reaches) if reaches else min(best, key=lambda k: details[k]['max_chord_error_mm'])
    ids = list(best[k])
    nodes = curve[ids]
    return (nodes, dict(count=k, indices=ids, **details[k], partitions_examined=examined, thickness_rejections=rejected, minimum_count_reaches_2mm=bool(reaches), errors_by_count={str(i): details.get(i, {}).get('max_chord_error_mm') for i in range(1, 4)}, length_mm=float(np.linalg.norm(np.diff(nodes, axis=0), axis=1).sum()), operation='finite cylinder miter extents are a constraint inside candidate search; not posthoc rejection'))
