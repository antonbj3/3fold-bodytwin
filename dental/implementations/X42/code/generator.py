"""Participant-only constrained roof generator; all geometry is public.

Reuses the published V3 task representation and legacy LP feasible set.
No references, source archives or scorer are imported here.
"""
import numpy as np
from scipy.sparse import csr_matrix, eye, hstack, vstack
from scipy.optimize import linprog
from legacy.generators import submit
from legacy.geometry import section_weights

def structural_refusal(t):
    r = t['requirements']
    if t['family'] == 'implant_crown':
        a = np.array(t['implant_axis'])
        ang = np.degrees(np.arccos(a[2] / np.linalg.norm(a)))
        if ang > r['channel_max_deg'] + 1e-06:
            return 'fixed implant channel'
    if t['family'] == 'lattice_onlay' and t['minimum_strut_mm'] < r['strut_mm'] - 1e-06:
        return 'fixed strut width'
    return None

def project(t, target):
    if structural_refusal(t):
        return {'status': 'ABSTAIN', 'reason': structural_refusal(t), 'task_id': t['task_id']}
    n = len(target)
    r = t['requirements']
    inner = t['preparation_z'] + r['film_min_mm']
    lower = inner + r['wall_mm']
    z = np.maximum(target, lower)
    A = t['A']
    b = t['obstacle_b']
    rows = []
    rhs = []
    if A.shape[0] and np.min(b - A @ lower) < -1e-06:
        return {'status': 'ABSTAIN', 'reason': 'wall-film obstacle obstruction', 'task_id': t['task_id']}
    feasible = not A.shape[0] or np.max(A @ z - b) <= 0
    if t['family'] == 'bridge3':
        for x in t['connector_x']:
            w = section_weights(t['xy'], t['faces'], x)
            rows.append(csr_matrix(np.r_[-w, np.zeros(n)][None, :]))
            rhs.append(np.array([-r['connector_mm2'] - w @ inner]))
            feasible &= w @ (z - inner) >= r['connector_mm2']
    if feasible:
        return submit(t, z, inner, {'projection': 'separable lower-bound certificate'})
    I = eye(n, format='csr')
    rows = [hstack([A, csr_matrix((A.shape[0], n))]), hstack([I, -I]), hstack([-I, -I])] + rows
    rhs = [b, target, -target] + rhs
    sol = linprog(np.r_[np.zeros(n), t['weights']], A_ub=vstack(rows), b_ub=np.concatenate(rhs), bounds=[(v, None) for v in lower] + [(0, None)] * n, method='highs', options={'threads': 1})
    if not sol.success:
        return {'status': 'ABSTAIN', 'reason': 'LP ' + str(sol.status), 'task_id': t['task_id']}
    return submit(t, sol.x[:n], inner, {'projection': 'exact constrained target LP', 'iterations': int(sol.nit)})

def features(t):
    p = np.array(t['prior'])
    c = np.array(t['ceiling'])
    ok = np.isfinite(c)
    gap = np.where(ok, c - p, 0.0)
    half = np.max(np.abs(t['xy']), axis=0)
    return np.r_[half, np.median(p), p - np.median(p), np.clip(gap, -10, 10), ok.astype(float)]

def predict(t, model):
    x = (features(t) - model['xmean']) / model['xscale']
    kind = str(model['kind'])
    if kind == 'kernel':
        xf = model['xtrain']
        d = np.mean((xf - x) ** 2, axis=1)
        k = np.exp(-float(model['gamma']) * d)
        res = model['ymean'] + k @ model['coef'] @ model['basis']
    elif kind == 'ridge':
        res = model['ymean'] + x @ model['coef']
    elif kind == 'mean':
        res = model['ymean'].copy()
    elif kind == 'forest':
        x = x.astype(np.float32)
        pred = []
        for i in range(int(model['ntrees'])):
            prefix = 'tree' + str(i) + '_'
            left = model[prefix + 'left']
            right = model[prefix + 'right']
            feature = model[prefix + 'feature']
            threshold = model[prefix + 'threshold']
            values = model[prefix + 'values']
            j = 0
            while left[j] != -1:
                j = int(left[j] if x[int(feature[j])] <= threshold[j] else right[j])
            pred.append(values[j])
        res = np.mean(pred, axis=0)
    else:
        raise ValueError(kind)
    return np.array(t['prior']) + res

def conditional(t, model):
    return project(t, predict(t, model))

def branch_target(t, shape_model, contact_model, beta):
    z = predict(t, shape_model)
    if beta == 0:
        return z
    p = np.clip(predict(t, contact_model) - t['prior'], 0, 1)
    w = t['weights']
    ceil = t['ceiling']
    ok = np.isfinite(ceil)
    targetarea = float(np.sum(w * p * (p >= 0.5) * ok))
    order = np.flatnonzero(ok)
    order = order[np.argsort(-p[order], kind='stable')]
    accum = np.r_[0.0, np.cumsum(w[order])]
    k = int(np.argmin(np.abs(accum - targetarea)))
    mask = np.zeros(len(z), bool)
    mask[order[:k]] = True
    shaped = z.copy()
    shaped[ok & mask] = ceil[ok & mask] - 0.06
    shaped[ok & ~mask] = np.minimum(z[ok & ~mask], ceil[ok & ~mask] - 0.10001)
    return (1 - beta) * z + beta * shaped

def branch(t, shape_model, contact_model, beta):
    return project(t, branch_target(t, shape_model, contact_model, beta))
