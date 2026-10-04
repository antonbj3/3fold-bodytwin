"""Frozen phenomenological log-torque model; no clinical or pressure inference."""
import numpy as np
from collections import Counter

def feature(rows, kind):
    values = []
    for r in rows:
        v = [np.log(r['rho_gcc'])]
        if kind in ['cortex', 'both']:
            v.append(r['cortex_mm'] / r['L_mm'] * (r['cortex_rho_gcc'] / r['rho_gcc']))
        if kind in ['drilling', 'both']:
            v.append(1 - r['Df_mm'] / r['D_mm'])
        values.append(v)
    return np.asarray(values)

def fit(train, kind, penalty=1.0):
    raw = feature(train, kind)
    counts = Counter((r['study'] for r in train))
    w = np.asarray([1 / counts[r['study']] for r in train])
    w = w / w.sum() * len(set((r['study'] for r in train)))
    mean = np.average(raw, axis=0, weights=w)
    std = np.sqrt(np.average((raw - mean) ** 2, axis=0, weights=w))
    std = np.where(std > 1e-12, std, 1.0)
    X = np.column_stack([np.ones(len(train)), (raw - mean) / std])
    y = np.log([r['torque_Ncm'] / (r['D_mm'] * r['L_mm']) for r in train])
    A = X.T @ (w[:, None] * X)
    b = X.T @ (w * y)
    reg = np.eye(X.shape[1]) * penalty
    reg[0, 0] = 0
    beta = np.linalg.lstsq(A + reg, b, rcond=None)[0]
    return dict(kind=kind, beta=beta, mean=mean, std=std, rank=int(np.linalg.matrix_rank(X)), columns=X.shape[1], condition=float(np.linalg.cond(A + reg)))

def predict(model, rows):
    X = np.column_stack([np.ones(len(rows)), (feature(rows, model['kind']) - model['mean']) / model['std']])
    return np.exp(X @ model['beta']) * np.asarray([r['D_mm'] * r['L_mm'] for r in rows])

def json_model(m):
    return {k: v.tolist() if hasattr(v, 'tolist') else v for (k, v) in m.items()}
