import numpy as np

def closest(P, T):
    a = T[:, 0]
    e = T[:, 1] - a
    f = T[:, 2] - a
    w = P - a
    ee = np.einsum('ij,ij->i', e, e)
    ff = np.einsum('ij,ij->i', f, f)
    ef = np.einsum('ij,ij->i', e, f)
    we = np.einsum('ij,ij->i', w, e)
    wf = np.einsum('ij,ij->i', w, f)
    den = ee * ff - ef * ef
    ok = den > 1e-24
    u = np.divide(we * ff - wf * ef, den, out=np.zeros(len(P)), where=ok)
    v = np.divide(wf * ee - we * ef, den, out=np.zeros(len(P)), where=ok)
    Q = a + u[:, None] * e + v[:, None] * f
    d = np.linalg.norm(P - Q, axis=1)
    valid = ok & (u >= 0) & (v >= 0) & (u + v <= 1)
    best = np.where(valid, d, np.inf)
    bestQ = Q.copy()
    for (i, j) in [(0, 1), (1, 2), (2, 0)]:
        aa = T[:, i]
        bb = T[:, j]
        ed = bb - aa
        ed2 = np.einsum('ij,ij->i', ed, ed)
        t = np.clip(np.divide(np.einsum('ij,ij->i', P - aa, ed), ed2, out=np.zeros(len(P)), where=ed2 > 0), 0, 1)
        q = aa + t[:, None] * ed
        dd = np.linalg.norm(P - q, axis=1)
        take = dd < best
        best[take] = dd[take]
        bestQ[take] = q[take]
    return (best, bestQ)
