"""Feasible active-set contact solve in full force/compliance coordinates."""
import numpy as np

def solve(A, C, g, w, h, fstart):
    n = len(g)
    current = np.array(fstart, dtype=float).copy()
    blocked = set()
    assert np.max(np.abs(A.T @ current - w)) < 1e-07 and current.min() > -1e-08
    for iteration in range(1000):
        free = np.array([i for i in range(n) if i not in blocked], dtype=int)
        AF = A[free]
        CF = C[np.ix_(free, free)]
        mat = np.block([[CF, -AF], [AF.T, np.zeros((3, 3))]])
        sol = np.linalg.solve(mat, np.r_[h[free] - g[free], w])
        trial = np.zeros(n)
        trial[free] = sol[:len(free)]
        q = sol[len(free):]
        bad = free[trial[free] < -1e-10]
        if len(bad):
            ratios = current[bad] / (current[bad] - trial[bad])
            p = int(np.argmin(ratios))
            alpha = float(ratios[p])
            current = current + alpha * (trial - current)
            current[int(bad[p])] = 0.0
            blocked.add(int(bad[p]))
            continue
        current = trial
        slack = C @ current + g - h - A @ q
        inactive = np.array(sorted(blocked), dtype=int)
        if len(inactive) and slack[inactive].min() < -1e-12:
            blocked.remove(int(inactive[np.argmin(slack[inactive])]))
            continue
        meta = {'iterations': iteration + 1, 'min_force_N': float(current.min()), 'min_slack_mm': float(slack.min()), 'balance_residual_N': float(np.max(np.abs(A.T @ current - w))), 'complementarity_N_mm': float(np.max(np.abs(current * slack))), 'active_stationarity_mm': float(np.max(np.abs(slack[free]))), 'active_count': int(np.sum(current > 1e-08))}
        return (current, q, meta)
    raise RuntimeError('Active-set iteration cap; no fallback silently substituted')
