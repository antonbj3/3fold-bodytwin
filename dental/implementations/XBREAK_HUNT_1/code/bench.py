"""Declared spring/platen instrument. No patient support is measured here."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
import hashlib, json, numpy as np
from scipy.optimize import minimize
ROOT = Path(__file__).resolve().parents[1]
X21 = Path(_release_expand('@DENTAL_IMPLEMENTATIONS@/X21'))

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def write(path, obj):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False))

def check_frozen(tag):
    for line in (ROOT / 'raw' / f'FROZEN_{tag}.sha256').read_text().splitlines():
        (digest, p) = line.split(maxsplit=1)
        assert sha(ROOT / p) == digest, ('Frozen artifact drift', p)

def geometry(case):
    path = X21 / 'raw/cases' / f'{case}.json'
    s = json.loads(path.read_text())
    selected = {}
    for pair in s['pairs']:
        t = pair['upper_fdi']
        if t not in selected or pair['minimum_projected_gap_mm'] < selected[t]['minimum_projected_gap_mm']:
            selected[t] = pair
    ids = sorted(selected)
    xy = np.array([selected[t]['witness']['xy_frame_mm'] for t in ids])
    A = np.column_stack((np.ones(len(ids)), xy / 30.0))
    M = np.array([A[:, i] * A[:, j] for (i, j) in [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]])
    return (ids, A, M, {'path': str(path), 'sha256': sha(path), 'xy_mm': xy.tolist(), 'predicted_fdi': ids, 'paired_lower_fdi': [selected[t]['lower_fdi'] for t in ids], 'geometry_scope': 'Real projected witness positions only; original gaps replaced by a leveled bench g=0', 'original_gaps_mm': [selected[t]['minimum_projected_gap_mm'] for t in ids]})

def fixed(A, k, w, h):
    K = A.T @ (k[:, None] * A)
    q = np.linalg.solve(K, w - A.T @ (k * h))
    c = A @ q + h
    if c.min() <= 0:
        raise ValueError('Fixed-active formula outside full contact')
    return (q, k * c, K)

def unilateral(A, k, w, h, q_start):

    def fun(q):
        c = np.maximum(A @ q + h, 0)
        return 0.5 * np.dot(k * c, c) - w @ q

    def jac(q):
        c = np.maximum(A @ q + h, 0)
        return A.T @ (k * c) - w

    def hess(q):
        active = A @ q + h > 0
        return A.T @ ((k * active)[:, None] * A)
    sol = minimize(fun, q_start, jac=jac, hess=hess, method='trust-exact', options={'gtol': 1e-10, 'maxiter': 200})
    f = k * np.maximum(A @ sol.x + h, 0)
    res = np.max(np.abs(A.T @ f - w))
    return (sol.x, f, {'success': bool(sol.success), 'message': str(sol.message), 'equilibrium_residual_N': float(res), 'iterations': int(sol.nit)})

def base_w(A, k):
    d = 300.0 / k.sum()
    return (A.T @ (k * np.full(len(k), d)), np.array([d, 0, 0]))
