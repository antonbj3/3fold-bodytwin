"""Numerical rank and observable guards for a fixed contact branch.

These are local, floating-point statements at the declared relative cutoff.
They do not prove global uniqueness of a nonlinear friction problem. A right
nullspace changes the force response; a left nullspace tests compatibility of
parameter perturbations. A row Q is invariant only when Q @ N vanishes.
"""
from __future__ import annotations
import numpy as np

def rank_factor(A, rcond=1e-12):
    """Return SVD factors and numerical right/left nullspaces of square A."""
    A = np.asarray(A, dtype=np.float64)
    if A.ndim != 2 or A.shape[0] != A.shape[1] or (not np.isfinite(A).all()):
        raise ValueError('A must be a finite square matrix')
    if not 0.0 < rcond < 1.0:
        raise ValueError('rcond must be between zero and one')
    (U, s, Vt) = np.linalg.svd(A, full_matrices=True)
    cutoff = rcond * s[0] if s.size else 0.0
    rank = int(np.sum(s > cutoff))
    (N, L) = (Vt[rank:].T.copy(), U[:, rank:].copy())
    cond = float(s[0] / s[-1]) if s.size and s[-1] > 0 else float('inf') if s.size else 0.0
    info = {'status': 'REDUNDANT_CONTACTS' if rank < A.shape[1] else 'UNIQUE_LINEAR_SYSTEM', 'rank': rank, 'nullity': A.shape[1] - rank, 'nullspace': N, 'left_nullspace': L, 'singular_values': s.copy(), 'rank_cutoff': cutoff, 'rcond': rcond, 'cond': cond, 'singular': rank < A.shape[1], 'scope': 'NUMERICAL_FIXED_BRANCH_AT_SUPPLIED_IMPULSE', 'selection': None}
    return (U, s, Vt, info)

def guarded_response(A, B, factors, selection=None, observable=None, tol=1e-10):
    """Solve A X = -B without representing an ambiguous X as unique.

    selection='minimum_norm' chooses the Euclidean minimum norm *linear
    response*, not a compliant physical model or a nonlinear force selection.
    Inconsistent columns remain NaN even when a selection is requested.
    Observable rows are returned only for compatible columns and Q N = 0.
    """
    if not np.isfinite(tol) or not 0.0 < tol < 1.0:
        raise ValueError('tol must be finite and between zero and one')
    if selection not in (None, 'minimum_norm'):
        raise ValueError("selection must be None or 'minimum_norm'")
    (U, s, Vt, info) = factors
    B = np.asarray(B, dtype=np.float64)
    if B.ndim != 2 or B.shape[0] != A.shape[0] or (not np.isfinite(B).all()):
        raise ValueError('B must be a finite matrix with A.shape[0] rows')
    rank = info['rank']
    bscale = np.max(np.abs(B), axis=0, initial=0.0)
    scaled_B = B / np.where(bscale > 0, bscale, 1.0)
    scale = np.linalg.norm(scaled_B, axis=0)
    relative_mismatch = np.linalg.norm(info['left_nullspace'].T @ scaled_B, axis=0)
    compatible = relative_mismatch <= tol * scale
    with np.errstate(over='ignore'):
        mismatch = relative_mismatch * bscale
    if info['singular']:
        projected = -(Vt[:rank].T @ (U[:, :rank].T @ B / s[:rank, None]))
        X = projected.copy() if selection == 'minimum_norm' else np.full_like(projected, np.nan)
        X[:, ~compatible] = np.nan
    else:
        projected = -np.linalg.solve(A, B)
        X = projected.copy()
    finite_response = np.isfinite(projected).all(axis=0)
    X[:, ~finite_response] = np.nan
    result = {'response_finite': finite_response, 'compatible': compatible, 'compatibility_residual': mismatch, 'selection': 'MINIMUM_NORM_LINEAR_RESPONSE' if selection else None}
    if observable is not None:
        Q = np.atleast_2d(np.asarray(observable, dtype=np.float64))
        if Q.shape[1] != A.shape[1] or not np.isfinite(Q).all():
            raise ValueError('observable must be finite with A.shape[1] columns')
        qscale = np.max(np.abs(Q), axis=1, initial=0.0)
        scaled_Q = Q / np.where(qscale > 0, qscale, 1.0)[:, None]
        relative_leakage = np.linalg.norm(scaled_Q @ info['nullspace'], axis=1)
        invariant = relative_leakage <= tol * np.linalg.norm(scaled_Q, axis=1)
        with np.errstate(over='ignore', invalid='ignore'):
            leakage = relative_leakage * qscale
            Y = Q @ projected
        valid = invariant[:, None] & compatible[None, :] & np.isfinite(Y)
        Y[~valid] = np.nan
        result.update(observable=Y, observable_valid=valid, observable_invariant=invariant, observable_nullspace_residual=leakage)
    return (X, result)
