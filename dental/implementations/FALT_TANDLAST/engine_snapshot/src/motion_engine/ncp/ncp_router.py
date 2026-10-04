"""Routing an exact-cone contact problem to a solver, and the spectral rho rules.

Three pieces, in the order the router uses them:

  1. `bilateral_certificate` -- solve the active closed set as if it were bilateral,
     dropping contacts with lam_n <= 0 until the set stops changing, then ask whether
     every remaining contact is strictly inside its friction cone and no dropped contact
     is penetrating. If it is interior, that solve IS the answer and no iteration runs.
  2. `spectral_rho` -- lambda_min^+ and lambda_max of G restricted to the active rows,
     and the two candidate step sizes built from them:
         rule A   rho = lambda_min^+
         rule B   rho = sqrt(lambda_min^+ lambda_max)
     with the structure choice "A when the contact graph is a chain (max body degree <= 2
     and the pencil null space is trivial), otherwise B".
  3. `route` -- the decision tree over those two.

HOW WELL THE SPECTRAL RULE PREDICTS rho*, MEASURED AGAINST A REAL SWEEP.
A 21-scene by 25-rho sweep (`data/ncp/rho_oracle_sweep.csv`, bit-identical over two runs,
525/525) puts the prediction within 2x of the measured rho* as follows:

    rule A            2/21
    rule B            7/21
    structure choice  5/21
    hindsight best    8/21   (picking the better of A and B after the fact; not a rule)

So the spectral rule does not hold. The earlier "18/21" was the prediction pasted in as
the oracle on 11 of the rows; with a real sweep it is 5/21. What does stand is the
spectrum identity itself -- spec(G) equals the spectrum of the (patch) line-graph pencil
to 1e-15 -- and the lumping that makes it exact; only the step-size rule built on top
of it fails. 223 of 525 grid points are censored at the 5000-iteration cap and 6 scenes
have their minimum on the edge of the grid, so a censored scene counts as a miss, not a
hit, and an edge minimum identifies no optimum outside the grid.

THE ROUTER'S OWN AUDIT RESERVATION.
Routing 21 benchmark scenes takes the mis-routed count from 4 to 1 with the articulated
branch and to 0 with the generalised "n_c <= 4 goes to PGS" rule, leave-one-out 0/21 --
but "articulated" was detected on the scene NAME. With the name test removed, a real
articulated quadruped is mis-routed by 28-39x and a six-contact articulated scene by 52x:
the CG sweep the lane used had no discriminating power (CG 6 is bit-identical to CG 256
on one scene, CG 4 equals CG 16 on a packing at a third of the time). `route` therefore
takes `articulated` as an explicit argument; it does not guess from a key.

Measured against the sweep oracle rather than a hard-coded table, the router's median
cost ratio to the per-scene best is 1.005x with a worst case of 9.2x, and it converges
where a fixed 100-sweep PGS converges on only 8 of 21 scenes.
"""
from __future__ import annotations
import numpy as np
__all__ = ['bilateral_certificate', 'spectral_rho', 'structure_class', 'route', 'RULE_HIT_RATES']
RULE_HIT_RATES = {'A': (2, 21), 'B': (7, 21), 'structure': (5, 21), 'hindsight_best': (8, 21)}

def bilateral_certificate(Geta: np.ndarray, b: np.ndarray, mu: np.ndarray, max_rounds: int=12, pen_tol: float=1e-07) -> dict:
    """Solve the closed set bilaterally, dropping lam_n <= 0, and certify interiority.

    Redundant closed forces use a stated Euclidean minimum-norm selection;
    the certificate proves feasibility of that representative, not uniqueness.
    An inconsistent least-squares system is never certified.

    Returns `is_interior`: every closed contact strictly inside its cone AND no open
    contact penetrating. When true, `lam` is the solution and no cone iteration is needed.
    """
    import scipy.linalg as sla
    from .contact_rank import rank_factor, guarded_response
    n_c = len(mu)
    Geta = np.asarray(Geta, float)
    b = np.asarray(b, float)
    mu = np.asarray(mu, float)
    closed = np.ones(n_c, bool)
    (rounds, L, x, idx) = (0, np.zeros((n_c, 3)), None, None)
    rank_info = rank_factor(np.zeros((0, 0)))[3]
    (compatible, stable) = (True, False)
    selection = 'UNIQUE_BILATERAL'
    for rounds in range(1, max_rounds + 1):
        idx = np.repeat(np.where(closed)[0] * 3, 3) + np.tile([0, 1, 2], closed.sum())
        if idx.size == 0:
            x = np.zeros(0)
            L = np.zeros((n_c, 3))
            rank_info = rank_factor(np.zeros((0, 0)))[3]
            (compatible, stable) = (True, True)
            selection = 'UNIQUE_BILATERAL'
            break
        A = Geta[np.ix_(idx, idx)]
        factors = rank_factor(A)
        rank_info = factors[3]
        if rank_info['singular']:
            (response, response_info) = guarded_response(A, b[idx, None], factors, selection='minimum_norm')
            compatible = bool(response_info['compatible'].all())
            if not compatible:
                x = None
                break
            x = response[:, 0]
            selection = 'MINIMUM_NORM_BILATERAL'
        else:
            compatible = True
            selection = 'UNIQUE_BILATERAL'
            try:
                c = sla.cho_factor(np.array(A, float, order='F'), lower=True, check_finite=False)
                x = sla.cho_solve(c, -b[idx], check_finite=False)
            except np.linalg.LinAlgError:
                x = np.linalg.solve(A, -b[idx])
        L = np.zeros((n_c, 3))
        L[closed] = x.reshape(-1, 3)
        neg = closed & (L[:, 0] <= 0)
        if not neg.any():
            stable = True
            break
        closed = closed & ~neg
    (ln, lt) = (L[:, 0], np.hypot(L[:, 1], L[:, 2]))
    out_cone = closed & (lt >= mu * ln)
    lam_full = np.zeros(3 * n_c)
    if idx is not None and x is not None:
        lam_full[idx] = x
    u_n_open = (Geta @ lam_full + b)[0::3][~closed]
    penetrating_open = bool(np.any(u_n_open < -pen_tol)) if u_n_open.size else False
    N = np.zeros((3 * n_c, rank_info['nullity']))
    if idx is not None:
        N[idx] = rank_info['nullspace']
    closed_residual = float(np.max(np.abs((Geta @ lam_full + b)[np.repeat(closed, 3)]))) if closed.any() else 0.0
    equation_ok = compatible and closed_residual <= pen_tol
    return {'status': rank_info['status'] if compatible else 'INCOMPATIBLE_BILATERAL', 'nullspace': N, 'force_selection': selection if compatible else None, 'rank_cutoff': rank_info['rank_cutoff'], 'compatible': compatible, 'closed_equation_residual': closed_residual, 'n_closed': int(closed.sum()), 'n_out': int(out_cone.sum()), 'is_interior': bool(stable and equation_ok and (out_cone.sum() == 0) and (not penetrating_open)), 'rounds': int(rounds), 'lam': lam_full, 'closed': closed}

def spectral_rho(Geta: np.ndarray, pos_tol_rel: float=1e-11) -> dict:
    """lambda_min^+, lambda_max and the two candidate rho.

    lambda_min^+ is the smallest eigenvalue above `pos_tol_rel * lambda_max`, i.e. the
    smallest one that is not in the numerical null space of the pencil.
    """
    ev = np.linalg.eigvalsh(np.asarray(Geta, float))
    lmax = float(ev[-1])
    pos = ev[ev > pos_tol_rel * lmax]
    lmin = float(pos[0]) if pos.size else lmax
    return {'lambda_min_positive': lmin, 'lambda_max': lmax, 'rho_A': lmin, 'rho_B': float(np.sqrt(lmin * lmax))}

def structure_class(J: np.ndarray, blk: int=3, dof: int=6) -> dict:
    """The structure test the rho rule branches on, read off J.

    C1 is true when all three hold:
      * no body carries more than two contacts (`degree_max <= 2`),
      * no two contacts share more than one body (`max_shared <= 1`),
      * on every body with exactly two contacts the two Jacobian blocks are exact
        negatives of each other (`opposite`).
    Degree alone is not enough: a two-contact scene whose blocks are not opposite is
    not a chain in the sense the Laplacian identity needs, and the sweep's `humanoid`
    row is exactly that case -- degree 2 and C1 = 0.
    """
    J = np.asarray(J, float)
    (n_c, n_b) = (J.shape[0] // blk, J.shape[1] // dof)
    supports = [set(np.flatnonzero(np.abs(J[r * blk:(r + 1) * blk]).sum(axis=0)) // dof) for r in range(n_c)]
    contacts = [[c for (c, s) in enumerate(supports) if b in s] for b in range(n_b)]
    degree = [len(cs) for cs in contacts]
    shared = max((len(supports[c] & supports[d]) for c in range(n_c) for d in range(c)), default=0)
    opposite = True
    for (b, cs) in enumerate(contacts):
        if len(cs) == 2:
            (c, d) = cs
            x = J[c * blk:(c + 1) * blk, b * dof:(b + 1) * dof]
            y = J[d * blk:(d + 1) * blk, b * dof:(b + 1) * dof]
            if np.any(np.abs(x + y) > 1e-14 + 1e-12 * np.abs(y)):
                opposite = False
    deg = max(degree, default=0)
    return {'n_c': n_c, 'degree_max': int(deg), 'max_shared': int(shared), 'opposite': bool(opposite), 'C1': int(deg <= 2 and shared <= 1 and opposite)}

def route(Geta: np.ndarray, b: np.ndarray, mu: np.ndarray, J: np.ndarray, dof: int=6, eta_positive: bool=False, chain: bool=False, articulated: bool=False) -> dict:
    """Pick a branch. `articulated` is an argument, never inferred from a name."""
    n_c = len(mu)
    cert = bilateral_certificate(Geta, b, mu)
    if cert['is_interior']:
        return {'branch': 'exact_bilateral', 'rho': None, 'certificate': cert}
    spec = spectral_rho(Geta)
    struct = structure_class(J, 3, dof)
    if eta_positive:
        return {'branch': 'ruiz_admm', 'rho': 0.1, 'certificate': cert, 'spectrum': spec, 'structure': struct}
    if chain:
        return {'branch': 'pgs_arc_freeze', 'rho': None, 'certificate': cert, 'spectrum': spec, 'structure': struct}
    if articulated and n_c <= 4:
        return {'branch': 'pgs_certificate_stop', 'rho': None, 'certificate': cert, 'spectrum': spec, 'structure': struct}
    rho = spec['rho_A'] if struct['C1'] else spec['rho_B']
    branch = 'spectral_admm' if cert['n_out'] == 0 else 'ladder_admm'
    return {'branch': branch, 'rho': rho, 'certificate': cert, 'spectrum': spec, 'structure': struct, 'rho_rule': 'A' if struct['C1'] else 'B'}
