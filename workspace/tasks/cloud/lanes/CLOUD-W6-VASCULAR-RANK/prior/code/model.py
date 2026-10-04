"""Original two-compartment (plasma/interstitium) synthetic model.

States: Vp, Vi [L]; Mp, Mi [g] (solute mass). Closed system:
  J  = kf * ((Vp-Vp_eq)/Vp_eq - (Vi-Vi_eq)/Vi_eq)      fluid filtration p->i [L/min]
  F  = PS*(Cp-Ci) + (1-sigma)*J*C_up                   solute flux p->i [g/min]
       (C_up = Cp if J>=0 else Ci; upwind convection)
  dVp=-J, dVi=+J, dMp=-F, dMi=+F   => Vp+Vi and Mp+Mi are conserved exactly.
Scenario: 1.0 L solute-free bolus into plasma at t=0 (Vp0 = Vp_eq+1), known
pre-bolus plasma concentration; Vi starts at Vi_eq; Ci0 unknown.
Unknown theta (log-parametrised): kf [L/min], PS [L/min], Vi_eq [L], Ci0 [g/L].
All numbers are synthetic, chosen to be physiologically plausible orders of
magnitude only; none are fitted to real data.
"""
import numpy as np
from scipy.integrate import solve_ivp

VP_EQ = 3.0        # L, known
BOLUS = 1.0        # L solute-free, known
CP_PRE = 40.0      # g/L, known pre-bolus plasma concentration
SIGMA = 0.9        # reflection coefficient, known
NAMES = ["kf", "PS", "Vi_eq", "Ci0"]
THETA_NOM = np.array([0.03, 0.003, 12.0, 20.0])
SD_P = 0.4         # g/L known additive noise, plasma samples
SD_I = 0.4         # g/L known additive noise, interstitial samples
CAND = np.arange(5.0, 240.0 + 1e-9, 5.0)   # candidate sample times [min]


def rhs(t, y, kf, PS, Vi_eq):
    Vp, Vi, Mp, Mi = y
    Cp, Ci = Mp / Vp, Mi / Vi
    J = kf * ((Vp - VP_EQ) / VP_EQ - (Vi - Vi_eq) / Vi_eq)
    F = PS * (Cp - Ci) + (1 - SIGMA) * J * (Cp if J >= 0 else Ci)
    return [-J, J, -F, F]


def simulate(theta, t):
    kf, PS, Vi_eq, Ci0 = theta
    Vp0 = VP_EQ + BOLUS
    y0 = [Vp0, Vi_eq, CP_PRE * VP_EQ, Ci0 * Vi_eq]
    sol = solve_ivp(rhs, (0, t.max()), y0, t_eval=t, args=(kf, PS, Vi_eq),
                    method="LSODA", rtol=1e-11, atol=1e-12)
    Vp, Vi, Mp, Mi = sol.y
    return Mp / Vp, Mi / Vi, sol.y


def observe(logtheta, tp, ti):
    """Concatenated plasma Cp(tp) and interstitial Ci(ti) predictions."""
    theta = np.exp(logtheta)
    t = np.unique(np.concatenate([tp, ti, [1.0]]))
    Cp, Ci, _ = simulate(theta, t)
    idx = {v: k for k, v in enumerate(t)}
    return np.concatenate([Cp[[idx[x] for x in tp]], Ci[[idx[x] for x in ti]]])


def jac(logtheta, tp, ti, h=1e-4):
    """Central finite-difference sensitivities wrt log-parameters."""
    cols = []
    for k in range(len(logtheta)):
        e = np.zeros_like(logtheta); e[k] = h
        cols.append((observe(logtheta + e, tp, ti) - observe(logtheta - e, tp, ti)) / (2 * h))
    return np.array(cols).T


def fisher(logtheta, tp, ti):
    J = jac(logtheta, tp, ti)
    w = np.concatenate([np.full(len(tp), 1 / SD_P**2), np.full(len(ti), 1 / SD_I**2)])
    return J.T @ (J * w[:, None])


RANK_TOL = 1e-6  # singular value / max singular value (preregistered in code before runs)


def fim_rank(F):
    s = np.linalg.svd(F, compute_uv=False)
    return int(np.sum(s > RANK_TOL * s[0])), s
