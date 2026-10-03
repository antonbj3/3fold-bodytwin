"""Minimal coupled ATP/ADP reaction-diffusion toy cell (original code, stdlib only).

Geometry (spherically symmetric, nondimensional, cell radius R = 1):
    r in [r_nuc, 1]     cytoplasm; r = r_nuc is an impermeable nuclear envelope,
                        r = 1 an impermeable plasma membrane (no adenylate flux).
    [m_in, m_out]       mitochondrial shell: ATP synthesis  P -> T
    [c_in, 1]           sub-membrane shell: ATP consumption T -> P (e.g. pumps)

Nondimensional equations (length / R, concentration / A0 = initial total
adenylate, time / (R^2 / D_T)):
    dT/dt = div(grad T)      + Da_p*chi_m*P/(kappa+P) - Da_c*chi_c*T
    dP/dt = delta*div(grad P) - Da_p*chi_m*P/(kappa+P) + Da_c*chi_c*T
Groups: Da_c = k_c R^2/D_T, Da_p = V_p R^2/(D_T A0), kappa = K_m/A0,
        delta = D_P/D_T.  All values in this file are ASSUMED, not calibrated.

Discretisation: cell-centred finite volumes on spherical shells (exact shell
volumes, face areas r_f^2; the 4*pi factor is dropped everywhere), region
indicators chi are exact volume fractions of overlap, backward Euler in time,
Picard (lagged-denominator) linearisation of the Michaelis-Menten term
iterated to 1e-12 (Newton optional; see counterexample), 2x2 block-tridiagonal solve.
"""
import math

BASE = dict(Da_c=5.0, Da_p=50.0, kappa=0.1, delta=1.2,
            r_nuc=0.3, m_in=0.35, m_out=0.55, c_in=0.85, T0=0.9)


def shell_vol(a, b):
    return (b ** 3 - a ** 3) / 3.0


def overlap_frac(a, b, lo, hi):
    x, y = max(a, lo), min(b, hi)
    return shell_vol(x, y) / shell_vol(a, b) if y > x else 0.0


class Mesh:
    def __init__(self, N, p, wrong_area=False):
        self.N = N
        h = (1.0 - p["r_nuc"]) / N
        self.rf = [p["r_nuc"] + i * h for i in range(N + 1)]
        self.rc = [0.5 * (self.rf[i] + self.rf[i + 1]) for i in range(N)]
        self.V = [shell_vol(self.rf[i], self.rf[i + 1]) for i in range(N)]
        # conductance per unit diffusivity at interior faces
        if wrong_area:  # deliberate non-conservative variant (counterexample)
            self.Gl = [0.0] + [self.rc[i] ** 2 / h for i in range(1, N)]
            self.Gr = [self.rc[i] ** 2 / h for i in range(N - 1)] + [0.0]
        else:
            g = [self.rf[i] ** 2 / h for i in range(1, N)]
            self.Gl = [0.0] + g
            self.Gr = g + [0.0]
        self.chim = [overlap_frac(self.rf[i], self.rf[i + 1], p["m_in"], p["m_out"]) for i in range(N)]
        self.chic = [overlap_frac(self.rf[i], self.rf[i + 1], p["c_in"], 1.0) for i in range(N)]


def solve_block(A, B, C, d):
    """Block Thomas: A[i] sub, B[i] diag, C[i] super (2x2 diagonal-or-full), d[i] 2-vectors."""
    n = len(B)
    def inv(m):
        det = m[0][0] * m[1][1] - m[0][1] * m[1][0]
        return [[m[1][1] / det, -m[0][1] / det], [-m[1][0] / det, m[0][0] / det]]
    def mm(a, b):
        return [[a[0][0]*b[0][0]+a[0][1]*b[1][0], a[0][0]*b[0][1]+a[0][1]*b[1][1]],
                [a[1][0]*b[0][0]+a[1][1]*b[1][0], a[1][0]*b[0][1]+a[1][1]*b[1][1]]]
    def mv(a, v):
        return [a[0][0]*v[0]+a[0][1]*v[1], a[1][0]*v[0]+a[1][1]*v[1]]
    Cp, dp = [None] * n, [None] * n
    Bi = inv(B[0]); Cp[0] = mm(Bi, C[0]); dp[0] = mv(Bi, d[0])
    for i in range(1, n):
        m = [[B[i][r][c] - (A[i][r][0]*Cp[i-1][0][c] + A[i][r][1]*Cp[i-1][1][c]) for c in range(2)] for r in range(2)]
        mi = inv(m)
        Cp[i] = mm(mi, C[i])
        av = mv(A[i], dp[i-1])
        dp[i] = mv(mi, [d[i][0] - av[0], d[i][1] - av[1]])
    x = [None] * n
    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        cv = mv(Cp[i], x[i+1])
        x[i] = [dp[i][0] - cv[0], dp[i][1] - cv[1]]
    return x


def rates(m, p, T, P):
    prod = sum(m.V[i] * m.chim[i] * p["Da_p"] * P[i] / (p["kappa"] + P[i]) for i in range(m.N))
    cons = sum(m.V[i] * m.chic[i] * p["Da_c"] * T[i] for i in range(m.N))
    return prod, cons


def step(m, p, T, P, dt, tol=1e-12, maxit=500, newton=False):
    N, dl = m.N, p["delta"]
    Pk = P[:]
    for it in range(maxit):
        A, B, C, d = [], [], [], []
        for i in range(N):
            V = m.V[i]
            kc = V * m.chic[i] * p["Da_c"]
            kap = p["kappa"]
            if newton:  # Newton linearisation of Da_p*P/(kappa+P) about Pk (not positivity-safe)
                f = Pk[i] / (kap + Pk[i]); fp = kap / (kap + Pk[i]) ** 2
            else:       # Picard: Da_p*P^{new}/(kappa+Pk); M-matrix, keeps T, P >= 0
                fp = 1.0 / (kap + Pk[i]); f = fp * Pk[i]
            kp = V * m.chim[i] * p["Da_p"] * fp
            s0 = V * m.chim[i] * p["Da_p"] * (f - fp * Pk[i])
            gl, gr = m.Gl[i], m.Gr[i]
            B.append([[V / dt + gl + gr + kc, -kp], [-kc, V / dt + dl * (gl + gr) + kp]])
            A.append([[-gl, 0.0], [0.0, -dl * gl]])
            C.append([[-gr, 0.0], [0.0, -dl * gr]])
            d.append([V * T[i] / dt + s0, V * P[i] / dt - s0])
        x = solve_block(A, B, C, d)
        Tn, Pn = [v[0] for v in x], [v[1] for v in x]
        err = max(abs(Pn[i] - Pk[i]) for i in range(N))
        Pk = Pn
        if err < tol:
            break
    return Tn, Pn


def total(m, X):
    return sum(m.V[i] * X[i] for i in range(m.N))


def transient(N, p, dt=0.01, nsteps=300, wrong_area=False):
    """Returns max relative adenylate-drift and max relative ATP-budget residual."""
    m = Mesh(N, p, wrong_area)
    T = [p["T0"]] * N; P = [1.0 - p["T0"]] * N
    A0 = total(m, T) + total(m, P)
    T_prev_tot = total(m, T); budget_int = 0.0; T_init = T_prev_tot
    drift = bud = 0.0
    for n in range(nsteps):
        T, P = step(m, p, T, P, dt)
        pr, co = rates(m, p, T, P)              # backward Euler: rates at t^{n+1}
        budget_int += dt * (pr - co)
        drift = max(drift, abs(total(m, T) + total(m, P) - A0) / A0)
        bud = max(bud, abs(total(m, T) - T_init - budget_int) / A0)
    return dict(N=N, dt=dt, nsteps=nsteps, t_end=dt * nsteps, max_rel_adenylate_drift=drift,
                max_rel_ATP_budget_residual=bud, final_T_mean=total(m, T) / total(m, [1.0] * N))


def steady(N, p, tol=1e-9, dt0=1e-2, dtmax=1e2, maxsteps=400, newton=False):
    """Pseudo-transient continuation: dt grows x2 per step from dt0 to dtmax.
    dt0 = dtmax = 1e2 (a single huge first step) is the unsafe variant kept as a
    counterexample: Newton can overshoot to P < 0 and converge to a spurious,
    still perfectly conservative, steady state."""
    m = Mesh(N, p)
    T = [p["T0"]] * N; P = [1.0 - p["T0"]] * N
    dt = dt0
    for k in range(maxsteps):
        Tn, Pn = step(m, p, T, P, dt, newton=newton)
        ch = max(abs(Tn[i] - T[i]) for i in range(N))
        T, P = Tn, Pn
        if ch < tol and dt >= dtmax:
            break
        dt = min(2 * dt, dtmax)
    return m, T, P


def metrics(m, p, T, P):
    Vm = sum(m.V[i] * m.chim[i] for i in range(m.N))
    Vc = sum(m.V[i] * m.chic[i] for i in range(m.N))
    Tm = sum(m.V[i] * m.chim[i] * T[i] for i in range(m.N)) / Vm
    Tc = sum(m.V[i] * m.chic[i] * T[i] for i in range(m.N)) / Vc
    Pc = sum(m.V[i] * m.chic[i] * P[i] for i in range(m.N)) / Vc
    pr, co = rates(m, p, T, P)
    Tw = well_mixed_T(p, Vm, Vc)
    return dict(T_mito=Tm, T_cons=Tc, ATP_ADP_cons=Tc / Pc, gradient_index=1.0 - Tc / Tm,
                flux=co, effectiveness=co / (p["Da_c"] * Vc * Tw),
                steady_balance_rel=abs(pr - co) / max(co, 1e-300), T_wellmixed=Tw,
                min_conc=min(min(T), min(P)))


def well_mixed_T(p, Vm, Vc):
    """Uniform-concentration (infinitely fast diffusion) steady state, T+P = 1."""
    f = lambda T: p["Da_p"] * Vm * (1 - T) / (p["kappa"] + 1 - T) - p["Da_c"] * Vc * T
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0: lo = mid
        else: hi = mid
    return 0.5 * (lo + hi)


def run_metrics(N, p, **kw):
    m, T, P = steady(N, p, **kw)
    return metrics(m, p, T, P)
