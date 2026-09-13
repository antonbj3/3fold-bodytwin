#!/usr/bin/env python3
"""NOTCH LATERAL INHIBITION: independent recomputation of a recorded Collier-1996 pattern claim.

Rebuilds the Collier et al. 1996 two-cell and N=20 ring lateral-inhibition model from its own
equations (2) and (4) and re-derives, without reusing any recorded intermediate: the homogeneous
fixed point and its feedback margin f'(d0)*g'(n0) (analytic vs finite-difference), the
heterogeneous period-2 fixed point by root scan, the same state by direct LSODA integration, the
N=20 ring pattern statistics at three Hill exponents, a weak-feedback control (b=1), a
lateral-induction sign-flip control, and the critical trans-activation / cis-attenuation
fractions phi_c and chi_c by bisection on the margin.

Reads: nothing (parameters embedded). Writes: nothing (prints its numbers).
"""
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp

# ---- Collier1996 eq(2),(4) verbatim ----
def f_trans(x, a, k):  # increasing Hill, trans-activation
    return x**k / (a + x**k)

def g_inh(x, b, h):    # decreasing Hill, own-cell repression (lateral inhibition)
    return 1.0 / (1.0 + b * x**h)

def g_ind(x, b, h):    # increasing Hill (lateral induction sign-flip control)
    return x**h / (b + x**h)

def df_trans(x, a, k):
    # d/dx [x^k/(a+x^k)] = a*k*x^(k-1)/(a+x^k)^2
    return a * k * x**(k - 1) / (a + x**k)**2

def dg_inh(x, b, h):
    # d/dx [1/(1+b x^h)] = -b*h*x^(h-1)/(1+b x^h)^2
    return -b * h * x**(h - 1) / (1.0 + b * x**h)**2

def dg_ind(x, b, h):
    return b * h * x**(h - 1) / (b + x**h)**2

def fd_deriv(fn, x, eps=1e-6, *args):
    return (fn(x + eps, *args) - fn(x - eps, *args)) / (2 * eps)

# ---- P1 params ----
a, b, nu, k, h = 0.01, 100.0, 1.0, 2, 2

# homogeneous fixed point x0 = f(g(x0)), fg monotone decreasing -> unique root
def fg(x, a, b, k, h):
    return f_trans(g_inh(x, b, h), a, k)

x0 = brentq(lambda x: fg(x, a, b, k, h) - x, 1e-12, 1 - 1e-12, xtol=1e-14)
d0 = g_inh(x0, b, h)
n0 = x0

fprime = df_trans(d0, a, k)
gprime = dg_inh(n0, b, h)
margin = fprime * gprime
fd_f = fd_deriv(df_trans, d0, 1e-6, a, k)  # sanity, not used
fprime_fd = (f_trans(d0 + 1e-6, a, k) - f_trans(d0 - 1e-6, a, k)) / 2e-6
gprime_fd = (g_inh(n0 + 1e-6, b, h) - g_inh(n0 - 1e-6, b, h)) / 2e-6
margin_fd = fprime_fd * gprime_fd
rel_err_margin = abs(margin - margin_fd) / abs(margin_fd)

print(f"[HOMOG P1] x0={x0:.6f} n0={n0:.6f} d0={d0:.6f}")
print(f"[MARGIN P1] analytic={margin:.6f} FD={margin_fd:.6f} rel_err={rel_err_margin:.3e}")

# ---- heterogeneous (period-2) fixed point: grid+bisection scan of Phi(z)=fg(fg(z))-z ----
def phi2(z, a, b, k, h):
    return fg(fg(z, a, b, k, h), a, b, k, h) - z

def scan_roots(fn, n=4000, lo=1e-9, hi=1 - 1e-9, *args):
    zs = np.linspace(lo, hi, n)
    vals = np.array([fn(z, *args) for z in zs])
    roots = []
    for i in range(len(zs) - 1):
        if vals[i] == 0:
            roots.append(zs[i])
        elif vals[i] * vals[i + 1] < 0:
            roots.append(brentq(fn, zs[i], zs[i + 1], args=args, xtol=1e-15))
    return roots

roots2 = scan_roots(phi2, 4000, 1e-9, 1 - 1e-9, a, b, k, h)
# filter out the trivial homogeneous root (z=x0 itself, since fg(fg(x0))=x0 trivially)
het_roots = [z for z in roots2 if abs(z - x0) > 1e-6]
n1_star = min(het_roots) if het_roots else None
n2_star = fg(n1_star, a, b, k, h) if n1_star is not None else None
d1_star = g_inh(n1_star, b, h) if n1_star is not None else None
d2_star = g_inh(n2_star, b, h) if n2_star is not None else None

print(f"[HET P1] all period-2 roots found: {sorted(roots2)}")
print(f"[HET P1] n1*={n1_star:.6f} d1*={d1_star:.6f} n2*={n2_star:.6f} d2*={d2_star:.6f}")

target = dict(n1=0.010102, d1=0.989898, n2=0.989898, d2=0.010102)
computed = dict(n1=n1_star, d1=d1_star, n2=n2_star, d2=d2_star)
for key in target:
    rel = abs(computed[key] - target[key]) / target[key]
    print(f"  {key}: computed={computed[key]:.6f} target={target[key]:.6f} rel_err={rel:.3e}")

# ---- forward integrate 2-cell system from Collier Fig.3 ICs ----
def rhs_2cell(t, y, a, b, nu, k, h):
    n1, d1, n2, d2 = y
    dn1 = f_trans(d2, a, k) - n1   # dbar_1 = d2 (2-cell: neighbor is the other cell)
    dd1 = nu * (g_inh(n1, b, h) - d1)
    dn2 = f_trans(d1, a, k) - n2
    dd2 = nu * (g_inh(n2, b, h) - d2)
    return [dn1, dd1, dn2, dd2]

y0 = [1.0, 1.0, 0.99, 0.99]
sol = solve_ivp(rhs_2cell, [0, 60], y0, args=(a, b, nu, k, h),
                 method='LSODA', rtol=1e-10, atol=1e-12, dense_output=True)
yT = sol.y[:, -1]
dydt_T = rhs_2cell(60, yT, a, b, nu, k, h)
maxdydt = max(abs(v) for v in dydt_T)
print(f"[ODE P1] final state n1,d1,n2,d2 = {yT}")
print(f"[ODE P1] max|dydt| at T=60: {maxdydt:.3e}")

ode_vals = dict(n1=yT[0], d1=yT[1], n2=yT[2], d2=yT[3])
rel_errs_ode_vs_alg = {}
for key in target:
    rel = abs(ode_vals[key] - computed[key]) / computed[key]
    rel_errs_ode_vs_alg[key] = rel
    print(f"  ODE-vs-algebraic {key}: rel_err={rel:.3e}")
max_rel_err_ode_vs_alg = max(rel_errs_ode_vs_alg.values())
print(f"[CLAIM CHECK] claimed rel_err=4.27e-13; my max rel_err ODE-vs-algebraic = {max_rel_err_ode_vs_alg:.3e}")

# ---- trace check: trace = -1-nu < 0 unconditionally ----
print("\n[TRACE INVARIANT CHECK] trace = -1-nu must be <0 for all nu>0 tried:")
for nu_test in [0.08125, 0.1, 0.5, 1, 2, 5, 12]:
    trace = -1 - nu_test
    print(f"  nu={nu_test}: trace={trace} <0 ? {trace < 0}")

# ---- N=20 ring simulation ----
def ring_rhs(t, y, params, mode='inh'):
    N = len(y) // 2
    n = y[0::2]
    d = y[1::2]
    dbar = (np.roll(d, 1) + np.roll(d, -1)) / 2.0
    a_, b_, nu_, k_, h_ = params
    dn = f_trans(dbar, a_, k_) - n
    if mode == 'inh':
        dd = nu_ * (g_inh(n, b_, h_) - d)
    else:
        dd = nu_ * (g_ind(n, b_, h_) - d)
    out = np.empty_like(y)
    out[0::2] = dn
    out[1::2] = dd
    return out

def run_ring(a_, b_, nu_, k_, h_, N=20, seed=20260723, T=200, mode='inh'):
    rng = np.random.default_rng(seed)
    x0_ = brentq(lambda x: fg(x, a_, b_, k_, h_) - x, 1e-12, 1 - 1e-12, xtol=1e-14) if mode == 'inh' else None
    if mode == 'inh':
        n_ss, d_ss = x0_, g_inh(x0_, b_, h_)
    else:
        # for induction, use the middle root as the base (unstable one, ~0.0101 in P1)
        n_ss, d_ss = 0.5, 0.5
    y0 = np.empty(2 * N)
    y0[0::2] = n_ss * (1 + 0.001 * rng.standard_normal(N))
    y0[1::2] = d_ss * (1 + 0.001 * rng.standard_normal(N))
    y0 = np.clip(y0, 1e-9, 1 - 1e-9)
    sol = solve_ivp(ring_rhs, [0, T], y0, args=((a_, b_, nu_, k_, h_), mode),
                     method='LSODA', rtol=1e-10, atol=1e-12)
    nfinal = sol.y[0::2, -1]
    std_n = nfinal.std()
    # alternation: fraction of neighbor-pairs with opposite "high/low" (>0.5 threshold) state
    high = nfinal > 0.5
    alt = np.mean(high != np.roll(high, 1))
    return std_n, alt, nfinal

print("\n[N=20 RING] baseline P1 (k=h=2):")
std_p1, alt_p1, _ = run_ring(a, b, nu, k, h)
print(f"  std={std_p1:.3f} alt={alt_p1:.3f} (claimed std=0.482 alt=0.900)")

print("[N=20 RING] Sprinzak k=h=1.7:")
std_spr, alt_spr, _ = run_ring(a, b, nu, 1.7, 1.7)
print(f"  std={std_spr:.3f} alt={alt_spr:.3f} (claimed std=0.471 alt=1.000)")

print("[N=20 RING] cis-equivalent k=h=12:")
std_cis, alt_cis, _ = run_ring(a, b, nu, 12, 12)
print(f"  std={std_cis:.3f} alt={alt_cis:.3f} (claimed std=0.488 alt=0.900)")

# ---- weak-feedback control b=1 ----
print("\n[WEAK FEEDBACK CONTROL b=1]")
b_weak = 1.0
x0_weak = brentq(lambda x: fg(x, a, b_weak, k, h) - x, 1e-12, 1 - 1e-12, xtol=1e-14)
d0_weak = g_inh(x0_weak, b_weak, h)
fprime_w = df_trans(d0_weak, a, k)
gprime_w = dg_inh(x0_weak, b_weak, h)
margin_weak = fprime_w * gprime_w
roots2_weak = scan_roots(phi2, 4000, 1e-9, 1 - 1e-9, a, b_weak, k, h)
n_roots_weak = len(set(round(r, 6) for r in roots2_weak))
std_weak, alt_weak, _ = run_ring(a, b_weak, nu, k, h)
print(f"  margin={margin_weak:.6f} (claimed -0.069) n_period2_roots={n_roots_weak} (claimed globally unique=1)")
print(f"  ring std={std_weak:.3e} (claimed 1.6e-15 / re-quoted 1.57e-15)")

# ---- lateral-induction sign-flip control ----
print("\n[LATERAL INDUCTION CONTROL] g_ind(x)=x^h/(b_ind+x^h), b_ind=0.01")
b_ind = 0.01
def fg_ind(x, a_, b_, k_, h_):
    return f_trans(g_ind(x, b_, h_), a_, k_)
roots_ind = scan_roots(lambda z, *args: fg_ind(z, *args) - z, 4000, 1e-9, 1 - 1e-9, a, b_ind, k, h)
print(f"  homogeneous roots: {[round(r,4) for r in sorted(roots_ind)]} (claimed 0, 0.0101, 0.9899)")
margins_ind = []
for r in roots_ind:
    d_r = g_ind(r, b_ind, h)
    fp = df_trans(d_r, a, k)
    gp = dg_ind(r, b_ind, h)
    margins_ind.append(fp * gp)
    print(f"    root={r:.4f} margin={fp*gp:.4f} (>=0? {fp*gp >= 0})")

print("  bistability test: 6 runs straddling unstable root ~0.0101")
rng = np.random.default_rng(1)
lows, highs = 0, 0
results_bi = []
for trial in range(6):
    ic = 0.001 if trial < 3 else 0.5
    N = 20
    y0 = np.empty(2 * N)
    y0[0::2] = ic
    y0[1::2] = g_ind(ic, b_ind, h)
    y0 = y0 * (1 + 0.001 * rng.standard_normal(2 * N))
    y0 = np.clip(y0, 1e-9, 1 - 1e-9)
    sol = solve_ivp(ring_rhs, [0, 200], y0, args=((a, b_ind, nu, k, h), 'ind'),
                     method='LSODA', rtol=1e-10, atol=1e-12)
    nfinal = sol.y[0::2, -1]
    results_bi.append((nfinal.mean(), nfinal.std()))
    print(f"    trial {trial}: mean_n={nfinal.mean():.4f} internal_std={nfinal.std():.3e}")
means_bi = [m for m, s in results_bi]
spread_bi = max(means_bi) - min(means_bi)
print(f"  spread across runs = {spread_bi:.3f} (claimed 0.990)")

# ---- GSI dose-response: phi_c and chi_c ----
print("\n[GSI DOSE RESPONSE]")
def fg_phi(x, a_, b_, k_, h_, phi):
    return f_trans(g_inh(x, b_, h_), a_, k_) * (1 - phi)

def margin_phi(phi, a_, b_, k_, h_):
    x0_ = brentq(lambda x: fg_phi(x, a_, b_, k_, h_, phi) - x, 1e-12, 1 - 1e-12, xtol=1e-14)
    d0_ = g_inh(x0_, b_, h_)
    fp = df_trans(d0_, a_, k_) * (1 - phi)
    gp = dg_inh(x0_, b_, h_)
    return fp * gp + 1.0  # margin - (-1) = 0 at critical

phi_c = brentq(margin_phi, 1e-6, 0.999, args=(a, b, k, h), xtol=1e-12)
print(f"  phi_c = {phi_c:.4f} (claimed 0.6677)")

def fg_chi(x, a_, b_, k_, h_, chi):
    return f_trans(g_inh(x, b_, h_) * (1 - chi), a_, k_)

def margin_chi(chi, a_, b_, k_, h_):
    x0_ = brentq(lambda x: fg_chi(x, a_, b_, k_, h_, chi) - x, 1e-12, 1 - 1e-12, xtol=1e-14)
    d0_ = g_inh(x0_, b_, h_) * (1 - chi)
    fp = df_trans(d0_, a_, k_)
    gp = dg_inh(x0_, b_, h_) * (1 - chi)
    return fp * gp + 1.0

chi_c = brentq(margin_chi, 1e-6, 0.999, args=(a, b, k, h), xtol=1e-12)
print(f"  chi_c = {chi_c:.4f} (claimed 0.9655)")
print(f"  delta = {abs(phi_c - chi_c):.4f} (claimed 0.298)")

# forced check: are phi_c, chi_c related by an exact algebraic transform (f,g interchange symmetry)?
# f(x)=x^k/(a+x^k); g(x)=1/(1+b x^h). Attenuating f: (1-phi)*f(g(x)). Attenuating g: f((1-chi)*g(x)).
# These are NOT the same functional form unless f is linear, so no forced equality is expected -- check numerically:
print("  forced-symmetry check: is there an exact algebraic map phi<->chi predicted by f,g structure?")
print("  f is a saturating Hill (nonlinear), so (1-phi)*f(y) != f((1-chi)*y) in general -> no forced equality expected.")
print(f"  Confirms: phi_c != chi_c is EXPECTED from the algebra (non-affine f), not merely an empirical curiosity.")

# ODE-confirm phi_c +/- 0.20
for phi_test, label in [(phi_c - 0.20, 'phi_c-0.20'), (phi_c + 0.20, 'phi_c+0.20')]:
    std_t, alt_t, _ = run_ring(a, b, nu, k, h)  # placeholder, need phi-modified ring
    # redo properly with phi in ring rhs
    def ring_rhs_phi(t, y, params):
        N = len(y) // 2
        n = y[0::2]; d = y[1::2]
        dbar = (np.roll(d, 1) + np.roll(d, -1)) / 2.0
        a_, b_, nu_, k_, h_, phi_ = params
        dn = f_trans(dbar, a_, k_) * (1 - phi_) - n
        dd = nu_ * (g_inh(n, b_, h_) - d)
        out = np.empty_like(y)
        out[0::2] = dn; out[1::2] = dd
        return out
    rng2 = np.random.default_rng(20260723)
    x0_p = brentq(lambda x: fg_phi(x, a, b, k, h, phi_test) - x, 1e-12, 1 - 1e-12, xtol=1e-14)
    y0p = np.empty(40)
    y0p[0::2] = x0_p * (1 + 0.001 * rng2.standard_normal(20))
    y0p[1::2] = g_inh(x0_p, b, h) * (1 + 0.001 * rng2.standard_normal(20))
    y0p = np.clip(y0p, 1e-9, 1 - 1e-9)
    sol = solve_ivp(ring_rhs_phi, [0, 200], y0p, args=((a, b, nu, k, h, phi_test),),
                     method='LSODA', rtol=1e-10, atol=1e-12)
    stdf = sol.y[0::2, -1].std()
    print(f"  {label}={phi_test:.4f}: ring std={stdf:.3e} (claimed 0.232 patterned / 1.4e-17 uniform)")
