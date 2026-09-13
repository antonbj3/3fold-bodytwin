#!/usr/bin/env python3
"""somite_hes7_lewis_zeiser_dde_rebuild.py., resolves
MODEL-SOMITE-SEGMENTATION-CLOCK (cluster: developmental -- the embryonic segmentation clock that
sets somite/vertebra periodicity).

STATUS BEFORE THIS SCRIPT: node's cert_design.verify said "SEED-DESIGN ... cited literature, not
executed". The node names a fully-specifying delay-differential-equation (DDE) source (Zeiser et al.
2006, Theor Biol Med Model 3:11, PMID16504083, built on the Lewis 2003/Monk 2003 Hes7 autorepression
form) and a held-out anchor (Masamizu et al. 2006 PNAS, PMID16432209, in-utero period ~120min), and
reports a numeric result (period=120.94min) but no script in the project implements the DDE (grep-
checked: no scripts/*/*.py references Zeiser, Hes7-DDE, or this node id before this file).

SOURCE (live-fetched when this cell was written, PMC1402261, Zeiser et al. 2006 Methods + reference-case
parameters, itself built on Lewis 2003 Curr Biol 13:1398 / Monk 2003 Curr Biol 13:1409):
  dm/dt = -c*m(t) + k * f_h( p(t - Tm) )      [Hes7 mRNA, molecules/cell]
  dp/dt = -b*p(t) + a * m(t - Tp)             [Hes7 protein, molecules/cell]
  f_h(p) = 1 / (1 + (p/H)^h)                  [Hill repression]
Reference-case parameters (live-fetched): protein half-life tau_p=20min -> b=ln2/20 /min; mRNA
half-life tau_m=3min -> c=ln2/3 /min; translation rate a=4.5 molecules/mRNA/min; basal transcription
k=4.5 mRNA/min; Hill threshold H=40 protein molecules/cell; Hill coefficient h=2 (reference case);
transcriptional delay Tm=30min, translational delay Tp=7min (Tm+Tp=37min total, matches the paper's
own reported "for a protein half-life of tau_p=20min and h=2, the system shows undamped oscillations
with a period of about 120min").

QUESTION (pre-registered before running): does a from-scratch numerical DDE integration (explicit
method-of-steps RK4 with linear-interpolated history, NOT copying the node's stored 120.94min
figure) at these exact live-fetched parameters produce a stable limit-cycle oscillation whose period
lands near the Masamizu et al. 2006 in-utero anchor (~120min, decorrelated: Kageyama-lab live-imaging
wet lab vs this Zeiser-parameter numeric DDE)? Cross-checked against an independent ANALYTIC
characteristic-equation (linearized-DDE, Hopf-adjacent) root-find for the same frequency.

GATE (pre-registered):
  G1 numeric period in [100,140]min (the node's honest_gaps disclosed bracket: the SAME anchor
     paper's in-utero-vs-ex-vivo condition spread 120-160min, used here as a symmetric +/-20min
     sanity band around the in-vivo point, not a fabricated sigma).
  G2 analytic-vs-numeric cross-check: the linearized characteristic-equation frequency must agree
     with the full-nonlinear numeric limit-cycle frequency to within 15% (two independent methods on
     the SAME parameters -- a correctness check on this rebuild, not a biological claim).
VOID FLOOR (pre-registered, two variants, per the discipline "measure your void floor before
  trusting it"):
  V1 zero delay (Tm=Tp=0): the same nonlinear negative-feedback ODE (no delay) should NOT sustain
     oscillation (relax to a stable fixed point) -- delay is THE mechanism, not incidental.
  V2 weak cooperativity (h=1, linear-ish repression, everything else unchanged): should also fail to
     sustain a robust limit cycle at these rate constants (Hill steepness is load-bearing, not free).
  Both floors are measured, not assumed; if either fails to discriminate it is reported, not hidden.
"""
import numpy as np
from scipy.optimize import brentq

# ---- live-fetched reference-case parameters (Zeiser et al. 2006, PMC1402261) ----
TAU_P = 20.0
TAU_M = 3.0
B = np.log(2) / TAU_P
C = np.log(2) / TAU_M
A_TRANSL = 4.5
K_TRANSCR = 4.5
H_THRESH = 40.0
HILL_H = 2.0
T_M_DELAY = 30.0
T_P_DELAY = 7.0

ANCHOR_MIN = 120.0
ANCHOR_BAND = (100.0, 140.0)  # +/-20min around in-vivo point = the paper's condition spread

def f_hill(p, H=H_THRESH, h=HILL_H):
    return 1.0 / (1.0 + (max(p, 0.0) / H) ** h)

def steady_state(H=H_THRESH, h=HILL_H, a=A_TRANSL, k=K_TRANSCR, b=B, c=C):
    # m* = k*f(p*)/c ; p* = a*m*/b  => solve self-consistent fixed point for p*
    def g(p):
        m = k * f_hill(p, H, h) / c
        return a * m / b - p
    p_star = brentq(g, 0.0, 10 * H + 1000)
    m_star = k * f_hill(p_star, H, h) / c
    return m_star, p_star

def integrate_dde(Tm, Tp, H=H_THRESH, h=HILL_H, a=A_TRANSL, k=K_TRANSCR, b=B, c=C,
                   t_end=2000.0, dt=0.02, perturb=1.05):
    """Explicit method-of-steps RK4 for the 2-delay Hes7 DDE, fixed-step, linear-interp history."""
    m_star, p_star = steady_state(H, h, a, k, b, c)
    n_steps = int(t_end / dt) + 1
    t = np.linspace(0.0, t_end, n_steps)
    m = np.empty(n_steps)
    p = np.empty(n_steps)
    m[0] = m_star * perturb
    p[0] = p_star

    max_delay = max(Tm, Tp)
    n_hist = max(1, int(round(max_delay / dt)))

    def hist_val(arr, i, delay_steps):
        # arr defined for indices 0..i; delayed index i-delay_steps, clamp to history (t<0 -> ICs)
        idx = i - delay_steps
        if idx < 0:
            return arr[0]
        return arr[idx]

    dTm = int(round(Tm / dt))
    dTp = int(round(Tp / dt))

    def deriv(i_now, m_now, p_now):
        p_delay = hist_val(p, i_now, dTm) if i_now - dTm >= 0 else p[0]
        m_delay = hist_val(m, i_now, dTp) if i_now - dTp >= 0 else m[0]
        dm = -c * m_now + k * f_hill(p_delay, H, h)
        dp = -b * p_now + a * m_delay
        return dm, dp

    for i in range(n_steps - 1):
        # RK4 with delayed terms evaluated at the (frozen, already-known) history -- standard
        # method-of-steps: delayed args at t_i, t_i+dt/2, t_i+dt all lie in the ALREADY-COMPUTED
        # past (since Tm,Tp >> dt), so a plain RK4 on the non-delayed state is valid here.
        dm1, dp1 = deriv(i, m[i], p[i])
        dm2, dp2 = deriv(i, m[i] + 0.5 * dt * dm1, p[i] + 0.5 * dt * dp1)
        dm3, dp3 = deriv(i, m[i] + 0.5 * dt * dm2, p[i] + 0.5 * dt * dp2)
        dm4, dp4 = deriv(i, m[i] + dt * dm3, p[i] + dt * dp3)
        m[i + 1] = m[i] + (dt / 6.0) * (dm1 + 2 * dm2 + 2 * dm3 + dm4)
        p[i + 1] = p[i] + (dt / 6.0) * (dp1 + 2 * dp2 + 2 * dp3 + dp4)
        if m[i + 1] < 0:
            m[i + 1] = 0.0
        if p[i + 1] < 0:
            p[i + 1] = 0.0
    return t, m, p

def measure_period(t, p, transient_frac=0.5):
    """Peak-to-peak period in the late (post-transient) portion of the trace."""
    n0 = int(len(t) * transient_frac)
    tail_t, tail_p = t[n0:], p[n0:]
    peaks = []
    for i in range(1, len(tail_p) - 1):
        if tail_p[i] > tail_p[i - 1] and tail_p[i] >= tail_p[i + 1]:
            peaks.append(tail_t[i])
    if len(peaks) < 3:
        return None, 0.0
    gaps = np.diff(peaks)
    return float(np.mean(gaps)), float(np.std(gaps))

def analytic_hopf_frequency(Tm, Tp, H=H_THRESH, h=HILL_H, a=A_TRANSL, k=K_TRANSCR, b=B, c=C,
                             omega0=2 * np.pi / 120.0):
    """Complex Newton root-find of the linearized-DDE characteristic equation
    (lambda+c)(lambda+b) - k*a*f'(p*)*exp(-lambda*(Tm+Tp)) = 0
    near lambda = i*omega0, independent cross-check of the numeric limit-cycle frequency."""
    m_star, p_star = steady_state(H, h, a, k, b, c)
    # f_h(p) = (1+(p/H)^h)^-1 ; f_h'(p) = -h/H * (p/H)^(h-1) * (1+(p/H)^h)^-2
    x = p_star / H
    fprime = -h / H * x ** (h - 1) * (1.0 + x ** h) ** (-2)
    T = Tm + Tp

    def F(lam):
        return (lam + c) * (lam + b) - k * a * fprime * np.exp(-lam * T)

    lam = 1j * omega0
    for _ in range(200):
        eps = 1e-6
        Fval = F(lam)
        Fp = (F(lam + eps) - F(lam - eps)) / (2 * eps)
        step = Fval / Fp
        lam = lam - step
        if abs(step) < 1e-12:
            break
    return lam

def main():
    print("=== somite_hes7_lewis_zeiser_dde_rebuild: MODEL-SOMITE-SEGMENTATION-CLOCK ===")
    print(f"params: b={B:.5f}/min c={C:.5f}/min a={A_TRANSL} k={K_TRANSCR} H={H_THRESH} "
          f"h={HILL_H} Tm={T_M_DELAY}min Tp={T_P_DELAY}min")

    t, m, p = integrate_dde(T_M_DELAY, T_P_DELAY, t_end=2400.0, dt=0.02)
    period_mean, period_std = measure_period(t, p, transient_frac=0.5)
    print(f"\nNUMERIC DDE (method-of-steps RK4, fresh integration): period = {period_mean:.3f} "
          f"+/- {period_std:.3f} min (limit-cycle regularity, n peaks-derived)")

    lam = analytic_hopf_frequency(T_M_DELAY, T_P_DELAY)
    period_analytic = 2 * np.pi / lam.imag if lam.imag != 0 else float("nan")
    print(f"ANALYTIC characteristic-equation root: lambda={lam:.6f}  "
          f"(Re={lam.real:.5f}/min -> {'growing/unstable fixed pt (expected, feeds limit cycle)' if lam.real>0 else 'damped'}), "
          f"period_analytic={period_analytic:.3f} min")

    cross_check_pct = abs(period_mean - period_analytic) / period_analytic * 100
    g2_pass = cross_check_pct < 15.0
    print(f"G2 analytic-vs-numeric cross-check: {cross_check_pct:.2f}% difference "
          f"(gate <15%) -> {'PASS' if g2_pass else 'FAIL'}")

    g1_pass = ANCHOR_BAND[0] <= period_mean <= ANCHOR_BAND[1]
    print(f"\nAnchor: Masamizu et al. 2006 PNAS (PMID16432209) in-utero period ~120min "
          f"(band {ANCHOR_BAND} min, the paper's in-utero-vs-ex-vivo 120-160min condition "
          f"spread used as the honest uncertainty proxy).")
    print(f"G1 anchor: numeric period {period_mean:.3f}min in band? -> {'PASS' if g1_pass else 'FAIL'} "
          f"(residual = {period_mean - ANCHOR_MIN:+.3f} min, {100*(period_mean-ANCHOR_MIN)/ANCHOR_MIN:+.2f}%)")

    # --- void floor V1: zero delay ---
    t0, m0, p0 = integrate_dde(Tm=0.001, Tp=0.001, t_end=1200.0, dt=0.02)
    per0, _ = measure_period(t0, p0, transient_frac=0.5)
    amp0 = np.std(p0[int(len(p0) * 0.7):])
    ss_m, ss_p = steady_state(H_THRESH, HILL_H)
    print(f"\nVOID FLOOR V1 (zero delay, Tm=Tp~0): late-window p std={amp0:.4f} "
          f"(steady-state p*={ss_p:.4f}); period_detected={per0} "
          f"-> {'CORRECTLY NO SUSTAINED OSCILLATION' if (per0 is None or amp0 < 0.01*ss_p) else 'STILL OSCILLATES (floor does not discriminate)'}")

    # --- void floor V2: weak cooperativity h=1 ---
    t1, m1, p1 = integrate_dde(T_M_DELAY, T_P_DELAY, H=H_THRESH, h=1.0, t_end=2400.0, dt=0.02)
    per1, per1sd = measure_period(t1, p1, transient_frac=0.5)
    amp1 = np.std(p1[int(len(p1) * 0.7):])
    ss_m1, ss_p1 = steady_state(H_THRESH, 1.0)
    print(f"VOID FLOOR V2 (weak cooperativity h=1, delays unchanged): late-window p std={amp1:.4f} "
          f"(steady-state p*={ss_p1:.4f}); period_detected={per1} "
          f"-> {'CORRECTLY NO SUSTAINED OSCILLATION (disclosed)' if (per1 is None or amp1 < 0.01*ss_p1) else 'STILL OSCILLATES -- floor does NOT discriminate, reported informationally not gated'}")

    verdict = "CONFIRMED" if (g1_pass and g2_pass) else "PARTIAL"
    print(f"\nVERDICT: {verdict}")
    print(f"DATAPOINTS: period_numeric_min={period_mean:.3f}, period_numeric_sd={period_std:.3f}, "
          f"period_analytic_min={period_analytic:.3f}, cross_check_pct={cross_check_pct:.2f}, "
          f"anchor_min={ANCHOR_MIN}, residual_min={period_mean-ANCHOR_MIN:+.3f}, "
          f"void_v1_amp={amp0:.4f}, void_v1_ss_p={ss_p:.4f}, "
          f"void_v2_amp={amp1:.4f}, void_v2_ss_p={ss_p1:.4f}")

if __name__ == "__main__":
    main()
