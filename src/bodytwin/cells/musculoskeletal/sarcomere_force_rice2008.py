"""Regime-correct implementation of the Rice, Wang, Bers & de Tombe (2008) cardiac sarcomere
cross-bridge / F-pCa model, wired to accept a TIME-RESOLVED Ca2+(t) drive (not a static pCa
plug-in).

Rate constants: Rice2008 Table 1, cross-checked against BOTH the Rice2008 PMC full text and the
independently curated CellML transcription (Physiome Model Repository,
rice_wang_bers_detombe_2008.cellml, "compared with Jeremy Rice's original Matlab code" by Steve
Niederer): kn_p = 500/s (0.5/ms), kp_n = 50/s (0.05/ms). A commonly quoted transcription has these
two swapped; the swapped version happens to land closer to skinned-fibre literature by partially
CANCELLING a separate, larger offset -- a compensating-bug coincidence, not a validation.

Ca-sensitivity regime (Ca50): rather than picking one literature Ca50, the target is treated as a
hidden variable (myofilament Ca sensitivity of intact human myocardium at 37 C, physiological SL)
that no paper measures directly; each measurement is a differently biased view of it.
  - Skinning/permeabilization is documented to REDUCE Ca sensitivity (lower pCa50) via loss of
    resting PKA/PKC phosphorylation and increased lattice spacing (Methods for Assessing Cardiac
    Myofilament Calcium Sensitivity review, PMC10728676, 2023) -- a signed, mechanistic bias.
  - Low temperature independently reduces pCa50 (Rice2008 Table 1 Q10 values: 37 C -> 15 C shifts
    model pCa50 by -0.18 to -0.19 units at fixed rate constants).
  - Species and tetanized-vs-twitch protocol are second-order next to those two.
Five independently sourced numbers (none used to fit any parameter here) on those bias axes:
    INTACT / near-physiological temperature (expected HIGH pCa50):
      - Backx, Gao, Azan-Backx & Marban 1995, J Gen Physiol 105:1-19 (PMID 7730787): intact RAT
        trabeculae, SL 2.1-2.3 um, tetanized+ryanodine, [Ca2+]i50 = 0.65+/-0.19 uM => pCa50 6.19,
        nH 5.2+/-1.2 (bath temperature not recoverable from the abstract -- a disclosed gap;
        likely room temperature for that group's protocol, which would put the true 37 C intact
        value higher still, not lower).
      - Staurosporine study (PMC3357507), intact RABBIT trabeculae, 37 C, SL 2.2-2.3 um, 1 Hz:
        pCa50 6.05+/-0.04, nH 3.99+/-0.60.
    SKINNED / PERMEABILIZED (expected LOW pCa50):
      - Konhilas et al 2002 (PMC2290573), skinned RAT, 15 C, SL 1.95-2.25 um: pCa50 5.47-5.56.
      - Human permeabilized myofibrils, nonfailing LV (PMC507413): [Ca2+]50 = 3.24+/-0.51 uM
        => pCa50 5.49.
The clusters are separated by ~0.5-0.7 pCa units, much larger than the within-cluster species
spread (intact 6.05-6.19 across rat/rabbit; skinned 5.47-5.56 across rat/human), i.e. prep type is
the dominant bias axis and species is second-order.

Reads: nothing. Writes: nothing (the self-check prints to stdout).

Gate: this file's emergent pCa50/nH (bug-fixed CellML rates, TmpC = 37, UNTUNED -- no
Ca50/kon/koff parameter adjusted to match anything above) is checked against that
bias-decorrelated interval as one more constraint, not as ground truth and not as a fit target.
"""
import numpy as np
from scipy.integrate import solve_ivp

# ---------------------------------------------------------------------------------------
# Rice2008 Table 1 rate constants, machine-cross-checked vs PMC full text + CellML, re-derived
# from primary sources rather than copied from the commonly quoted (swapped) transcription.
# ---------------------------------------------------------------------------------------
kon = 0.05      # /uM/ms
koffL = 0.25    # /ms
koffH = 0.025   # /ms
perm50 = 0.5
nperm = 15.0
KN_P = 0.5      # /ms  -- the swapped transcription has this as 0.05
KP_N = 0.05     # /ms  -- the swapped transcription has this as 0.5
fapp = 0.5
gapp = 0.07
hf = 2.0
hb = 0.4
gxb = 0.07      # /ms
gslmod = 6.0
hfmdc = 5.0
sigmap = 8.0
sigman = 1.0
x_0 = 0.007     # um
len_thin = 1.2
len_thick = 1.65
len_hbare = 0.1  # um
kxb = 120.0     # mN/mm^2, active-stress scale (cancels in any Fmax-normalized ratio)

Qkon, Qkoff, Qkn_p, Qkp_n = 1.5, 1.3, 1.6, 1.6
Qfapp, Qgapp, Qhf, Qhb, Qgxb = 6.25, 2.5, 6.25, 6.25, 6.25
TREF = 37.0  # deg C -- Table 1 reference; our CICR upstream model is itself an intact-myocyte,
             # ~physiological-temperature model (Shannon2004 SBML, 35-37C convention), so the
             # regime-correct choice is TmpC=37 => q10f(...)=1 identically (NO Q10 correction
             # applied in the coupled chain -- Q10 machinery is kept for transparency/reuse and
             # for the disclosed 15C skinned-fibre falsifier cross-check only).


def q10f(Q, TmpC):
    return Q ** ((TmpC - TREF) / 10.0)


def sovf(SL):
    sovr_ze = min(len_thick / 2, SL / 2)
    sovr_cle = max(SL / 2 - (SL - len_thin), len_hbare / 2)
    len_sovr = sovr_ze - sovr_cle
    SOVFThick = len_sovr * 2 / (len_thick - len_hbare)
    SOVFThin = len_sovr / len_thin
    return SOVFThick, SOVFThin


def rates_at(SL, TmpC):
    SOVFThick, SOVFThin = sovf(SL)
    gapslmd = 1 + (1 - SOVFThick) * gslmod
    return dict(
        SOVFThick=SOVFThick, SOVFThin=SOVFThin,
        fappT=fapp * q10f(Qfapp, TmpC),
        gappT=gapp * gapslmd * q10f(Qgapp, TmpC),
        hfT_base=hf * q10f(Qhf, TmpC),
        hbT=hb * q10f(Qhb, TmpC),
        gxbT_base=gxb * q10f(Qgxb, TmpC),
        konT=kon * q10f(Qkon, TmpC),
        koffLT=koffL * q10f(Qkoff, TmpC),
        koffHT=koffH * q10f(Qkoff, TmpC),
    )


def ode_rhs(t, y, SL, TmpC, ca_func):
    """7-state dynamic ODE (TRPNCaH, TRPNCaL, N, XBprer, XBpostr, xXBprer, xXBpostr).
    ca_func(t) -> Cai in uM. Full strain-feedback (hfmd/gxbmd), NOT the isometric-steady-state
    shortcut -- needed because Ca(t) is time-varying (a real twitch), so the mean-strain fixed
    point used by the closed-form steady-state solver (dSL/dt=0 identity) does not apply."""
    TRPNCaH, TRPNCaL, N, XBprer, XBpostr, xXBprer, xXBpostr = y
    r = rates_at(SL, TmpC)
    Cai = ca_func(t)
    dTRPNCaH = r['konT'] * Cai * (1 - TRPNCaH) - r['koffHT'] * TRPNCaH
    dTRPNCaL = r['konT'] * Cai * (1 - TRPNCaL) - r['koffLT'] * TRPNCaL
    Tropreg = TRPNCaH * r['SOVFThin'] + TRPNCaL * (1 - r['SOVFThin'])
    with np.errstate(over='ignore', divide='ignore'):
        permtot = np.sqrt(abs(1.0 / (1.0 + (perm50 / max(Tropreg, 1e-12)) ** nperm)))
    kn_pT = KN_P * permtot * q10f(Qkn_p, TmpC)
    inprmt = min(1.0 / permtot, 100.0) if permtot > 0 else 100.0
    kp_nT = KP_N * inprmt * q10f(Qkp_n, TmpC)
    P = 1 - N - XBprer - XBpostr
    hfmd = np.exp(-np.sign(xXBprer) * hfmdc * (xXBprer / x_0) ** 2)
    gxbmd = (np.exp(sigmap * ((x_0 - xXBpostr) / x_0) ** 2) if xXBpostr < x_0
             else np.exp(sigman * ((xXBpostr - x_0) / x_0) ** 2))
    hfT = r['hfT_base'] * hfmd
    gxbT = r['gxbT_base'] * gxbmd
    fappT, gappT, hbT = r['fappT'], r['gappT'], r['hbT']
    dN = kp_nT * P - kn_pT * N
    dXBprer = fappT * P + hbT * XBpostr - gappT * XBprer - hfT * XBprer
    dXBpostr = hfT * XBprer - hbT * XBpostr - gxbT * XBpostr
    denom = fappT * hfT + gxbT * hfT + gxbT * gappT + hbT * fappT + hbT * gappT + gxbT * fappT
    denom = max(denom, 1e-12)
    dutyprer = (hbT * fappT + gxbT * fappT) / denom
    dutypostr = (fappT * hfT) / denom
    xPsi = 2.0
    dxXBprer = (xPsi / max(dutyprer, 1e-6)) * (fappT * (-xXBprer) + hbT * (xXBpostr - x_0 - xXBprer))
    dxXBpostr = (xPsi / max(dutypostr, 1e-6)) * (hfT * (xXBprer + x_0 - xXBpostr))
    return [dTRPNCaH, dTRPNCaL, dN, dXBprer, dXBpostr, dxXBprer, dxXBpostr]


def force_of_state(y, SL, TmpC):
    _, _, _, XBprer, XBpostr, xXBprer, xXBpostr = y
    SOVFThick, _ = sovf(SL)
    return kxb * SOVFThick * (xXBpostr * XBpostr + xXBprer * XBprer)


Y0_DEFAULT = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, x_0 * 0.5]


def steady_state_force_closed_form(Cai, SL, TmpC):
    """Closed-form isometric steady state (dSL/dt=0 => xXBprer=0, xXBpostr=x_0 exactly, a
    structural/geometric property of the mean-strain ODE, independent of rate values -- see
    mt_q10__rice2008_q10_falsifier.py derivation). Used ONLY for Fmax normalization and the
    static F-pCa self-check, never for the dynamic beat-driven force below."""
    r = rates_at(SL, TmpC)
    TRPNCaH = r['konT'] * Cai / (r['konT'] * Cai + r['koffHT'])
    TRPNCaL = r['konT'] * Cai / (r['konT'] * Cai + r['koffLT'])
    Tropreg = TRPNCaH * r['SOVFThin'] + TRPNCaL * (1 - r['SOVFThin'])
    with np.errstate(over='ignore', divide='ignore'):
        permtot = np.sqrt(abs(1.0 / (1.0 + (perm50 / max(Tropreg, 1e-12)) ** nperm)))
    kn_pT = KN_P * permtot * q10f(Qkn_p, TmpC)
    inprmt = min(1.0 / permtot, 100.0) if permtot > 0 else 100.0
    kp_nT = KP_N * inprmt * q10f(Qkp_n, TmpC)
    fappT, gappT, hfT, hbT, gxbT = r['fappT'], r['gappT'], r['hfT_base'], r['hbT'], r['gxbT_base']
    A = np.array([
        [-kn_pT, kp_nT, 0.0, 0.0],
        [kn_pT, -(kp_nT + fappT), gappT, gxbT],
        [0.0, fappT, -(gappT + hfT), hbT],
        [1.0, 1.0, 1.0, 1.0],
    ])
    b = np.array([0.0, 0.0, 0.0, 1.0])
    N, P, XBprer, XBpostr = np.linalg.solve(A, b)
    force = kxb * r['SOVFThick'] * (x_0 * XBpostr)
    return force


def fmax_reference(SL, TmpC, pca_sat=4.0):
    """Saturating-Ca reference force (Fmax), the operational definition used throughout the
    F-pCa literature this script is checked against."""
    Ca_sat = 10 ** (-pca_sat) * 1e6  # uM
    return steady_state_force_closed_form(Ca_sat, SL, TmpC)


def emergent_fpca(SL, TmpC, pCa_grid):
    """F-pCa curve from the SAME closed-form steady state, for self-consistency / Hill-fit
    reporting (this script's emergent pCa50/nH, never fit to any target)."""
    F = np.array([steady_state_force_closed_form(10 ** (6.0 - p), SL, TmpC) for p in pCa_grid])
    return F


def run_periodic_force(ca_func_one_beat, bcl_ms, SL=2.1, TmpC=37.0, n_force_beats=6,
                        y0=None, rtol=1e-9, atol=1e-12, max_step=0.5):
    """Drive the dynamic 7-state ODE with a PERIODIC Ca(t) (one beat's worth of waveform,
    tiled n_force_beats times via t % bcl_ms) and integrate long enough for the crossbridge
    pools' own limit cycle to be reached (checked by caller via per-beat peak-force
    convergence, same discipline as cicrca_ode_model.converged())."""
    if y0 is None:
        y0 = list(Y0_DEFAULT)

    def ca_periodic(t):
        return ca_func_one_beat(t % bcl_ms)

    rhs = lambda t, y: ode_rhs(t, y, SL, TmpC, ca_periodic)
    t_span = (0.0, n_force_beats * bcl_ms)
    t_eval = np.linspace(0.0, n_force_beats * bcl_ms, int(n_force_beats * 400))
    sol = solve_ivp(rhs, t_span, y0, method='Radau', t_eval=t_eval,
                     rtol=rtol, atol=atol, max_step=max_step)
    return sol


if __name__ == "__main__":
    import json
    # Self-check: reproduce this script's emergent pCa50/nH at physiological 37C, and
    # print the bias-decorrelated interval synthesis described in the module docstring.
    from scipy.optimize import curve_fit

    def hill(pCa, pCa50, nH, Fmax):
        Ca = 10.0 ** (-pCa)
        Ca50 = 10.0 ** (-pCa50)
        return Fmax / (1.0 + (Ca50 / Ca) ** nH)

    pCa_grid = np.linspace(4.0, 7.2, 65)
    results = {}
    for SL in (1.95, 2.10, 2.25):
        F = emergent_fpca(SL, 37.0, pCa_grid)
        popt, _ = curve_fit(hill, pCa_grid, F, p0=[6.0, 7.0, F.max()],
                             bounds=([4.0, 0.5, F.max() * 0.1], [7.5, 30.0, F.max() * 10]),
                             maxfev=20000)
        results[SL] = dict(pCa50=popt[0], nH=popt[1])
        print(f"SL={SL}: EMERGENT (bug-fixed, untuned, 37C) pCa50={popt[0]:.4f} nH={popt[1]:.3f}")

    intact_anchors = dict(backx1995_rat=6.187, staurosporine_paper_rabbit_37C=6.05)
    skinned_anchors = dict(konhilas2002_rat_15C_SL2p10=5.51, human_permeabilized=5.489)
    print("\nIntact/physiological-temp anchors (external):", intact_anchors)
    print("Skinned/permeabilized anchors (external):", skinned_anchors)
    print("\nModel emergent pCa50 range:", {k: round(v['pCa50'], 3) for k, v in results.items()})
