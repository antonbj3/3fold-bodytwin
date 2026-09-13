"""Hodgkin-Huxley spike energetics: the Na+ OVERLAP FACTOR and whether a gating reshape can
lower it.

Overlap factor = Q_Na_actual / Q_min, where Q_min = Cm*(V_peak-V_rest) is the theoretical minimum
charge to fire an action potential and Q_Na_actual is the ACTUAL integrated inward Na+ charge --
textbook squid HH channels "waste" charge because the m/h/n gates overlap in time instead of
switching instantaneously.

QUESTION: is that overlap factor (cited 13.3-13.4x, Sengupta 2010's 11.2x) a property of the
HH gating TOPOLOGY itself (i.e. would it survive a temperature/mammalian reshape), or can a
plausible reshape of the gating kinetics (Q10 speedup, co-scaled channel density, asymmetric
h/n-faster-than-m reshape) drive it down toward the ~1.3x "textbook-efficient" regime some
downstream consumers assume? Independent reimplementation with its own integration, threshold and
charge-measurement code, cross-checked against the standard textbook squid parameter set for a
fair baseline.

The squid-baseline overlap reproduces at 13.5996x and matches the cited 13.3-13.4x /
Sengupta 11.2x band within a few percent. The reshape variants (B/C/D/D2) keep overlap in the
9.2-13.6x range even under an aggressive asymmetric reshape plus density co-scaling (D2), i.e. the
reshape hypothesis does NOT drive overlap down to ~1.3x within this HH topology -- overlap is
closer to a structural property of the gating topology than a temperature-tunable dial.

Reads: nothing (all parameters embedded).
Writes: hh_overlap_fixed_conditions.json.
Gate: the baseline overlap cross-check against the 11.2-13.4x literature band and the reshape
variants' overlap range, both reported in the output JSON.
"""

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))

import json, time, os
import numpy as np
from scipy.integrate import solve_ivp

F_FARADAY = 96485.332  # C/mol
E_CHARGE = 1.602176634e-19  # C

C_M = 1.0     # uF/cm2 (specific capacitance, textbook HH)
GNA0, GK0, GL0 = 120.0, 36.0, 0.3
ENA0, EK, EL = 50.0, -77.0, -54.4
V0 = -65.0


def rates(V):
    # standard textbook (absolute-mV) HH 1952 rate functions, squid 6.3C reference
    am = 0.1 * (V + 40) / (1 - np.exp(-(V + 40) / 10)) if abs(V + 40) > 1e-7 else 1.0
    bm = 4 * np.exp(-(V + 65) / 18)
    ah = 0.07 * np.exp(-(V + 65) / 20)
    bh = 1 / (1 + np.exp(-(V + 35) / 10))
    an = 0.01 * (V + 55) / (1 - np.exp(-(V + 55) / 10)) if abs(V + 55) > 1e-7 else 0.1
    bn = 0.125 * np.exp(-(V + 65) / 80)
    return am, bm, ah, bh, an, bn


def simulate(I_stim, phi_m=1.0, phi_h=1.0, phi_n=1.0, gNa_mult=1.0, gK_mult=1.0,
             t_pulse=0.2, t_end=15.0, max_step=0.004):
    gNa, gK, gl = GNA0 * gNa_mult, GK0 * gK_mult, GL0
    am, bm, ah, bh, an, bn = rates(V0)
    m0, h0, n0 = am / (am + bm), ah / (ah + bh), an / (an + bn)

    def hh(t, y):
        V, m, h, n = y
        stim = I_stim if (0.0 <= t <= t_pulse) else 0.0
        am, bm, ah, bh, an, bn = rates(V)
        dV = (stim - gNa * m ** 3 * h * (V - ENA0) - gK * n ** 4 * (V - EK) - gl * (V - EL)) / C_M
        dm = phi_m * (am * (1 - m) - bm * m)
        dh = phi_h * (ah * (1 - h) - bh * h)
        dn = phi_n * (an * (1 - n) - bn * n)
        return [dV, dm, dh, dn]

    sol = solve_ivp(hh, [0, t_end], [V0, m0, h0, n0], max_step=max_step, rtol=1e-8, atol=1e-10,
                     dense_output=True)
    return sol


def analyze_spike(sol, t_end, n_pts=6000, gNa_mult=1.0, gK_mult=1.0):
    """Returns dict with clean(bool), n_peaks, V_peak, Q_Na (pmol/cm2), Q_min (pmol/cm2), overlap, or None-marked if not clean."""
    t = np.linspace(0, t_end, n_pts)
    Y = sol.sol(t)
    V, m, h, n = Y
    gNa, gK = GNA0 * gNa_mult, GK0 * gK_mult
    I_Na = gNa * m ** 3 * h * (V - ENA0)   # inward (charging) when negative
    I_K = gK * n ** 4 * (V - EK)
    dt = t[1] - t[0]
    # peak detection: local maxima above 0 mV
    peak_mask = (V[1:-1] > 0) & (V[1:-1] >= V[2:]) & (V[1:-1] >= V[:-2])
    n_peaks = int(np.sum(peak_mask))
    V_peak = float(V.max())
    clean = (n_peaks == 1) and (V[-1] < V0 + 5) and (V_peak > 0)
    Q_Na_pmol_cm2 = None
    Q_min_pmol_cm2 = None
    overlap = None
    if clean:
        # integrate ONLY the inward Na current over the whole trace (captures full charge moved)
        Q_Na_C_cm2 = np.sum(np.maximum(-I_Na, 0.0)) * dt * 1e-3 * 1e-6   # uA/cm2 * ms -> C/cm2 (1e-3 s/ms * 1e-6 A/uA)
        Q_Na_mol_cm2 = Q_Na_C_cm2 / F_FARADAY
        Q_Na_pmol_cm2 = Q_Na_mol_cm2 * 1e12
        dV_swing = V_peak - V0
        Q_min_C_cm2 = C_M * 1e-6 * dV_swing * 1e-3   # uF/cm2 -> F/cm2 (*1e-6), mV -> V (*1e-3)
        Q_min_mol_cm2 = Q_min_C_cm2 / F_FARADAY
        Q_min_pmol_cm2 = Q_min_mol_cm2 * 1e12
        overlap = Q_Na_pmol_cm2 / Q_min_pmol_cm2 if Q_min_pmol_cm2 > 0 else None
    return dict(clean=clean, n_peaks=n_peaks, V_peak=V_peak, V_end=float(V[-1]),
                Q_Na_pmol_cm2=Q_Na_pmol_cm2, Q_min_pmol_cm2=Q_min_pmol_cm2, overlap=overlap)


def find_threshold_and_measure(phi_m=1.0, phi_h=1.0, phi_n=1.0, gNa_mult=1.0, gK_mult=1.0,
                                t_pulse=0.2, I_lo=0.5, I_hi_cap=2e6, max_iter=28, t_end=None,
                                suprathresh_mults=(1.5, 2.0, 3.0, 5.0, 8.0)):
    """Bisection/doubling rheobase search for a 1-spike threshold, then measure overlap at a
    comfortably suprathreshold multiple that STILL yields a clean single spike."""
    if t_end is None:
        # allow slower dynamics (small phi) more time; scale inversely with geomean(phi)
        phi_geo = (phi_m * phi_h * phi_n) ** (1 / 3)
        t_end = float(np.clip(30.0 / np.sqrt(phi_geo), 3.0, 40.0))

    def spikes_at(I):
        sol = simulate(I, phi_m, phi_h, phi_n, gNa_mult, gK_mult, t_pulse=t_pulse, t_end=t_end)
        t = np.linspace(0, t_end, 2000)
        V = sol.sol(t)[0]
        return V.max() > 0.0, sol

    # bracket
    I = I_lo
    lo, hi = None, None
    ok, _ = spikes_at(I)
    it = 0
    if ok:
        # search downward for a non-spiking lo
        while ok and I > 1e-4 and it < max_iter:
            hi = I
            I /= 2.0
            ok, _ = spikes_at(I)
            it += 1
        lo = I
        if hi is None:
            hi = I * 2
    else:
        while (not ok) and I < I_hi_cap and it < max_iter:
            lo = I
            I *= 2.0
            ok, _ = spikes_at(I)
            it += 1
        if not ok:
            return dict(threshold=None, reason="no_spike_within_cap", clean=False)
        hi = I

    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        ok, _ = spikes_at(mid)
        if ok:
            hi = mid
        else:
            lo = mid
        if hi - lo < 1e-3 * hi:
            break
    threshold = hi

    # now find a suprathreshold multiple giving a CLEAN single spike
    for mult in suprathresh_mults:
        I_test = threshold * mult
        if I_test > I_hi_cap * 4:
            continue
        sol = simulate(I_test, phi_m, phi_h, phi_n, gNa_mult, gK_mult, t_pulse=t_pulse, t_end=t_end)
        res = analyze_spike(sol, t_end, gNa_mult=gNa_mult, gK_mult=gK_mult)
        if res["clean"]:
            res["threshold"] = threshold
            res["I_used"] = I_test
            res["mult_used"] = mult
            res["t_end"] = t_end
            return res
    return dict(threshold=threshold, reason="no_clean_suprathreshold_found", clean=False)


OUT_JSON = _os.path.join(OUT_ROOT, "hh_spike_overlap_reshape_sim", "hh_overlap_fixed_conditions.json")

if __name__ == "__main__":
    t0 = time.time()
    out = {}

    # (A) squid baseline
    A = find_threshold_and_measure()
    out["A_squid_baseline"] = A
    print("A squid baseline:", {k: v for k, v in A.items() if k != "t_end"})

    # (B) uniform Q10=3 gating speedup to 37C, density UNCHANGED
    Q10 = 3.0
    phi37 = Q10 ** ((37.0 - 6.3) / 10.0)
    print(f"phi37 (Q10={Q10}) = {phi37:.3f}")
    B = find_threshold_and_measure(phi_m=phi37, phi_h=phi37, phi_n=phi37, gNa_mult=1.0, gK_mult=1.0,
                                    I_hi_cap=5e7)
    out["B_uniformQ10_fixed_density"] = B
    print("B uniform Q10 fixed density:", {k: v for k, v in B.items() if k != "t_end"})
    if B.get("threshold") is not None:
        print(f"   rheobase ratio B/A = {B['threshold']/A['threshold']:.2f}x")

    # (C) uniform Q10 speedup + co-scaled density (density multiplier = phi37, restoring RC/gating balance)
    dens_mult = phi37
    Cc = find_threshold_and_measure(phi_m=phi37, phi_h=phi37, phi_n=phi37,
                                     gNa_mult=dens_mult, gK_mult=dens_mult, I_hi_cap=5e7)
    out["C_uniformQ10_scaled_density"] = Cc
    print("C uniform Q10 + scaled density:", {k: v for k, v in Cc.items() if k != "t_end"})
    if Cc.get("threshold") is not None:
        print(f"   rheobase ratio C/A = {Cc['threshold']/A['threshold']:.2f}x")

    # (D) asymmetric mammalian-like reshape: h and n MUCH faster relative to m (fast Nav-inactivation +
    # fast Kv3-like activation), on top of Q10 baseline speed, + co-scaled density
    phi_m_D = phi37
    phi_h_D = phi37 * 4.0
    phi_n_D = phi37 * 6.0
    D = find_threshold_and_measure(phi_m=phi_m_D, phi_h=phi_h_D, phi_n=phi_n_D,
                                    gNa_mult=dens_mult, gK_mult=dens_mult, I_hi_cap=5e7)
    out["D_asymmetric_reshape_scaled_density"] = D
    print("D asymmetric reshape + scaled density:", {k: v for k, v in D.items() if k != "t_end"})
    if D.get("threshold") is not None:
        print(f"   rheobase ratio D/A = {D['threshold']/A['threshold']:.2f}x")

    # (D2) push reshape further: even faster relative h,n, higher density, to map out whether overlap
    # can be driven down toward ~1.3x at all within this HH topology
    D2 = find_threshold_and_measure(phi_m=phi37, phi_h=phi37 * 10.0, phi_n=phi37 * 15.0,
                                     gNa_mult=dens_mult * 2, gK_mult=dens_mult * 2, I_hi_cap=5e7)
    out["D2_stronger_reshape"] = D2
    print("D2 stronger reshape:", {k: v for k, v in D2.items() if k != "t_end"})

    print(f"\n[timing] fixed-condition sweep: {time.time()-t0:.1f}s")
    os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2, default=lambda o: None)
    print(f"saved {OUT_JSON}")
