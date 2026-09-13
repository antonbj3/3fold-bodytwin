"""Urea-cycle flux model: recalibration of the nitrogen input Jin.

Independent from-scratch reimplementation of an 8-state urea-cycle ODE (NH3, CP, ORN, CIT, ASA,
ARG + UREA/OROTATE accumulators); every equation and numeric parameter used is embedded below.

Purpose: run the named falsifier -- recalibrate Jin to 500 mmol/day and check whether the 1-3%
OTC-crisis bifurcation band shifts. The baseline reproduction is the fidelity gate; the
Jin=500 mmol/day results are the new finding.

Results printed:
  1. baseline (Jin=272.8mmol/day=0.19mM/min) reproduces the cited NH3_ss=284.8umol/L
     exactly, from three different initial conditions, with the exact conservation law
     ORN+CIT+ASA+ARG=10.0mM holding to <1e-6.
  2. naive Jin recalibration to 500mmol/day (holding every other parameter, incl. CPS1
     capacity Vmax1, fixed) BREAKS the model: NH3_ss=652.2umol/L even at 100% (fully healthy)
     OTC activity -- already crisis-range hyperammonemia with NO unhealthy enzyme assumed.
     The model has NO OTC-activity level at which NH3 < 500umol/L under this naive test.
  3. forced adversary to its strongest form (OODA, not a one-shot negative): if CPS1 capacity
     (Vmax1) is assumed to co-scale with nitrogen intake (Jin_new/Jin_old = 1.827x, the same
     ratio as the RDA-floor-vs-typical-protein-intake regime choice), a healthy baseline is
     restored (NH3_ss=287.6umol/L, matching the original 284.8 to 1%) and the OTC-crisis
     residual-activity band shifts from [2.21%,3.32%] (old Jin) to [4.48%,6.44%] (new Jin +
     co-scaled capacity) -- roughly DOUBLING the % residual activity needed to cross the same
     absolute NH3 crisis thresholds.

Reads: nothing. Writes: nothing (prints only).
Gate: baseline NH3_ss must reproduce 284.8 umol/L from three initial conditions with the
ORN+CIT+ASA+ARG=10.0mM conservation law holding.
"""
import numpy as np
from scipy.integrate import solve_ivp

P = dict(
    Vmax1=1.0, Km1=1.2, Ki1=0.1,          # CPS1 (NH3 -> carbamoyl-phosphate), product-inhibited by CP
    Vmax2=10.0, Km2c=0.05, Km2o=2.0,      # OTC (CP + ornithine -> citrulline)
    Vmax3=1.5, Km3=0.15,                  # ASS1 (citrulline -> argininosuccinate)
    Vmax4=8.0, Km4=0.05,                  # ASL (argininosuccinate -> arginine)
    Vmax5=15.0, Km5=5.0,                  # ARG1 (arginine -> urea + ornithine, closes cycle)
    Vmax_oro=0.05, Km_oro=0.1,            # CP -> orotate spillover (pyrimidine-pathway leak)
    ORN_pool=10.0,                        # conserved ORN+CIT+ASA+ARG (mM)
)

JIN_BASELINE_MM_MIN = 0.19               # calibrated value (272.8 mmol/day)
JIN_TASK_ANCHOR_MM_MIN = 500.0 / 1440.0  # Rolfe-Brown-consistent task anchor (500 mmol/day)


def rhs(t, y, Jin, Vmax2_scale=1.0, Vmax1_scale=1.0):
    NH3, CP, ORN, CIT, ASA, ARG = y
    J1 = (P["Vmax1"] * Vmax1_scale) * NH3 / (P["Km1"] * (1 + CP / P["Ki1"]) + NH3)
    Vmax2_eff = P["Vmax2"] * Vmax2_scale
    J2 = Vmax2_eff * (CP / (P["Km2c"] + CP)) * (ORN / (P["Km2o"] + ORN))
    J3 = P["Vmax3"] * CIT / (P["Km3"] + CIT)
    J4 = P["Vmax4"] * ASA / (P["Km4"] + ASA)
    J5 = P["Vmax5"] * ARG / (P["Km5"] + ARG)
    Joro = P["Vmax_oro"] * CP / (P["Km_oro"] + CP)
    return [Jin - J1, J1 - J2 - Joro, J5 - J2, J2 - J3, J3 - J4, J4 - J5]


def steady_state(Jin, Vmax2_scale=1.0, Vmax1_scale=1.0, T=200000.0, y0=None):
    if y0 is None:
        y0 = [0.3, 0.05, 2.5, 2.5, 2.5, 2.5]
    sol = solve_ivp(rhs, [0, T], y0, args=(Jin, Vmax2_scale, Vmax1_scale),
                     method="LSODA", rtol=1e-10, atol=1e-13)
    yf = sol.y[:, -1]
    max_resid = np.max(np.abs(rhs(T, yf, Jin, Vmax2_scale, Vmax1_scale)))
    return yf, max_resid


def find_otc_crossing(Jin, target_mM, Vmax1_scale=1.0, lo=1e-5, hi=1.0, y0=None):
    """Bisection on residual OTC-activity fraction (Vmax2_scale) for NH3_ss == target_mM."""
    def f(scale):
        yf, _ = steady_state(Jin, scale, Vmax1_scale, y0=y0)
        return yf[0] - target_mM

    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        return None
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if fm * flo > 0:
            lo, flo = mid, fm
        else:
            hi = mid
    return 0.5 * (lo + hi)


if __name__ == "__main__":
    y0 = [0.28, 0.05, 2.5, 2.5, 2.5, 2.5]

    print("=== FIDELITY CHECK: baseline Jin=0.19mM/min (272.8mmol/day) ===")
    for ic_label, y0test in [("mid", y0), ("low", [0.1, 0.02, 1, 1, 1, 7]), ("high", [0.5, 0.1, 7, 1, 1, 1])]:
        yf, resid = steady_state(JIN_BASELINE_MM_MIN, 1.0, 1.0, y0=y0test)
        print(f"  IC={ic_label}: NH3_ss={yf[0]*1000:.2f}umol/L (cited: 284.8) "
              f"max|dydt|={resid:.2e} conserved={sum(yf[2:]):.6f} (target 10.0)")

    print(f"\n=== NAIVE FALSIFIER: Jin -> 500mmol/day ({JIN_TASK_ANCHOR_MM_MIN:.5f}mM/min), all else fixed ===")
    yf_naive, _ = steady_state(JIN_TASK_ANCHOR_MM_MIN, 1.0, 1.0, y0=y0)
    print(f"  NH3_ss @ 100% OTC activity = {yf_naive[0]*1000:.1f}umol/L "
          f"(vs Haberle healthy <80umol/L, vs 500umol/L crisis threshold)")
    c500 = find_otc_crossing(JIN_TASK_ANCHOR_MM_MIN, 0.5, 1.0, y0=y0)
    print(f"  500umol/L crossing exists at any OTC activity in [0,100%]? "
          f"{'NO -- model predicts crisis-range NH3 at every activity level' if c500 is None else f'yes, at {c500*100:.2f}%'}")

    scale_factor = JIN_TASK_ANCHOR_MM_MIN / JIN_BASELINE_MM_MIN
    print(f"\n=== FORCED-ADVERSARY RESCUE: Vmax1 (CPS1 capacity) co-scaled {scale_factor:.3f}x with Jin ===")
    yf_rescue, _ = steady_state(JIN_TASK_ANCHOR_MM_MIN, 1.0, scale_factor, y0=y0)
    print(f"  NH3_ss @ 100% OTC activity = {yf_rescue[0]*1000:.1f}umol/L (baseline 284.8umol/L)")
    c500_old = find_otc_crossing(JIN_BASELINE_MM_MIN, 0.5, 1.0, y0=y0)
    c1000_old = find_otc_crossing(JIN_BASELINE_MM_MIN, 1.0, 1.0, y0=y0)
    c500_new = find_otc_crossing(JIN_TASK_ANCHOR_MM_MIN, 0.5, scale_factor, y0=y0)
    c1000_new = find_otc_crossing(JIN_TASK_ANCHOR_MM_MIN, 1.0, scale_factor, y0=y0)
    print(f"  OLD OTC-crisis band (Jin=272.8mmol/day): [{c1000_old*100:.2f}%, {c500_old*100:.2f}%]")
    print(f"  NEW OTC-crisis band (Jin=500mmol/day, capacity co-scaled): [{c1000_new*100:.2f}%, {c500_new*100:.2f}%]")
