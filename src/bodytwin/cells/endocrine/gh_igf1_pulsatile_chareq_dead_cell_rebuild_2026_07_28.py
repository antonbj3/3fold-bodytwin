#!/usr/bin/env python3
"""GH/IGF-1 PULSATILE AXIS: characteristic-equation check of a recorded two-loop delay model.

Recomputes, from the delayed two-pole negative-feedback characteristic equation
(lambda+a1)(lambda+a2) + K*exp(-lambda*tau) = 0 evaluated at the marginal-stability point
lambda = i*w, the delays tau1 (fast R-A loop, 36 min target period) and tau2 (slow H-P loop,
150 min target period) that a recorded parameter set claims, verifies the complex residual of
the characteristic equation at those delays, and re-derives the tau2/tau1 mode-locking ratios of
the recorded multiplier sweep against the recorded bifurcation bracket.

Reads: nothing (parameters embedded). Writes: nothing (prints its numbers).
Finding: the recorded "~25% above the nearest mode-locked ratio" figure is not reproducible from
the same two numbers; both independent routes give ~33%.
"""
import numpy as np
from scipy.optimize import brentq

# ---- parameters of the recorded two-loop model ----
aR, aA, aP, aH = 0.09902, 0.06931, 0.01733, 0.04077
tau1_target_period = 36.0   # min, design target (fast loop, R-A)
tau2_target_period = 150.0  # min, design target (slow loop, H-P)
tau1_claimed = 5.123
tau2_claimed = 27.792

# characteristic eq for a delayed 2-pole negative feedback loop:
# (lambda+a1)(lambda+a2) + K*exp(-lambda*tau) = 0
# at lambda=i*w (Hopf/marginal stability):
#   K*cos(w*tau) = a1*a2 + w^2         (matching real part; sign convention verified below)
#   K*sin(w*tau) = w*(a1+a2)           (imaginary part)
# => tan(w*tau) = w*(a1+a2) / (a1*a2 + w^2)   [cell states w^2-a1a2 in the denominator; test both signs]

def solve_tau_for_period(a1, a2, period, denom_sign=-1):
    w = 2 * np.pi / period
    denom = w**2 + denom_sign * a1 * a2   # denom_sign=-1 -> w^2-a1a2 (cell's stated form); +1 -> w^2+a1a2
    rhs = w * (a1 + a2) / denom
    theta = np.arctan(rhs)  # principal branch in (-pi/2,pi/2)
    # tau must be positive; scan branches theta+n*pi for smallest positive tau
    for n in range(0, 5):
        tau_try = (theta + n * np.pi) / w
        if tau_try > 0:
            return tau_try, w, denom
    return None, w, denom

print("=== FAST LOOP (R-A), target period=36min ===")
for sign in (-1, 1):
    tau1, w1, denom1 = solve_tau_for_period(aR, aA, tau1_target_period, denom_sign=sign)
    print(f"  denom_sign={sign:+d} (w^2{'+' if sign>0 else '-'}a1a2): tau1={tau1:.4f}min  w={w1:.6f}")
print(f"  claimed tau1 = {tau1_claimed}min")

print("\n=== SLOW LOOP (H-P), target period=150min ===")
for sign in (-1, 1):
    tau2, w2, denom2 = solve_tau_for_period(aH, aP, tau2_target_period, denom_sign=sign)
    print(f"  denom_sign={sign:+d} (w^2{'+' if sign>0 else '-'}a1a2): tau2={tau2:.4f}min  w={w2:.6f}")
print(f"  claimed tau2 = {tau2_claimed}min")

# pick the branch matching the cell's stated form (denom = w^2-a1a2) and report exact match
tau1_fwd, w1, _ = solve_tau_for_period(aR, aA, tau1_target_period, denom_sign=-1)
tau2_fwd, w2, _ = solve_tau_for_period(aH, aP, tau2_target_period, denom_sign=-1)
print(f"\n[MATCH CHECK, cell's stated tan(w*tau)=w(a1+a2)/(w^2-a1a2) form]")
print(f"  tau1: computed={tau1_fwd:.4f} claimed={tau1_claimed} rel_err={abs(tau1_fwd-tau1_claimed)/tau1_claimed:.4%}")
print(f"  tau2: computed={tau2_fwd:.4f} claimed={tau2_claimed} rel_err={abs(tau2_fwd-tau2_claimed)/tau2_claimed:.4%}")

# ---- verify the char-eq residual directly (complex arithmetic, lambda=i*w) ----
def char_eq_residual(w, a1, a2, tau, K):
    lam = 1j * w
    return (lam + a1) * (lam + a2) + K * np.exp(-lam * tau)

def K_from_w_tau(a1, a2, w, tau):
    # from K*cos(w tau)=w^2-a1a2, K*sin(w tau)=w(a1+a2) (derived from (iw+a1)(iw+a2)+K*e^-iwtau=0)
    return np.hypot(w**2 - a1 * a2, w * (a1 + a2))

K1 = K_from_w_tau(aR, aA, w1, tau1_fwd)
K2 = K_from_w_tau(aH, aP, w2, tau2_fwd)
res1 = char_eq_residual(w1, aR, aA, tau1_fwd, K1)
res2 = char_eq_residual(w2, aH, aP, tau2_fwd, K2)
print(f"\n[CHAR EQ RESIDUAL] fast loop |residual|={abs(res1):.3e} (K1={K1:.5f})")
print(f"[CHAR EQ RESIDUAL] slow loop |residual|={abs(res2):.3e} (K2={K2:.5f})")

# ---- mult sweep: tau2 *= {0.5,0.75,1.0,1.25,1.5}, using claimed tau1=5.123, tau2_baseline=27.792 ----
print("\n=== MULT SWEEP tau2/tau1 ratios (baseline tau1=5.123, tau2=27.792) ===")
tau1_b, tau2_b = tau1_claimed, tau2_claimed
mults = [0.5, 0.75, 1.0, 1.25, 1.5]
ratios = []
for mlt in mults:
    tau2_m = tau2_b * mlt
    ratio = tau2_m / tau1_b
    ratios.append(ratio)
    zone = "MODE-LOCKED (claimed)" if mlt in (0.5, 0.75) else "clean/monotonic (claimed)"
    print(f"  mult={mlt}: tau2={tau2_m:.4f}min  tau2/tau1={ratio:.4f}  -> {zone}")

lo_bracket, hi_bracket = 4.07, 5.43
print(f"\n[BRACKET CONSISTENCY CHECK] claimed bifurcation bracket tau2/tau1 in ({lo_bracket},{hi_bracket})")
print(f"  mult=0.75 ratio = {ratios[1]:.4f}  vs bracket lower edge {lo_bracket}  (match to 3dp: {abs(ratios[1]-lo_bracket)<0.005})")
print(f"  mult=1.0  ratio = {ratios[2]:.4f}  vs bracket upper edge {hi_bracket}  (match to 3dp: {abs(ratios[2]-hi_bracket)<0.005})")
print(f"  mult=1.0 ratio {ratios[2]:.4f} >= 5.4 threshold (honest_gaps 'monotonic only for tau2/tau1>=~5.4')? {ratios[2] >= 5.4}")
margin_over_54 = (ratios[2] - 5.4) / 5.4
print(f"  margin over 5.4 threshold: {margin_over_54:.4%} (thin, as honest_gaps states '~25% above nearest locked point' -- check that separately)")
nearest_locked = ratios[1]  # mult=0.75, ratio=4.070, claimed locked
pct_above_locked = (ratios[2] - nearest_locked) / nearest_locked
print(f"  mult=1.0 ratio is {pct_above_locked:.1%} above nearest-confirmed-locked ratio (mult=0.75); honest_gaps claims '~25% above'")

print("\n[SELF-CONSISTENCY VERDICT]")
print(f"  Locked set claimed = {{mult 0.5({ratios[0]:.3f}), mult 0.75({ratios[1]:.3f})}}")
print(f"  Bracket claimed = ({lo_bracket},{hi_bracket})")
print(f"  mult=0.75's ratio ({ratios[1]:.4f}) sits AT/BELOW the bracket's lower edge ({lo_bracket}) -- consistent (locked point below/at threshold).")
print(f"  mult=1.0's ratio ({ratios[2]:.4f}) sits AT/ABOVE the bracket's upper edge ({hi_bracket}) -- consistent (unlocked point above threshold).")
print(f"  => the two claims (discrete locked-set and bracket range) are MUTUALLY CONSISTENT, not contradictory.")

print("\n[SELF-CONTRADICTION CHECK: honest_gaps' own '~25% above nearest confirmed-locked point']")
pct_vs_nearest_locked = (ratios[2] - ratios[1]) / ratios[1]
pct_vs_bracket_edges = (hi_bracket - lo_bracket) / lo_bracket
print(f"  (baseline ratio - nearest-locked ratio)/nearest-locked = ({ratios[2]:.4f}-{ratios[1]:.4f})/{ratios[1]:.4f} = {pct_vs_nearest_locked:.1%}")
print(f"  (bracket hi-lo)/lo = ({hi_bracket}-{lo_bracket})/{lo_bracket} = {pct_vs_bracket_edges:.1%}")
print(f"  honest_gaps' own stated figure = '~25%' -- both independent recomputations from the SAME two numbers give ~33%, not ~25%.")
print(f"  VERDICT: DIFFERS -- the cell's honest_gaps '~25% above nearest confirmed-locked point' does not match its own bracket/mult-sweep numbers (33.3-33.4% by both routes). Minor, but a genuine self-inconsistency, not merely imprecise rounding (33 vs 25 is an 8-point gap).")
