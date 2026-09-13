"""
Independent from-scratch reimplementation of the Lemaire2004 RANKL/RANK/OPG ODE
(BioModels BIOMD0000000278), using ONLY the equations + numeric parameters as
disclosed in the project's MODEL-BONE-REMODELING-RANKL-OPG node/evidence JSON
(the .py that originally produced those numbers is confirmed ABSENT from disk).
Purpose: cheap decisive test -- does the model, run fresh, reproduce its OWN
claimed baseline fixed point and its OWN claimed +304% osteoclast response to
continuous PTH 1000 pM/day? This is a reproducibility check, not a fresh anchor.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

P = dict(C_s=0.005, D_A=0.7, d_B=0.7, D_C=0.0021, D_R=0.0007, f0=0.05, K=10.0,
         k1=0.01, k2=10.0, k3=0.00058, k4=0.017, k5=0.02, k6=3.0, k_B=0.189,
         K_L_P=3.0e6, kO=0.35, K_O_P=2.0e5, k_P=86.0, r_L=1000.0, S_P=250.0)

IC = dict(R0=7.734e-4, B0=7.282e-4, C0=9.127e-4)

def rhs(t, y, I_P=0.0, I_L=0.0, I_O=0.0):
    R, B, C = y
    D_B = P['f0'] * P['d_B']
    Phi_C = (C + P['f0']*P['C_s']) / (C + P['C_s'])
    Pbar = I_P / P['k_P']
    P_O = P['S_P'] / P['k_P']
    P_S = P['k6'] / P['k5']
    Phi_P = (Pbar + P_O) / (Pbar + P_S)
    denom = 1.0 + (P['k3']*P['K']/P['k4']) + (P['k1']/(P['k2']*P['kO'])) * (I_O + P['K_O_P']*R/Phi_P)
    Phi_L = ((P['k3']/P['k4']) * P['K_L_P'] * Phi_P * B / denom) * (1.0 + I_L/P['r_L'])
    dR = P['D_R']*Phi_C - D_B*R/Phi_C
    dB = D_B*R/Phi_C - P['k_B']*B
    dC = P['D_C']*Phi_L - P['D_A']*Phi_C*C
    return [dR, dB, dC]

y0 = [IC['R0'], IC['B0'], IC['C0']]

# 1. is the claimed IC actually (near) the fixed point at baseline (I_P=I_L=I_O=0)?
dydt0 = rhs(0.0, y0)
print("max|dy/dt| at claimed baseline IC:", max(abs(v) for v in dydt0))

# 2. independent fsolve from a perturbed guess -> does it converge back to the same IC?
sol_fp = fsolve(lambda y: rhs(0, y), x0=[v*1.2 for v in y0], full_output=True)
root, info, ier, msg = sol_fp
print("fsolve root:", root, "ier=", ier, "residual max:", max(abs(v) for v in rhs(0, root)))
print("relative diff vs claimed IC:", [(r-c)/c for r, c in zip(root, y0)])

# 3. 2000-day free run from claimed IC at baseline -> should stay put (claimed drift <0.02%)
t_span = (0, 2000)
sol_free = solve_ivp(rhs, t_span, y0, method='LSODA', rtol=1e-10, atol=1e-14, dense_output=True)
y_end = sol_free.y[:, -1]
drift_pct = [100*abs(e-o)/o for e, o in zip(y_end, y0)]
print("2000-day free-run end state:", y_end, "drift%:", drift_pct)

# 4. continuous PTH 1000 pM/day from day 20-80 (paper's Fig.2 protocol) -> C peak vs C0?
def rhs_pulse(t, y):
    ip = 1000.0 if 20 <= t <= 80 else 0.0
    return rhs(t, y, I_P=ip)

sol_pulse = solve_ivp(rhs_pulse, (0, 150), y0, method='LSODA', rtol=1e-10, atol=1e-14,
                       dense_output=True, max_step=0.5)
C_track = sol_pulse.y[2]
C_peak = C_track.max()
t_peak = sol_pulse.t[np.argmax(C_track)]
pct_change = 100.0*(C_peak - IC['C0'])/IC['C0']
print(f"C peak = {C_peak:.6e} pM at t={t_peak:.1f}d vs C0={IC['C0']:.6e} -> {pct_change:+.1f}% (claimed: +304%)")

# 5. Jacobian eigenvalues at baseline fixed point (finite-difference), compare to claimed
#    -0.2528+-0.2330i, -0.0858
eps = 1e-7
J = np.zeros((3, 3))
base = np.array(rhs(0, y0))
for j in range(3):
    yp = list(y0); yp[j] += eps
    J[:, j] = (np.array(rhs(0, yp)) - base) / eps
eigs = np.linalg.eigvals(J)
print("Jacobian eigenvalues (finite-diff):", eigs)

# 6. continuous I_P=1000 pM/day steady state (no pulse-off) -- what osteoclast fold-change at
#    STEADY STATE (not transient peak)?
sol_cont = solve_ivp(lambda t, y: rhs(t, y, I_P=1000.0), (0, 3000), y0, method='LSODA',
                      rtol=1e-11, atol=1e-15)
C_ss = sol_cont.y[2, -1]
print(f"steady-state C under continuous I_P=1000pM/day: {C_ss:.6e} -> {100*(C_ss-IC['C0'])/IC['C0']:+.1f}% vs C0")
