"""Ultradian insulin oscillation: Hopf-margin test of the Sturis 1991 / BIOMD0000000382 model.

QUESTION: using the verbatim Sturis 1991 / BIOMD0000000382 6-state glucose-insulin-hepatic-delay
model (published parameters, NOT refit), is the equilibrium (the constant-infusion steady state
x*=y*=z* independent of the delay t3) linearly UNSTABLE at the model's physiological
hepatic-delay t3=36min -- i.e. does the model structurally support a sustained ultradian
oscillation (a Hopf-type instability), not just an assumed one? Cross-checked TWO ways: (1) an
analytically-derived 6x6 Jacobian at equilibrium vs (2) a finite-difference Jacobian of the full
nonlinear RHS -- if they disagree, the analytic derivation has a bug.

Gate: the analytic and finite-difference Jacobians must agree, and the dominant eigenvalue pair at
t3=36min must have a positive real part with an implied period inside the physiologically reported
ultradian insulin-oscillation range (~90-140 min).

Reads: nothing. Writes: insulin_equilibrium.json (equilibrium, Jacobian agreement, eigenvalues).
"""
import json, os
import numpy as np
from scipy.optimize import fsolve

# ===== Sturis 1991 / BIOMD0000000382 verbatim published parameters =====
v1, v2, v3 = 3.0, 11.0, 10.0
t1, t2 = 6.0, 100.0
I_infusion = 216.0
E = 0.21
T3_BASELINE = 36.0

def f1(z):  # glucose-stimulated secretion
    return 209.0/(1.0+np.exp(-z/(300*v3)+6.6))
def f2(z):  # insulin-independent utilization
    return 72.0*(1.0-np.exp(-z/(144*v3)))
def f3(z):
    return 0.01*z/v3
def f4(y):  # insulin-dependent utilization factor
    arg = y*(1.0/v2 + 1.0/(E*t2))
    return 90.0/(1.0+np.exp(-1.772*np.log(arg)+7.76))+4.0
def f5(h3):  # delayed inhibition of hepatic glucose production
    return 180.0/(1.0+np.exp(0.29*h3/v1-7.5))

# equilibrium: h1*=h2*=h3*=x* (independent of t3), solve (x,y,z)
def equations(vars):
    x, y, z = vars
    e1 = f1(z) - E*(x/v1 - y/v2) - x/t1
    e2 = E*(x/v1 - y/v2) - y/t2
    e3 = f5(x) + I_infusion - f2(z) - f3(z)*f4(y)
    return [e1, e2, e3]

x0_guess = [90.0, 180.0, 13000.0]
sol = fsolve(equations, x0_guess, full_output=True, xtol=1e-13)
xeq, infodict, ier, msg = sol
print("equilibrium solve:", xeq, "ier=", ier, "msg=", msg)
residual = equations(xeq)
print("residual at solution:", residual)
X_EQ, Y_EQ, Z_EQ = xeq
print(f"X_EQ={X_EQ:.6f} Y_EQ={Y_EQ:.6f} Z_EQ={Z_EQ:.6f}")

# sanity: compare to the paper's reported ICs (x0=90,y0=180,z0=13000) -- if IC is near
# equilibrium that's expected (typical for these published sims); but the true limit cycle
# ORBITS around this equilibrium (which is unstable, since we know a sustained oscillation exists).

# ===== Build full 6-state RHS for arbitrary t3 (for direct nonlinear simulation cross-check) =====
def rhs(t, s, t3):
    x, y, z, h1, h2, h3 = s
    dx = f1(z) - E*(x/v1 - y/v2) - x/t1
    dy = E*(x/v1 - y/v2) - y/t2
    dz = f5(h3) + I_infusion - f2(z) - f3(z)*f4(y)
    dh1 = 3*(x-h1)/t3
    dh2 = 3*(h1-h2)/t3
    dh3 = 3*(h2-h3)/t3
    return np.array([dx,dy,dz,dh1,dh2,dh3])

# ===== analytic Jacobian at equilibrium (h1=h2=h3=X_EQ), as function of t3 =====
def analytic_jacobian(t3):
    x,y,z = X_EQ, Y_EQ, Z_EQ
    # f1'(z)
    u = np.exp(-z/3000.0+6.6)
    f1p = 209.0*u/(3000.0*(1+u)**2)
    # f2'(z)
    f2p = 0.05*np.exp(-z/1440.0)
    # f3'(z) constant
    f3p = 0.001
    # f4'(y)
    c = 1.0/v2 + 1.0/(E*t2)
    w = np.exp(7.76)*(c*y)**(-1.772)
    f4p = 90.0*1.772*w/(y*(1+w)**2)
    # f5'(h3) at h3=x (since h3_eq = x_eq)
    vv = np.exp(0.29*x/v1 - 7.5)
    f5p = -180.0*0.29/v1*vv/(1+vv)**2

    J = np.zeros((6,6))
    # row0: dx/dt
    J[0,0] = -E/v1 - 1.0/t1
    J[0,1] = E/v2
    J[0,2] = f1p
    # row1: dy/dt
    J[1,0] = E/v1
    J[1,1] = -E/v2 - 1.0/t2
    # row2: dz/dt
    J[2,1] = -f3(z)*f4p
    J[2,2] = -f2p - f3p*f4(y)
    J[2,5] = f5p
    # row3: dh1/dt
    J[3,0] = 3.0/t3
    J[3,3] = -3.0/t3
    # row4: dh2/dt
    J[4,3] = 3.0/t3
    J[4,4] = -3.0/t3
    # row5: dh3/dt
    J[5,4] = 3.0/t3
    J[5,5] = -3.0/t3
    return J

# ===== finite-difference Jacobian cross-check (at t3=36) =====
def fd_jacobian(t3, h=1e-4):
    s0 = np.array([X_EQ, Y_EQ, Z_EQ, X_EQ, X_EQ, X_EQ])
    n = 6
    J = np.zeros((6,6))
    scales = np.array([1.0,1.0,10.0,1.0,1.0,1.0])  # z is ~O(10000), use relative step there
    for j in range(n):
        step = h*max(1.0, abs(s0[j]))
        sp = s0.copy(); sp[j]+=step
        sm = s0.copy(); sm[j]-=step
        fp = rhs(0, sp, t3)
        fm = rhs(0, sm, t3)
        J[:,j] = (fp-fm)/(2*step)
    return J

J_analytic_36 = analytic_jacobian(36.0)
J_fd_36 = fd_jacobian(36.0)
diff = np.abs(J_analytic_36 - J_fd_36)
print("\nmax abs diff analytic vs finite-difference Jacobian @t3=36:", diff.max())
print("analytic J:\n", np.array2string(J_analytic_36, precision=6, suppress_small=True))
print("FD J:\n", np.array2string(J_fd_36, precision=6, suppress_small=True))

eigs36 = np.linalg.eigvals(J_analytic_36)
print("\neigenvalues at t3=36:", eigs36)
idx = np.argsort(-eigs36.real)
print("sorted by Re desc:", eigs36[idx])

OUT_ROOT = os.environ.get("BODYTWIN_OUT", os.path.join(os.getcwd(), "outputs"))
OUT_JSON = os.path.join(OUT_ROOT, "insulin_ultradian_sturis_hopf", "insulin_equilibrium.json")
os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
with open(OUT_JSON, 'w') as f:
    json.dump({
        "X_EQ": float(X_EQ), "Y_EQ": float(Y_EQ), "Z_EQ": float(Z_EQ),
        "residual": [float(r) for r in residual],
        "max_diff_analytic_vs_fd_jacobian": float(diff.max()),
        "eigs_at_t3_36_real": eigs36.real.tolist(), "eigs_at_t3_36_imag": eigs36.imag.tolist(),
    }, f, indent=2)
print(f"saved {OUT_JSON}")
