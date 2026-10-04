import numpy as np
from scipy.optimize import linprog
from .regional import predict

def inverse_signed(model, B, weights, allow, extra_A=None, extra_b=None):
    J = np.array(model['J_N_per_mm'])
    f0 = np.array(model['f0'])
    eta = model['eta_N']
    beta = np.array(model['beta_N_per_mm'])
    D = np.c_[J, -J]
    br = np.r_[beta, beta]
    G = []
    b = []
    for i in range(2):
        G.extend([D[i] + br, -D[i] + br])
        b.extend([9 - f0[i] - eta, f0[i] - eta - 7])
    G.extend(-D + br)
    b.extend(f0 - eta - 1e-06)
    G.extend(np.c_[-B, B])
    b.extend(allow)
    if extra_A is not None:
        G.extend(np.c_[extra_A, -extra_A])
        b.extend(extra_b)
    sol = linprog(np.r_[weights @ B, weights @ B], A_ub=np.array(G), b_ub=np.array(b), bounds=[(0, 0.05)] * 4, method='highs', options={'threads': 1})
    if not sol.success:
        return dict(status='INFEASIBLE_OR_UNKNOWN', reason=sol.message)
    a = sol.x[:2] - sol.x[2:]
    p = predict(model, a)
    iv = np.array(p['force_interval_N'])
    good = p['positive_contact_path'] and np.all(iv[:2, 0] >= 7 - 1e-07) and np.all(iv[:2, 1] <= 9 + 1e-07)
    return dict(status='CONDITIONAL_SIGNED_TARGET' if good else 'UNKNOWN', coefficients_mm=a, at_selected_height=p, positive_height_requires_new_manufacture=bool(np.any(a > 0)), mean_absolute_design_change_mm=float(sol.fun / sum(weights)))
