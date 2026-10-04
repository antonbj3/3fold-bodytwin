"""Same measured information: conventional box-constrained affine response LP.
Supplementary stronger control, constructed after R4 primary prereg;
not used to change any primary metric or tolerance.
"""
import hashlib, json, time, warnings
from pathlib import Path
import numpy as np
from scipy.optimize import linprog, OptimizeWarning
R = Path(__file__).resolve().parents[1]

def solve(s, h, channel, corrupt=False):
    n = len(s['baseline_force_N'])
    A = np.array(s['A'])
    Y = np.array(s['raw_probe_force_N'])
    eps = s['raw_channel_error_N']
    step = s['probe_height_mm']
    if eps < 0 or not np.isfinite(eps):
        raise ValueError('Invalid raw error')
    if corrupt:
        Y = Y.copy()
        Y[0] += 3.0
    I = np.eye(n)
    D = np.r_[np.c_[I, np.zeros_like(I)], np.c_[I, -I], np.c_[I, I]]
    c = np.zeros(2 * n)
    c[channel] = 1
    c[n + channel] = h / step
    E = np.block([[A.T, np.zeros_like(A.T)], [np.zeros_like(A.T), A.T]])
    kwargs = dict(A_ub=np.r_[D, -D], b_ub=np.r_[Y.ravel() + eps, -Y.ravel() + eps], A_eq=E, b_eq=np.r_[s['w_N'], [0, 0, 0]], bounds=[(None, None)] * (2 * n), method='highs', options={'threads': 1})
    lo = linprog(c, **kwargs)
    if not lo.success:
        return None
    hi = linprog(-c, **kwargs)
    if not hi.success:
        return None
    return [float(lo.fun), float(-hi.fun)]

def main():
    warnings.filterwarnings('ignore', category=OptimizeWarning, message='Unrecognized options detected.*')
    t = time.perf_counter()
    rows = json.loads((R / 'raw/R4_CONTROLS.json').read_text())
    records = []
    calls = 0
    for row in rows:
        s = json.loads((R / 'exports' / f"{row['case']}_corridor_state.json").read_text())
        bounds = []
        for i in range(len(s['baseline_force_N'])):
            interval = solve(s, row['height_mm'], i)
            calls += 2
            if interval is None:
                raise ArithmeticError('Measured same-information affine response LP infeasible')
            bounds.append(interval)
        fc = np.array(row['true_force_N'])
        bounds = np.array(bounds)
        records.append({'case': row['case'], 'height_mm': row['height_mm'], 'intervals_N': bounds.tolist(), 'covered': bool(np.all(fc >= bounds[:, 0] - 1e-06) and np.all(fc <= bounds[:, 1] + 1e-06)), 'max_width_N': float(np.max(bounds[:, 1] - bounds[:, 0])), 'candidate_joint_radius_N': row['radius_N']})
    s = json.loads((R / 'exports/F5367_corridor_state.json').read_text())
    fault = solve(s, 0, 0, corrupt=True) is None
    calls += 1
    out = {'scope': 'Supplementary stronger control after primary prereg. Same raw measured force vectors, error box, wrench and linear full-contact closure; no support model provided.', 'claim_type': 'capability', 'queries': len(records), 'contained': sum((r['covered'] for r in records)), 'max_width_N': max((r['max_width_N'] for r in records)), 'injected_raw_force_corruption_detected': fault, 'LP_calls': calls, 'wall_s': time.perf_counter() - t, 'records': records}
    (R / 'raw/R4_STRONGER_LP_CONTROL.json').write_text(json.dumps(out, indent=2))
    print(json.dumps({k: out[k] for k in ['queries', 'contained', 'max_width_N', 'injected_raw_force_corruption_detected', 'LP_calls', 'wall_s']}))
if __name__ == '__main__':
    main()
