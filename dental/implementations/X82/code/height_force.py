"""X82 conditional force queries. No individual stiffness is fitted.
Projected energy and bound: XBREAK_HUNT_1/DERIVATION.md.
Numerical solve is not a formal floating-point certificate.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from active_contact import solve

def nonnegative_scalar(value, name):
    if isinstance(value, bool):
        raise ValueError(name + ' must be a finite nonnegative scalar')
    try:
        x = float(value)
    except (TypeError, ValueError):
        raise ValueError(name + ' must be a finite nonnegative scalar') from None
    if not np.isfinite(x) or x < 0:
        raise ValueError(name + ' must be finite and nonnegative')
    return x

def array(value, name):
    x = np.asarray(value, dtype=float)
    if not np.isfinite(x).all():
        raise ValueError(name + ' must be finite')
    return x

def validate(state):
    if state.get('schema') != 'conditional-edit-response-v1':
        raise ValueError('Unsupported response schema')
    rho = nonnegative_scalar(state['spectral_error_N_per_mm'], 'spectral_error_N_per_mm')
    eta = nonnegative_scalar(state['baseline_error_N'], 'baseline_error_N')
    for key in ['old_spectral_error_N_per_mm', 'activation_column_error_N_per_mm', 'load_response_error_l2_N_per_N', 'load_change_limit_N']:
        if key in state:
            nonnegative_scalar(state[key], key)
    A = array(state['A'], 'A')
    N = array(state['basis_N'], 'basis_N')
    H = array(state['H_N_per_mm'], 'H_N_per_mm')
    f0 = array(state['baseline_force_N'], 'baseline_force_N')
    w = array(state['w_N'], 'w_N')
    hb = array(state['baseline_height_mm'], 'baseline_height_mm')
    n = len(f0)
    if not 4 <= n <= 16 or A.shape != (n, 3) or N.shape != (n, n - 3) or (H.shape != (n - 3, n - 3)):
        raise ValueError('Response dimensions must match 4..16 supports and three wrench coordinates')
    if w.shape != (3,) or hb.shape != (n,):
        raise ValueError('Reference dimensions do not match')
    if len(state['predicted_fdi']) != n or len(set(state['predicted_fdi'])) != n:
        raise ValueError('FDI mapping must be complete and unique')
    if not np.all(A[:, 0] == 1.0) or np.linalg.matrix_rank(A) != 3:
        raise ValueError('A requires total axial force and two independent scaled moments')
    if np.linalg.norm(N.T @ A) > 1e-07 or np.linalg.norm(N.T @ N - np.eye(n - 3)) > 1e-07:
        raise ValueError('Invalid force-balanced orthonormal basis')
    if np.max(np.abs(H - H.T)) > 1e-07:
        raise ValueError('Response must be reciprocal')
    if f0.min() <= 0 or np.max(np.abs(A.T @ f0 - w)) > 1e-07:
        raise ValueError('Reference requires positive forces balanced against the recorded wrench')
    if not state.get('case') or not state.get('geometry_sha256'):
        raise ValueError('Case and geometry identity required')
    return (A, N, H, f0, w, hb, rho, eta)

def at_load(state, load_N):
    """Recenter only along the measured wrench ray, preserving its force energy."""
    (A, N, H, f0, w, hb, rho, eta) = validate(state)
    T = nonnegative_scalar(load_N, 'load_N')
    dT = T - w[0]
    if T <= 0:
        raise ValueError('Strictly positive applied load required')
    if dT == 0:
        return state
    if 'load_response_N_per_N' not in state:
        raise LookupError('Missing measured regional load-response column')
    limit = nonnegative_scalar(state['load_change_limit_N'], 'load_change_limit_N')
    ep = nonnegative_scalar(state['load_response_error_l2_N_per_N'], 'load_response_error_l2_N_per_N')
    if abs(dT) > limit + 1e-10:
        raise LookupError('Load outside measured load ray domain')
    p = array(state['load_response_N_per_N'], 'load_response_N_per_N')
    if p.shape != f0.shape or np.max(np.abs(A.T @ p - w / w[0])) > 1e-07:
        raise ValueError('Load-response column must balance recorded load ray')
    fref = f0 + p * dT
    if fref.min() <= 0:
        raise LookupError('Recentered full-contact reference is not positive; this implementation requires a new loaded reference')
    out = dict(state, baseline_force_N=fref.tolist(), w_N=(w * T / w[0]).tolist(), baseline_error_N=eta + ep * abs(dT), original_baseline_force_N=f0.tolist(), original_baseline_error_N=eta, original_load_N=float(w[0]), load_change_N=float(dT))
    return out

def predict(state, edit):
    if 'load_N' in edit:
        try:
            state = at_load(state, edit['load_N'])
        except LookupError as exc:
            return {'status': 'UNKNOWN', 'reason': str(exc), 'force_interval_N': None}
    (A, N, H, f0, w, hb, rho, eta) = validate(state)
    if edit.get('unit') != 'mm':
        raise ValueError('Height unit must be mm')
    if edit.get('case') != state['case'] or edit.get('geometry_sha256') != state['geometry_sha256']:
        raise ValueError('Case/geometry identity mismatch')
    if 'w_N' in edit and (not np.array_equal(array(edit['w_N'], 'w_N'), w)):
        return {'status': 'UNKNOWN', 'reason': 'A new wrench requires additional response measurements', 'force_interval_N': None}
    h = array(edit['height_change_mm'], 'height_change_mm')
    if h.shape != f0.shape or np.max(np.abs(h)) > 0.05 + 1e-12:
        raise ValueError('One height per support required within +/-0.05 mm')
    if not state.get('complete_active_calibration'):
        return {'status': 'UNKNOWN', 'reason': 'Candidate contact coverage is incomplete', 'force_interval_N': None}
    if state.get('radius_model') != 'full_response':
        return {'status': 'UNKNOWN', 'reason': 'Directional activation states require their separate verified error contract', 'force_interval_N': None}
    eig = np.linalg.eigvalsh(H)
    low = float(eig[0] - rho)
    high = float(eig[-1] + rho)
    if low <= 0:
        return {'status': 'UNKNOWN', 'reason': 'Uncertainty exceeds the positive response eigenvalue', 'force_interval_N': None}
    G = np.linalg.inv(H)
    C = N @ G @ N.T + (np.eye(len(f0)) - N @ N.T) / 750.0
    (f, pose, kkt) = solve(A, C, -C @ f0, w, h, f0)
    z = N.T @ (f - f0)
    R = float(np.sqrt(high / low) * (rho * np.linalg.norm(G @ z) + eta) + 1e-06)
    if not np.isfinite(R) or R < 0:
        raise ArithmeticError('Invalid error radius')
    numeric = kkt['balance_residual_N'] <= 1e-07 and kkt['min_force_N'] >= -1e-07 and (kkt['min_slack_mm'] >= -1e-07) and (kkt['complementarity_N_mm'] <= 1e-07) and (kkt['active_stationarity_mm'] <= 1e-07)
    if not numeric:
        return {'status': 'UNKNOWN', 'reason': 'Contact certificate failed', 'force_interval_N': None, 'KKT': kkt}
    f = np.maximum(f, 0)
    intervals = np.column_stack((np.maximum(0, f - R), f + R))
    original_f0 = np.array(state.get('original_baseline_force_N', f0))
    original_eta = state.get('original_baseline_error_N', eta)
    delta = f - original_f0
    delta_R = R + original_eta
    result = {'status': 'CONDITIONAL_PASS' if R <= 1 else 'UNKNOWN', 'case': state['case'], 'claim_type': 'capability', 'resolution': 'PER_TOOTH', 'timescale': 'SIMULTANEOUS', 'predicted_fdi': state['predicted_fdi'], 'force_N': f.tolist(), 'force_interval_N': intervals.tolist(), 'force_change_N': delta.tolist(), 'force_change_interval_N': np.column_stack((delta - delta_R, delta + delta_R)).tolist(), 'joint_l2_error_radius_N': R, 'delta_joint_l2_error_radius_N': delta_R, 'KKT': kkt, 'height_change_mm': h.tolist(), 'w_N': w.tolist(), 'joint_constraint': 'A.T f = w; marginal intervals are not independent force boxes', 'enclosure_kind': 'Analytic deterministic model bound; 1e-6 N numerical reserve. No formal floating-point enclosure.', 'physical_validity': 'UNKNOWN; ' + state.get('provenance_kind', 'missing provenance'), 'measurement_contract': state.get('measurement_contract', 'UNKNOWN')}
    if not np.isfinite(intervals).all() or np.any(intervals[:, 0] > intervals[:, 1]):
        raise ArithmeticError('Unordered or nonfinite interval')
    return result

def make_edit(state, index, height):
    h = np.zeros(len(state['baseline_force_N']))
    h[index] = height
    return {'case': state['case'], 'geometry_sha256': state['geometry_sha256'], 'unit': 'mm', 'height_change_mm': h.tolist()}

def inverse(state, fdi, target_N, target_halfwidth_N=1.0, load_N=None):
    if load_N is not None:
        try:
            state = at_load(state, load_N)
        except LookupError as exc:
            return {'status': 'UNKNOWN', 'reason': str(exc), 'query_calls': 0}
    validate(state)
    target = nonnegative_scalar(target_N, 'target_N')
    tol = nonnegative_scalar(target_halfwidth_N, 'target_halfwidth_N')
    if fdi not in state['predicted_fdi']:
        raise ValueError('Target tooth missing from measured contacts')
    i = state['predicted_fdi'].index(fdi)
    calls = 0

    def evaluate(h):
        nonlocal calls
        calls += 1
        ans = predict(state, make_edit(state, i, h))
        if ans.get('force_interval_N') is None:
            raise LookupError(ans.get('reason', 'Missing response'))
        return ans
    try:
        (left, right) = (-0.05, 0.05)
        al = evaluate(left)
        ar = evaluate(right)
        if target < al['force_interval_N'][i][0] or target > ar['force_interval_N'][i][1]:
            return {'status': 'INFEASIBLE_IN_CALIBRATED_BOX', 'query_calls': calls, 'target_N': target, 'fdi': fdi}
        if target < al['force_N'][i] or target > ar['force_N'][i]:
            return {'status': 'UNKNOWN', 'reason': 'Only uncertainty overlaps a height-box endpoint', 'query_calls': calls}
        while right - left > 1e-08:
            mid = (left + right) / 2
            am = evaluate(mid)
            if am['force_N'][i] < target:
                left = mid
            else:
                right = mid
        h = (left + right) / 2
        ans = evaluate(h)
        (lo, hi) = ans['force_interval_N'][i]
        guaranteed = lo >= target - tol and hi <= target + tol and (ans['status'] == 'CONDITIONAL_PASS')
        brackets = []
        for side in ['upper', 'lower']:
            (l, r) = (-0.05, 0.05)
            il = 1 if side == 'upper' else 0
            fl = evaluate(l)['force_interval_N'][i][il] - target
            fr = evaluate(r)['force_interval_N'][i][il] - target
            if fl > 0 or fr < 0:
                brackets.append(None)
                continue
            while r - l > 1e-08:
                m = (l + r) / 2
                fm = evaluate(m)['force_interval_N'][i][il] - target
                if fm <= 0:
                    l = m
                else:
                    r = m
            brackets.append(l if side == 'upper' else r)
        root_bracket = brackets if all((x is not None for x in brackets)) else None
        return {'status': 'CONDITIONAL_TARGET_BAND' if guaranteed else 'UNKNOWN', 'fdi': fdi, 'target_N': target, 'target_band_N': [max(0, target - tol), target + tol], 'nominal_height_change_mm': h, 'possible_exact_target_height_bracket_mm': root_bracket, 'height_bracket_kind': 'Conditional analytic model enclosure; numerical rounding not formally enclosed', 'at_selected_height': ans, 'query_calls': calls, 'resolution': 'PER_TOOTH', 'physical_validity': 'UNKNOWN; laboratory scenario, no clinical recommendation'}
    except LookupError as exc:
        return {'status': 'UNKNOWN', 'reason': str(exc), 'query_calls': calls}

def project_reaction(y0, yv, h, d, de, fe, matched_velocity=False, matched_history=False, position_closure=False):
    """Repaired NOSEARCH R6; preserve its conditional position/rate contract."""
    de = nonnegative_scalar(de, 'displacement_error_mm')
    fe = nonnegative_scalar(fe, 'reaction_error_N')
    for (name, value) in [('y0', y0), ('yv', yv), ('h', h), ('d', d)]:
        if not np.isfinite(value):
            raise ValueError(name + ' must be finite')
    if not (matched_velocity and matched_history and position_closure):
        return {'status': 'UNKNOWN', 'force_interval_N': None, 'reason': 'Unverified query-state/velocity/history or position closure'}
    if h <= 0 or d - de <= h or y0 - fe < 0:
        return {'status': 'UNKNOWN', 'force_interval_N': None, 'reason': 'Unsupported projection regime'}
    lo = ((d - de) / h - 1) * (y0 - fe) + (yv - fe)
    hi = ((d + de) / h - 1) * (y0 + fe) + (yv + fe)
    if not np.isfinite([lo, hi]).all() or lo > hi:
        raise ArithmeticError('Invalid interval')
    return {'status': 'CONDITIONAL_POINT_BOUND', 'force_interval_N': [lo, hi]}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('state')
    ap.add_argument('--fdi', type=int, required=True)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--height-mm', type=float)
    g.add_argument('--target-N', type=float)
    ap.add_argument('--load-N', type=float)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    s = json.loads(Path(args.state).read_text())
    if args.height_mm is not None:
        edit = make_edit(s, s['predicted_fdi'].index(args.fdi), args.height_mm)
        if args.load_N is not None:
            edit['load_N'] = args.load_N
        ans = predict(s, edit)
    else:
        ans = inverse(s, args.fdi, args.target_N, load_N=args.load_N)
    Path(args.out).write_text(json.dumps(ans, indent=2, allow_nan=False))
    print(json.dumps({'status': ans['status'], 'out': args.out}))
if __name__ == '__main__':
    main()
