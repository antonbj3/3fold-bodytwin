"""Consume a calibrated force-response state; never fit hidden support values."""
import argparse, json, time, numpy as np
from pathlib import Path
from bench import write, sha
from active_contact import solve
from run_r6 import directional_radii

def predict(state, edit, policy='corner'):
    if state['schema'] != 'conditional-edit-response-v1' or edit['unit'] != 'mm':
        raise ValueError('Expected conditional-edit-response-v1 and height unit mm')
    if edit['case'] != state['case'] or edit['geometry_sha256'] != state['geometry_sha256']:
        raise ValueError('Case/geometry identity mismatch')
    if not state.get('complete_active_calibration'):
        return {'status': 'UNKNOWN', 'reason': 'Initial candidate contacts were not all measured on the calibrated loaded branch'}
    A = np.array(state['A'])
    N = np.array(state['basis_N'])
    H = np.array(state['H_N_per_mm'])
    fbase = np.array(state['baseline_force_N'])
    w = np.array(state['w_N'])
    h = np.array(edit['height_change_mm'])
    n = len(fbase)
    hb = np.array(state['baseline_height_mm'])
    rho = state['spectral_error_N_per_mm']
    if n > 16:
        return {'status': 'UNKNOWN', 'reason': 'State exceeds the frozen 16-support query/control domain'}
    if h.shape != (n,) or not np.isfinite(h).all():
        raise ValueError('One finite height value required per support')
    if np.max(np.abs(h)) > 0.05 + 1e-12:
        raise ValueError('Edit outside frozen +/-0.05 mm query domain')
    if A.shape != (n, 3) or N.shape != (n, n - 3) or H.shape != (n - 3, n - 3):
        raise ValueError('Force/pose dimension mismatch')
    if np.linalg.norm(N.T @ A) > 1e-07 or np.linalg.norm(N.T @ N - np.eye(n - 3)) > 1e-07:
        raise ValueError('Invalid balanced basis')
    if np.max(np.abs(H - H.T)) > 1e-07:
        raise ValueError('Response must be reciprocal within declared numerical tolerance')
    if fbase.min() <= 0 or np.max(np.abs(A.T @ fbase - w)) > 1e-07:
        raise ValueError('Loaded reference force must be positive and balance the registered wrench')
    ev = np.linalg.eigvalsh(H)
    lower = ev.min() - rho
    upper = ev.max() + rho
    if lower <= 0:
        return {'status': 'UNKNOWN', 'reason': 'Response noise exceeds positive identifiable tangent'}
    G = np.linalg.inv(H)
    C = N @ G @ N.T + (np.eye(n) - N @ N.T) / 750
    g = -C @ fbase
    (f, q, kkt) = solve(A, C, g, w, h - hb, fbase)
    z = N.T @ (f - fbase)
    if state['radius_model'] == 'full_response':
        radius = float(np.sqrt(upper / lower) * (rho * np.linalg.norm(G @ z) + state['baseline_error_N']) + 1e-06)
    else:
        v = np.array(state['activation_column_N_per_mm'])
        corners = ((np.arange(2 ** n, dtype=np.uint32)[:, None] >> np.arange(n) & 1).astype(float) * 2 - 1) * 0.5
        (rc, rs, detail) = directional_radii(N, H, f, fbase, v, state['old_spectral_error_N_per_mm'], rho, corners, corners @ v)
        radius = rs if policy == 'corner' else rc
    numeric = kkt['balance_residual_N'] <= 1e-07 and kkt['min_force_N'] >= -1e-07 and (kkt['min_slack_mm'] >= -1e-07)
    return {'status': 'CONDITIONAL_PASS' if radius <= 1 and numeric else 'UNKNOWN', 'case': state['case'], 'predicted_fdi': state['predicted_fdi'], 'force_N': f.tolist(), 'force_interval_N': [[max(0, float(x - radius)), float(x + radius)] for x in f], 'joint_l2_error_radius_N': radius, 'radius_limit_N': 1.0, 'KKT': kkt, 'physical_validity': 'UNKNOWN; state provenance is ' + state['provenance_kind'], 'measurement_contract': state['measurement_contract'], 'height_change_mm': h.tolist()}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('state')
    ap.add_argument('edit')
    ap.add_argument('--out', required=True)
    ap.add_argument('--radius-policy', choices=['corner', 'triangle'], default='corner')
    args = ap.parse_args()
    start = time.perf_counter()
    s = json.loads(Path(args.state).read_text())
    e = json.loads(Path(args.edit).read_text())
    result = predict(s, e, args.radius_policy)
    result.update(state_sha256=sha(args.state), edit_sha256=sha(args.edit), wall_s=time.perf_counter() - start)
    write(args.out, result)
    print(json.dumps({'status': result['status'], 'radius_N': result.get('joint_l2_error_radius_N'), 'out': args.out}))
if __name__ == '__main__':
    main()
