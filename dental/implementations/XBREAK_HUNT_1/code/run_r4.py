from dental_release.paths import expand as _release_expand
from bench import *
from force_operator import balanced_basis, dual_contact
from active_contact import solve
from run_r3 import queries
import time, resource

def main():
    check_frozen('R4')
    start = time.perf_counter()
    cpu = time.process_time()
    panel = json.loads((ROOT / 'PREREG_R1.json').read_text())['source_panel']
    rows = []
    eps = 0.01
    s = 0.02
    for case in panel:
        (ids, A, M, src) = geometry(case)
        n = len(ids)
        N = balanced_basis(A)
        r = N.shape[1]
        xy = np.array(src['xy_mm'])
        k = 750 * (1 + 0.2 * np.cos(2 * np.pi * np.arange(n) / n))
        dist2 = np.sum((xy[:, None] - xy[None, :]) ** 2, axis=2)
        C = np.diag(1 / k) + 0.8 / 750 * np.exp(-dist2 / (2 * 15 ** 2))
        f0 = k * 300 / k.sum()
        w = A.T @ f0
        g = -C @ f0
        Htrue = np.linalg.inv(N.T @ C @ N)
        Jtrue = N @ Htrue @ N.T
        Y = []
        probes = []
        for l in range(r):
            h = s * N[:, l]
            (fp, qp, mp) = solve(A, C, g, w, h, f0)
            (fm, qm, mm) = solve(A, C, g, w, -h, f0)
            noise = eps * np.sin((np.arange(n) + 1) * (l + 2))
            measured_plus = fp + noise
            measured_minus = fm - noise
            Y.append((measured_plus - measured_minus) / (2 * s))
            probes.append({'direction': l, 'h_plus_mm': h.tolist(), 'force_true_plus_N': fp.tolist(), 'force_true_minus_N': fm.tolist(), 'raw_noise_plus_N': noise.tolist(), 'force_observed_plus_N': measured_plus.tolist(), 'force_observed_minus_N': measured_minus.tolist(), 'minimum_true_force_N': float(min(fp.min(), fm.min())), 'plus_KKT': mp, 'minus_KKT': mm})
        Y = np.array(Y).T
        Hhat = N.T @ Y
        Hhat = (Hhat + Hhat.T) / 2
        rho = eps / s * np.sqrt(n * r)
        eta0 = eps * np.sqrt(n)
        evalues = np.linalg.eigvalsh(Hhat)
        lower = float(evalues.min() - rho)
        upper = float(evalues.max() + rho)
        e0 = eps * np.sin(np.arange(n) + 1)
        fhat0 = f0 + N @ (N.T @ e0)
        Ghat = np.linalg.inv(Hhat)
        Q = np.eye(n) - N @ N.T
        Chat = N @ Ghat @ N.T + Q / 750
        ghat = -Chat @ fhat0
        Cprime = C + 0.0002 * A @ A.T
        gprime = g - 0.0002 * A @ w
        acquisition_ok = lower > 0 and min((p['minimum_true_force_N'] for p in probes)) > 0
        actual_error = float(np.linalg.norm(Hhat - Htrue, 2))
        taskrows = []
        for (task, h) in queries(n):
            (fr, qr, oracle) = solve(A, C, g, w, h, f0)
            (fa, qa, meta) = solve(A, Chat, ghat, w, h, fhat0)
            (fc, control) = dual_contact(N, Hhat, fhat0, h)
            (fprime, qprime, prime) = solve(A, Cprime, gprime, w, h, f0)
            z = N.T @ (fa - fhat0)
            radius = float(np.sqrt(upper / lower) * (rho * np.linalg.norm(Ghat @ z) + eta0) + 1e-06) if lower > 0 else None
            error = float(np.max(np.abs(fa - fr)))
            l2 = float(np.linalg.norm(fa - fr))
            parity = float(np.max(np.abs(fa - fc)))
            invariance = float(np.max(np.abs(fr - fprime)))
            numeric = meta['balance_residual_N'] <= 1e-07 and meta['active_stationarity_mm'] <= 1e-07 and (meta['min_force_N'] >= -1e-07) and (meta['min_slack_mm'] >= -1e-07) and (oracle['balance_residual_N'] <= 1e-07)
            taskrows.append({'query': task, 'h_mm': h.tolist(), 'oracle_force_N': fr.tolist(), 'candidate_force_N': fa.tolist(), 'control_force_N': fc.tolist(), 'radius_N': radius, 'max_force_error_N': error, 'l2_force_error_N': l2, 'radius_contains_truth': bool(radius is not None and l2 <= radius + 1e-06), 'radius_useful': bool(radius is not None and radius <= 1), 'n_released': int(np.sum(fr < 1e-08)), 'numerical_gate': bool(numeric), 'control_parity_N': parity, 'unobservable_compliance_force_parity_N': invariance, 'candidate_KKT': meta, 'oracle_KKT': oracle, 'control_KKT': control})
        rows.append({'case': case, 'source': src, 'n': n, 'n_probes': 1 + 2 * r, 'A': A.tolist(), 'basis_N': N.tolist(), 'truth_C_mm_per_N': C.tolist(), 'truth_gap_mm': g.tolist(), 'truth_H_N_per_mm': Htrue.tolist(), 'observed_H_N_per_mm': Hhat.tolist(), 'true_f0_N': f0.tolist(), 'observed_f0_N': fhat0.tolist(), 'baseline_raw_error_N': e0.tolist(), 'w_N': w.tolist(), 'rho_N_per_mm': float(rho), 'eta0_N': float(eta0), 'spectral_lower_N_per_mm': lower, 'spectral_upper_N_per_mm': upper, 'actual_spectral_error_N_per_mm': actual_error, 'acquisition_gate': bool(acquisition_ok and actual_error <= rho + 1e-07), 'calibration': probes, 'queries': taskrows, 'different_C_relative_frobenius_norm': float(np.linalg.norm(Cprime - C) / np.linalg.norm(C))})
    tasks = [t for r in rows for t in r['queries']]
    gate = all((r['acquisition_gate'] for r in rows)) and all((t['radius_contains_truth'] and t['radius_useful'] and t['numerical_gate'] for t in tasks))
    out = {'round': 'R4', 'claim_type': 'information_link', 'realization': 'SIMULATED_PROSPECTIVE_ACQUISITION', 'prereg_sha256': sha(ROOT / 'PREREG_R4.json'), 'n_cases': len(rows), 'n_queries': len(tasks), 'acquisition_gate': all((r['acquisition_gate'] for r in rows)), 'capability_gate': 'PASS' if gate else 'FAIL', 'truth_containment_gate': all((t['radius_contains_truth'] for t in tasks)), 'numerical_gate': all((t['numerical_gate'] for t in tasks)), 'useful_radius_queries': sum((t['radius_useful'] for t in tasks)), 'max_radius_N': max((t['radius_N'] for t in tasks if t['radius_N'] is not None)), 'max_actual_force_error_N': max((t['max_force_error_N'] for t in tasks)), 'control_gate': 'PASS' if all((t['control_parity_N'] <= 0.001 for t in tasks)) else 'FAIL', 'max_control_parity_N': max((t['control_parity_N'] for t in tasks)), 'max_unobservable_force_parity_N': max((t['unobservable_compliance_force_parity_N'] for t in tasks)), 'queries_with_release': sum((t['n_released'] > 0 for t in tasks)), 'physical_measurements': 'NOT_RUN', 'empirical_validity': 'UNKNOWN; external matched facit absent', 'rows': rows, 'wall_s': time.perf_counter() - start, 'cpu_s': time.process_time() - cpu, 'peak_rss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024}
    write(ROOT / 'raw/R4.json', out)
    print(json.dumps({k: v for (k, v) in out.items() if k != 'rows'}, indent=2))
    write(ROOT / 'CURRENT_WORK_STATE.json', {'lane': _release_expand('XBREAK_HUNT_1'), 'status': 'R4_COMPLETE', 'latest_gate': out['capability_gate'], 'milestone': 'Coupled support force energy and query-dependent metrology bound executed', 'next_operation': 'Attack initially inactive support coverage; baseline full contact is still load-bearing', 'evidence': 'raw/R4.json', 'empirical_validity': 'UNKNOWN'})
if __name__ == '__main__':
    main()
