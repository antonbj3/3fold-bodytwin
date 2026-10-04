from dental_release.paths import expand as _release_expand
from bench import *
from force_operator import *
import time, resource

def queries(n):
    for i in range(n):
        for sign in [-1, 1]:
            h = np.zeros(n)
            h[i] = sign * 0.05
            yield ({'kind': 'single', 'site': i, 'sign': sign}, h)
        h = np.zeros(n)
        h[i] = 0.05
        h[(i + 1) % n] = -0.05
        yield ({'kind': 'paired', 'site': i, 'other': (i + 1) % n}, h)

def main():
    check_frozen('R3')
    start = time.perf_counter()
    cpu = time.process_time()
    panel = json.loads((ROOT / 'PREREG_R1.json').read_text())['source_panel']
    rows = []
    for case in panel:
        (ids, A, M, src) = geometry(case)
        n = len(ids)
        N = balanced_basis(A)
        k = 750 * (1 + 0.2 * np.cos(2 * np.pi * np.arange(n) / n))
        (w, q0) = base_w(A, k)
        (q, f0, K) = fixed(A, k, w, np.zeros(n))
        gold = np.diag(k) - k[:, None] * A @ np.linalg.solve(K, A.T * k)
        Y = []
        probes = []
        valid = True
        for l in range(N.shape[1]):
            h = 0.02 * N[:, l]
            (qp, fp, _) = fixed(A, k, w, h)
            (qm, fm, _) = fixed(A, k, w, -h)
            probes.append({'direction': l, 'height_plus_mm': h.tolist(), 'force_plus_N': fp.tolist(), 'force_minus_N': fm.tolist(), 'minimum_compression_mm': float(min((A @ qp + h).min(), (A @ qm - h).min()))})
            Y.append((fp - fm) / 0.04)
        Y = np.array(Y).T
        H = N.T @ Y
        H = (H + H.T) / 2
        J = N @ H @ N.T
        parity = float(np.max(np.abs(J - gold)) / np.max(np.abs(gold)))
        taskrows = []
        for (task, h) in queries(n):
            (qr, fr, meta) = unilateral(A, k, w, h, q0)
            linear = f0 + J @ h
            (fc, control) = dual_contact(N, H, f0, h)
            le = float(np.max(np.abs(linear - fr)))
            ce = float(np.max(np.abs(fc - fr)))
            balance = float(np.max(np.abs(A.T @ fc - w)))
            taskrows.append({'query': task, 'h_mm': h.tolist(), 'oracle_force_N': fr.tolist(), 'candidate_force_N': linear.tolist(), 'control_force_N': fc.tolist(), 'candidate_max_error_N': le, 'control_max_error_N': ce, 'candidate_min_force_N': float(linear.min()), 'n_released': int(np.sum(fr < 1e-07)), 'candidate_gate': bool(le <= 1 and linear.min() >= -1e-07 and (np.max(np.abs(A.T @ linear - w)) <= 1e-07)), 'control_gate': bool(ce <= 1 and fc.min() >= -1e-07 and (balance <= 1e-07)), 'oracle_numerical_gate': bool(meta['equilibrium_residual_N'] <= 1e-07), 'oracle': meta, 'control': control, 'control_balance_N': balance})
        rows.append({'case': case, 'n': n, 'source': src, 'basis_N': N.tolist(), 'H_N_per_mm': H.tolist(), 'J_N_per_mm': J.tolist(), 'A': A.tolist(), 'f0_N': f0.tolist(), 'w_N': w.tolist(), 'calibration': probes, 'tangent_relative_error': parity, 'calibration_gate': bool(parity <= 1e-08 and min((p['minimum_compression_mm'] for p in probes)) > 0), 'queries': taskrows})
    tasks = [t for r in rows for t in r['queries']]
    out = {'round': 'R3', 'claim_type': 'information_link', 'prereg_sha256': sha(ROOT / 'PREREG_R3.json'), 'n_cases': len(rows), 'n_queries': len(tasks), 'calibration_gate': all((r['calibration_gate'] for r in rows)), 'oracle_numerical_gate': all((t['oracle_numerical_gate'] for t in tasks)), 'candidate_gate': 'PASS' if all((t['candidate_gate'] for t in tasks)) else 'FAIL', 'candidate_failed_queries': sum((not t['candidate_gate'] for t in tasks)), 'control_gate': 'PASS' if all((t['control_gate'] for t in tasks)) else 'FAIL', 'max_candidate_error_N': max((t['candidate_max_error_N'] for t in tasks)), 'max_control_error_N': max((t['control_max_error_N'] for t in tasks)), 'queries_with_release': sum((t['n_released'] > 0 for t in tasks)), 'physical_measurement': 'NOT_RUN', 'matched_external_facit': 'ABSENT', 'rows': rows, 'wall_s': time.perf_counter() - start, 'cpu_s': time.process_time() - cpu, 'peak_rss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024}
    write(ROOT / 'raw/R3.json', out)
    print(json.dumps({k: v for (k, v) in out.items() if k != 'rows'}, indent=2))
    write(ROOT / 'CURRENT_WORK_STATE.json', {'lane': _release_expand('XBREAK_HUNT_1'), 'status': 'R3_COMPLETE', 'latest_gate': out['candidate_gate'], 'milestone': 'Measured force response reconstructs tangent; unrestricted edit model is tested against release', 'next_operation': 'Carry conventional dual contact into coupled supports with a metrology bound; preserve R3 failure', 'evidence': 'raw/R3.json'})
if __name__ == '__main__':
    main()
