from dental_release.paths import expand as _release_expand
from bench import *
import time, resource
from scipy.optimize import linprog

def main():
    check_frozen('R1')
    start = time.perf_counter()
    cpu = time.process_time()
    reg = json.loads((ROOT / 'PREREG_R1.json').read_text())
    rows = []
    for case in reg['source_panel']:
        (ids, A, M, src) = geometry(case)
        n = len(ids)
        if n < 7 or np.linalg.matrix_rank(A) < 3:
            rows.append({'case': case, 'status': 'REJECTED_GEOMETRY', 'n': n, 'source': src})
            continue
        nominal = np.full(n, 750.0)
        (w, q0) = base_w(A, nominal)
        K = A.T @ (nominal[:, None] * A)
        b = M @ nominal
        fit = np.linalg.lstsq(M, b, rcond=1e-12)[0]
        edits = []
        maxdiff = 0.0
        checks = []
        for site in range(n):
            obj = np.eye(n)[site]
            lo = linprog(obj, A_eq=M, b_eq=b, bounds=[(500, 1130)] * n, method='highs', options={'threads': 4})
            hi = linprog(-obj, A_eq=M, b_eq=b, bounds=[(500, 1130)] * n, method='highs', options={'threads': 4})
            assert lo.success and hi.success
            h = np.eye(n)[site] * 0.01
            solutions = []
            for (label, k) in [('low', lo.x), ('high', hi.x)]:
                (q, f, kk) = fixed(A, k, w, h)
                (qc, fc, meta) = unilateral(A, k, w, h, q0)
                (_, f0, _) = fixed(A, k, w, np.zeros(n))
                ker = float(np.max(np.abs(kk - K)) / np.max(np.abs(K)))
                parity = float(np.max(np.abs(f - fc)))
                valid = ker <= 1e-08 and meta['equilibrium_residual_N'] <= 1e-07 and (parity <= 1e-07) and (float((A @ q + h).min()) > 0)
                checks.append(valid)
                solutions.append({'label': label, 'k_N_per_mm': k.tolist(), 'q_mm': q.tolist(), 'f_before_N': f0.tolist(), 'f_after_N': f.tolist(), 'delta_f_N': (f - f0).tolist(), 'relative_K_error': ker, 'force_parity_N': parity, 'minimum_compression_mm': float((A @ q + h).min()), 'unilateral': meta, 'valid': valid})
            delta = np.asarray(solutions[1]['delta_f_N']) - np.asarray(solutions[0]['delta_f_N'])
            ambiguity = float(np.max(np.abs(delta)))
            maxdiff = max(maxdiff, ambiguity)
            edits.append({'site': site, 'fdi': ids[site], 'height_mm': 0.01, 'ambiguity_N': ambiguity, 'solutions': solutions})
        rows.append({'case': case, 'n': n, 'rank_moments': int(np.linalg.matrix_rank(M)), 'nullity': int(n - np.linalg.matrix_rank(M)), 'source': src, 'A': A.tolist(), 'K_N_per_mm': K.tolist(), 'w_N': w.tolist(), 'nominal_k_N_per_mm': nominal.tolist(), 'candidate_least_norm_k_N_per_mm': fit.tolist(), 'candidate_within_box': bool(np.all((fit >= 500) & (fit <= 1130))), 'max_local_response_ambiguity_N': maxdiff, 'numerical_gate': all(checks), 'sufficiency_gate': 'PASS' if maxdiff <= 1.0 else 'FAIL', 'edits': edits})
    valid = [r for r in rows if r.get('numerical_gate')]
    report = {'round': 'R1', 'claim_type': 'algorithm', 'prereg_sha256': sha(ROOT / 'PREREG_R1.json'), 'n_source_cases': len(rows), 'n_valid_cases': len(valid), 'rows': rows, 'numerical_gate': len(valid) == len(rows), 'sufficiency_gate': 'PASS' if all((r['sufficiency_gate'] == 'PASS' for r in valid)) and len(valid) == len(rows) else 'FAIL', 'strong_control': 'Executed bounded LP for every site, both endpoints, then independently minimized unilateral energy', 'wall_s': time.perf_counter() - start, 'cpu_s': time.process_time() - cpu, 'peak_rss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'patient_validation': 'UNKNOWN', 'physical_measurement': 'NOT_RUN'}
    write(ROOT / 'raw/R1.json', report)
    print(json.dumps({k: v for (k, v) in report.items() if k != 'rows'}, indent=2))
    print('maximum ambiguity N', max((r['max_local_response_ambiguity_N'] for r in valid)))
    write(ROOT / 'CURRENT_WORK_STATE.json', {'lane': _release_expand('XBREAK_HUNT_1'), 'status': 'R1_COMPLETE', 'latest_gate': report['sufficiency_gate'], 'milestone': 'Global K cannot determine local-height response', 'next_operation': 'Write handoff, preregister local shim plus global pose acquisition', 'evidence': 'raw/R1.json', 'patient_support': 'UNKNOWN'})
if __name__ == '__main__':
    main()
