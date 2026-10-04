from dental_release.paths import expand as _release_expand
from bench import *
from force_operator import balanced_basis, dual_contact
from active_contact import solve
from run_r3 import queries
import time, resource

def main():
    check_frozen('R5')
    start = time.perf_counter()
    cpu = time.process_time()
    panel = json.loads((ROOT / 'PREREG_R1.json').read_text())['source_panel']
    rows = []
    eps = 0.01
    s = 0.02
    delta_height = 0.04
    for case in panel:
        (ids, A, M, src) = geometry(case)
        n = len(ids)
        site = n - 1
        active = np.arange(n - 1)
        m = n - 1
        N = balanced_basis(A)
        Nold = np.zeros((n, m - 3))
        Nold[active] = balanced_basis(A[active])
        xy = np.array(src['xy_mm'])
        dist2 = np.sum((xy[:, None] - xy[None, :]) ** 2, axis=2)
        k = 750 * (1 + 0.2 * np.cos(2 * np.pi * np.arange(n) / n))
        f0 = k.copy()
        f0[site] = 0
        f0 *= 300 / f0.sum()
        w = A.T @ f0
        worlds = []
        for ki in [500.0, 1130.0]:
            kk = k.copy()
            kk[site] = ki
            C = np.diag(1 / kk) + 0.8 / 750 * np.exp(-dist2 / (2 * 15 ** 2))
            g = -C @ f0
            g[site] += 0.04
            oldY = []
            oldYclean = []
            oldprobe = []
            oldbranches = True
            all_oldforces = []
            for l in range(m - 3):
                h = s * Nold[:, l]
                (fp, qp, mp) = solve(A, C, g, w, h, f0)
                (fm, qm, mm) = solve(A, C, g, w, -h, f0)
                noise = eps * np.sin((np.arange(n) + 1) * (l + 2))
                noise[site] = 0
                oldY.append((fp + noise - fm + noise) / (2 * s))
                oldYclean.append((fp - fm) / (2 * s))
                oldbranches = oldbranches and fp[active].min() > 0 and (fm[active].min() > 0) and (abs(fp[site]) < 1e-07) and (abs(fm[site]) < 1e-07)
                all_oldforces.extend([fp.tolist(), fm.tolist()])
                oldprobe.append({'direction': l, 'f_plus_N': fp.tolist(), 'f_minus_N': fm.tolist(), 'min_old_active_force_N': float(min(fp[active].min(), fm[active].min()))})
            Y = np.array(oldY).T
            Yclean = np.array(oldYclean).T
            Ho = Nold.T @ Y
            Ho = (Ho + Ho.T) / 2
            Jold = Nold @ Ho @ Nold.T
            Hoc = Nold.T @ Yclean
            Hoc = (Hoc + Hoc.T) / 2
            Joldclean = Nold @ Hoc @ Nold.T
            small = []
            for sign in [-1, 1]:
                h = np.eye(n)[site] * s * sign
                (ff, qq, meta) = solve(A, C, g, w, h, f0)
                small.append(ff.tolist())
                all_oldforces.append(ff.tolist())
            points = []
            for (l, height) in enumerate([0.05, 0.09]):
                h = np.eye(n)[site] * height
                (ff, qq, meta) = solve(A, C, g, w, h, f0)
                noise = eps * np.sin((np.arange(n) + 1) * (l + 3))
                obs = ff + noise
                points.append({'height_mm': height, 'h': h, 'true_force_N': ff, 'observed_force_N': obs, 'KKT': meta, 'noise_N': noise})
            v = (points[1]['observed_force_N'] - points[0]['observed_force_N']) / delta_height
            vclean = (points[1]['true_force_N'] - points[0]['true_force_N']) / delta_height
            Jnew = Jold + np.outer(v, v) / v[site]
            Jnewclean = Joldclean + np.outer(vclean, vclean) / vclean[site]
            Hhat = N.T @ Jnew @ N
            Hhat = (Hhat + Hhat.T) / 2
            Htrue = np.linalg.inv(N.T @ C @ N)
            Jtrue = N @ Htrue @ N.T
            oldrho = eps / s * np.sqrt(m * (m - 3))
            dv = 2 * eps / delta_height * np.sqrt(n)
            sig = 2 * eps / delta_height
            nv = float(np.linalg.norm(v))
            den = float(v[site] - sig)
            rho = float(oldrho + (2 * nv + dv) * dv / den + nv * nv * sig / (v[site] * den)) if den > 0 else None
            eig = np.linalg.eigvalsh(Hhat)
            lower = float(eig.min() - rho) if rho else None
            upper = float(eig.max() + rho) if rho else None
            base_true = points[1]['true_force_N']
            base_obs = base_true + N @ (N.T @ points[1]['noise_N'])
            hbase = points[1]['h']
            newbranches = all((p['true_force_N'].min() > 0 for p in points))
            acquisition = oldbranches and newbranches and (den > 0) and (lower is not None) and (lower > 0)
            Ghat = np.linalg.inv(Hhat)
            Q = np.eye(n) - N @ N.T
            Chat = N @ Ghat @ N.T + Q / 750
            ghat = -Chat @ base_obs
            taskrows = []
            for (task, h) in queries(n):
                (fr, qr, oracle) = solve(A, C, g, w, h, f0)
                (fa, qa, meta) = solve(A, Chat, ghat, w, h - hbase, base_obs)
                (fc, control) = dual_contact(N, Hhat, base_obs, h - hbase)
                z = N.T @ (fa - base_obs)
                rad = float(np.sqrt(upper / lower) * (rho * np.linalg.norm(Ghat @ z) + eps * np.sqrt(n)) + 1e-06) if acquisition else None
                err = float(np.max(np.abs(fa - fr)))
                l2 = float(np.linalg.norm(fa - fr))
                taskrows.append({'query': task, 'h_mm': h.tolist(), 'oracle_force_N': fr.tolist(), 'candidate_force_N': fa.tolist(), 'radius_N': rad, 'max_force_error_N': err, 'l2_force_error_N': l2, 'radius_contains_truth': bool(rad is not None and l2 <= rad + 1e-06), 'radius_useful': bool(rad is not None and rad <= 1), 'numerical_gate': bool(meta['balance_residual_N'] <= 1e-07 and meta['min_force_N'] >= -1e-07 and (meta['active_stationarity_mm'] <= 1e-07) and (meta['min_slack_mm'] >= -1e-07)), 'control_parity_N': float(np.max(np.abs(fa - fc))), 'KKT': meta})
            worlds.append({'hidden_k_N_per_mm': ki, 'source_information_old_force_vectors_N': all_oldforces, 'old_probe_branches_valid': bool(oldbranches), 'small_invisible_probes_force_N': small, 'activation': [{'height_mm': p['height_mm'], 'true_force_N': p['true_force_N'].tolist(), 'observed_force_N': p['observed_force_N'].tolist(), 'KKT': p['KKT']} for p in points], 'new_probe_branches_valid': bool(newbranches), 'new_acquisition_gate': bool(acquisition), 'J_old_noiseless': Joldclean.tolist(), 'v_N_per_mm': v.tolist(), 'H_observed_N_per_mm': Hhat.tolist(), 'new_noiseless_tangent_relative_error': float(np.max(np.abs(Jnewclean - Jtrue)) / np.max(np.abs(Jtrue))), 'old_radius_N_per_mm': float(oldrho), 'new_radius_N_per_mm': rho, 'actual_new_spectral_error_N_per_mm': float(np.linalg.norm(Hhat - Htrue, 2)), 'spectral_lower_N_per_mm': lower, 'spectral_upper_N_per_mm': upper, 'baseline_force_N': base_obs.tolist(), 'baseline_height_mm': hbase.tolist(), 'queries': taskrows})
        old_parity = float(np.max(np.abs(np.array(worlds[0]['source_information_old_force_vectors_N']) - np.array(worlds[1]['source_information_old_force_vectors_N']))))
        hidden_diff = float(np.max(np.abs(np.array(worlds[0]['activation'][0]['true_force_N']) - np.array(worlds[1]['activation'][0]['true_force_N']))))
        rows.append({'case': case, 'source': src, 'n': n, 'invisible_site': site, 'invisible_fdi': ids[site], 'A': A.tolist(), 'w_N': w.tolist(), 'f0_N': f0.tolist(), 'old_information_parity_N': old_parity, 'hidden_50um_response_difference_N': hidden_diff, 'worlds': worlds})
    worlds = [v for r in rows for v in r['worlds']]
    tasks = [t for v in worlds for t in v['queries']]
    gate = all((v['new_acquisition_gate'] and v['new_noiseless_tangent_relative_error'] <= 1e-08 and (v['actual_new_spectral_error_N_per_mm'] <= v['new_radius_N_per_mm'] + 1e-07) for v in worlds)) and all((t['radius_useful'] and t['radius_contains_truth'] and t['numerical_gate'] for t in tasks))
    out = {'round': 'R5', 'claim_type': 'information_link', 'realization': 'SIMULATED_PROSPECTIVE_ACQUISITION', 'prereg_sha256': sha(ROOT / 'PREREG_R5.json'), 'n_cases': len(rows), 'n_worlds': len(worlds), 'n_queries': len(tasks), 'old_information_parity_gate': all((r['old_information_parity_N'] <= 1e-07 for r in rows)), 'max_hidden_response_difference_N': max((r['hidden_50um_response_difference_N'] for r in rows)), 'n_cases_hidden_response_difference_gt_1_N': sum((r['hidden_50um_response_difference_N'] > 1 for r in rows)), 'exact_rank_completion_gate': all((v['new_noiseless_tangent_relative_error'] <= 1e-08 for v in worlds)), 'acquisition_gate': all((v['new_acquisition_gate'] for v in worlds)), 'recovery_gate': 'PASS' if gate else 'FAIL', 'useful_radius_queries': sum((t['radius_useful'] for t in tasks)), 'truth_containment_gate': all((t['radius_contains_truth'] for t in tasks)), 'max_radius_N': max((t['radius_N'] for t in tasks if t['radius_N'] is not None)), 'max_force_error_N': max((t['max_force_error_N'] for t in tasks)), 'max_control_parity_N': max((t['control_parity_N'] for t in tasks)), 'matched_control_gate': 'PASS' if all((t['control_parity_N'] <= 0.001 for t in tasks)) else 'FAIL', 'rows': rows, 'physical_measurements': 'NOT_RUN', 'empirical_validity': 'UNKNOWN', 'wall_s': time.perf_counter() - start, 'cpu_s': time.process_time() - cpu, 'peak_rss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024}
    write(ROOT / 'raw/R5.json', out)
    print(json.dumps({k: v for (k, v) in out.items() if k != 'rows'}, indent=2))
    write(ROOT / 'CURRENT_WORK_STATE.json', {'lane': _release_expand('XBREAK_HUNT_1'), 'status': 'R5_COMPLETE', 'latest_gate': out['recovery_gate'], 'milestone': 'Initial inactive information obstruction and directed activation completion executed', 'next_operation': 'Export runnable edit-response consumer, verify metrology requirements, then hand off physical calibration', 'evidence': 'raw/R5.json', 'empirical_validity': 'UNKNOWN'})
if __name__ == '__main__':
    main()
