from dental_release.paths import expand as _release_expand
from bench import *
from force_operator import balanced_basis
from active_contact import solve
import time, resource

def directional_radii(N, H, f, fbase, v, rho_old, global_rho, corner_e, corner_v):
    n = len(f)
    G = np.linalg.inv(H)
    z = N.T @ (f - fbase)
    x = G @ z
    y = N @ x
    sigma = 0.5
    delta = sigma * np.sqrt(n)
    si = v[-1]
    nv = np.linalg.norm(v)
    a = float(v @ y)
    corr = y.copy()
    corr[-1] -= a / si
    rank = (delta * (abs(a) + sigma * np.linalg.norm(y, 1)) + nv * sigma * np.linalg.norm(corr, 1)) / (si - sigma)
    eig = np.linalg.eigvalsh(H)
    lower = eig.min() - global_rho
    upper = eig.max() + global_rho
    factor = np.sqrt(upper / lower)
    common = rho_old * np.linalg.norm(x) + 0.01 * np.sqrt(n)
    rc = float(factor * (common + rank) + 1e-06)
    tc = time.perf_counter()
    dots = corner_e @ corr
    norms2 = (a * a * n * sigma * sigma + 2 * a * dots * corner_v + dots * dots * nv * nv) / (si * si)
    linear = float(np.sqrt(max(float(norms2.max()), 0)))
    remainder = sigma / (si * (si - sigma)) * (abs(a) * delta + nv * sigma * np.linalg.norm(corr, 1)) + delta * sigma * np.linalg.norm(y, 1) / (si - sigma)
    rs = float(factor * (common + linear + remainder) + 1e-06)
    return (rc, rs, {'linear_corner_norm_N': linear, 'remainder_N': float(remainder), 'rank_triangle_action_N': float(rank), 'query_compliance_norm_mm': float(np.linalg.norm(x)), 'corner_control_s': time.perf_counter() - tc, 'vertices': len(corner_e)})

def main():
    check_frozen('R6')
    start = time.perf_counter()
    cpu = time.process_time()
    parent = json.loads((ROOT / 'raw/R5.json').read_text())
    rows = []
    control_s = 0
    vertices = 0
    for source in parent['rows']:
        case = source['case']
        (ids, A, M, src) = geometry(case)
        n = len(ids)
        N = balanced_basis(A)
        site = n - 1
        xy = np.array(src['xy_mm'])
        dist2 = np.sum((xy[:, None] - xy[None, :]) ** 2, axis=2)
        k = 750 * (1 + 0.2 * np.cos(2 * np.pi * np.arange(n) / n))
        f0 = np.array(source['f0_N'])
        w = A.T @ f0
        ce = ((np.arange(2 ** n, dtype=np.uint32)[:, None] >> np.arange(n) & 1).astype(float) * 2 - 1) * 0.5
        for world in source['worlds']:
            ki = world['hidden_k_N_per_mm']
            kk = k.copy()
            kk[site] = ki
            C = np.diag(1 / kk) + 0.8 / 750 * np.exp(-dist2 / (2 * 15 ** 2))
            g = -C @ f0
            g[site] += 0.04
            H = np.array(world['H_observed_N_per_mm'])
            v = np.array(world['v_N_per_mm'])
            fbase = np.array(world['baseline_force_N'])
            hbase = np.array(world['baseline_height_mm'])
            G = np.linalg.inv(H)
            Chat = N @ G @ N.T + (np.eye(n) - N @ N.T) / 750
            ghat = -Chat @ fbase
            cv = ce @ v
            tasks = []
            for i in range(n):
                for height in [-0.04, -0.03, 0.03, 0.04]:
                    h = np.eye(n)[i] * height
                    (ft, qt, oracle) = solve(A, C, g, w, h, f0)
                    (fa, qa, meta) = solve(A, Chat, ghat, w, h - hbase, fbase)
                    (rc, rs, details) = directional_radii(N, H, fa, fbase, v, world['old_radius_N_per_mm'], world['new_radius_N_per_mm'], ce, cv)
                    control_s += details['corner_control_s']
                    vertices += details['vertices']
                    err = float(np.linalg.norm(fa - ft))
                    numeric = meta['balance_residual_N'] <= 1e-07 and oracle['active_stationarity_mm'] <= 1e-07
                    tasks.append({'site': i, 'fdi': ids[i], 'height_mm': height, 'oracle_force_N': ft.tolist(), 'candidate_force_N': fa.tolist(), 'l2_error_N': err, 'candidate_radius_N': rc, 'strong_control_radius_N': rs, 'candidate_truth_enclosed': bool(err <= rc + 1e-06), 'control_truth_enclosed': bool(err <= rs + 1e-06), 'candidate_useful': bool(rc <= 1), 'control_useful': bool(rs <= 1), 'numerical_gate': bool(numeric), 'directional_details': details})
            rows.append({'case': case, 'hidden_k_N_per_mm': ki, 'n': n, 'queries': tasks})
    tasks = [q for r in rows for q in r['queries']]
    cg = all((t['candidate_truth_enclosed'] and t['candidate_useful'] and t['numerical_gate'] for t in tasks))
    sg = all((t['control_truth_enclosed'] and t['control_useful'] and t['numerical_gate'] for t in tasks))
    out = {'round': 'R6', 'claim_type': 'algorithm', 'prereg_sha256': sha(ROOT / 'PREREG_R6.json'), 'source_R5_sha256': sha(ROOT / 'raw/R5.json'), 'n_worlds': len(rows), 'n_new_queries': len(tasks), 'candidate_gate': 'PASS' if cg else 'FAIL', 'strong_control_gate': 'PASS' if sg else 'FAIL', 'candidate_useful_queries': sum((t['candidate_useful'] for t in tasks)), 'control_useful_queries': sum((t['control_useful'] for t in tasks)), 'candidate_truth_containment_gate': all((t['candidate_truth_enclosed'] for t in tasks)), 'control_truth_containment_gate': all((t['control_truth_enclosed'] for t in tasks)), 'max_candidate_radius_N': max((t['candidate_radius_N'] for t in tasks)), 'max_control_radius_N': max((t['strong_control_radius_N'] for t in tasks)), 'max_actual_l2_error_N': max((t['l2_error_N'] for t in tasks)), 'corner_vertices_evaluated': vertices, 'corner_control_s': control_s, 'rows': rows, 'wall_s': time.perf_counter() - start, 'cpu_s': time.process_time() - cpu, 'peak_rss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'physical_measurements': 'NOT_RUN', 'empirical_validity': 'UNKNOWN'}
    write(ROOT / 'raw/R6.json', out)
    print(json.dumps({k: v for (k, v) in out.items() if k != 'rows'}, indent=2))
    write(ROOT / 'CURRENT_WORK_STATE.json', {'lane': _release_expand('XBREAK_HUNT_1'), 'status': 'R6_COMPLETE', 'latest_gate': out['candidate_gate'], 'milestone': 'Joint rank-one directional uncertainty and exact corner control on new heights executed', 'next_operation': 'Package measured-response consumer and hand off physical acquisition requirements with failed radius gates', 'evidence': 'raw/R6.json', 'empirical_validity': 'UNKNOWN'})
if __name__ == '__main__':
    main()
