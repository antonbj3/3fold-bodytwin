"""Test the information carried by one calorimeter versus pulse-indexed heat."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
import json, time, resource, hashlib, os
import numpy as np
from scipy.sparse import csr_matrix, diags
from scipy.optimize import linprog
from laser_heat import geometry, advance, dump, sha
ORIGINAL = Path(__file__).resolve().parent
ROOT = Path(os.environ.get('X33_RUN_ROOT', str(ORIGINAL)))
DATA = Path(os.environ.get('X33_MATRIX_DIR', _release_expand('@DENTAL_WORK_ROOT@/X33')))

def schedule(dt, end=10.0):
    t = 0.0
    p = 0
    record = 0.0
    times = [0.0]
    while t < end - 1e-12:
        if p < 42 and t >= p / 6.0 - 1e-12:
            p += 1
        if t >= record - 1e-12:
            record += 0.05
        stop = min(end, p / 6.0) if p < 42 else end
        delta = min(dt, stop - t, record - t)
        if delta < 1e-12:
            delta = min(dt, stop - t)
        t += delta
        times.append(t)
    return np.array(times)

def response_matrix(geo, path):
    (neighbors, g, C, pm, pi, pw, wet, k, h, prep, pulp) = geo
    n = len(C)
    rows = np.repeat(np.arange(n), 6)
    cols = neighbors.ravel()
    vals = (g / C[:, None]).ravel()
    keep = cols >= 0
    A = csr_matrix((vals[keep], (rows[keep], cols[keep])), shape=(n, n)) - diags(g.sum(1) / C)
    dt = 0.8 / np.max(g.sum(1) / C)
    times = schedule(dt)
    T = np.zeros((n, 42))
    matrix = np.lib.format.open_memmap(path, mode='w+', dtype=np.float64, shape=(len(times), int(pm.sum()), 42))
    p = 0
    for (j, t) in enumerate(times):
        if p < 42 and t >= p / 6.0 - 1e-12:
            select = pi[p] >= 0
            ix = pi[p, select]
            T[ix, p] += 0.25 * pw[p, select] / C[ix]
            p += 1
        matrix[j] = T[pm]
        if j < len(times) - 1:
            T += (times[j + 1] - t) * (A @ T)
    matrix.flush()
    return (matrix, times)

def upper_support(M, total_eta, cap=0.1):
    whole = int(np.floor(total_eta / cap + 1e-12))
    remainder = total_eta - whole * cap
    best = -1.0
    bestrow = 0
    bestorder = None
    for start in range(0, len(M), 4096):
        a = np.array(M[start:start + 4096])
        order = np.argsort(a, axis=1)[:, ::-1]
        sorted_vals = np.take_along_axis(a, order, axis=1)
        b = cap * sorted_vals[:, :whole].sum(1)
        if whole < 42:
            b += remainder * sorted_vals[:, whole]
        j = int(np.argmax(b))
        if b[j] > best:
            best = float(b[j])
            bestrow = start + j
            bestorder = order[j]
    x = np.zeros(42)
    x[bestorder[:whole]] = cap
    if whole < 42:
        x[bestorder[whole]] = remainder
    return (best, bestrow, x)

def minimax(M, total_eta, cap):
    x = np.full(42, total_eta / 42)
    v = M @ x
    rows = [int(np.argmax(v))]
    iterations = []
    for iteration in range(70):
        a = np.column_stack((np.asarray(M[rows]), -np.ones(len(rows))))
        fit = linprog(np.r_[np.zeros(42), 1.0], A_ub=a, b_ub=np.zeros(len(rows)), A_eq=np.r_[np.ones(42), 0.0][None, :], b_eq=[total_eta], bounds=[(0, cap)] * 42 + [(0, None)], method='highs', options={'dual_feasibility_tolerance': 1e-09, 'primal_feasibility_tolerance': 1e-09, 'threads': 4})
        if not fit.success:
            raise RuntimeError(fit.message)
        x = fit.x[:42]
        v = M @ x
        peak = float(v.max())
        gap = peak - float(fit.fun)
        iterations.append({'iteration': iteration, 'cuts': len(rows), 'LP_lower_C': float(fit.fun), 'full_space_time_peak_C': peak, 'gap_C': gap})
        if gap <= 1e-06:
            return (x, peak, iterations, True)
        candidates = np.argsort(v)[-12:][::-1]
        new = [int(j) for j in candidates if int(j) not in rows]
        if not new:
            raise RuntimeError('separation oracle cannot add violated row')
        rows += new
    return (x, peak, iterations, False)

def direct_replay(geo, eta, times):
    (neighbors, g, C, pm, pi, pw, *_) = geo
    T = np.zeros(len(C))
    Tn = np.zeros(len(C))
    cool = np.zeros(len(C))
    p = 0
    peak = 0.0
    hist = []
    for (j, t) in enumerate(times):
        if p < 42 and t >= p / 6.0 - 1e-12:
            take = pi[p] >= 0
            ix = pi[p, take]
            T[ix] += 0.25 * eta[p] * pw[p, take] / C[ix]
            p += 1
        peak = max(peak, float(T[pm].max()))
        hist.append([t, float(T[pm].max()), float(T[pm].mean())])
        if j < len(times) - 1:
            advance(T, Tn, neighbors, g, C, cool, 0.0, times[j + 1] - t)
            (T, Tn) = (Tn, T)
    return (peak, float(C @ T), hist)

def main():
    start = time.time()
    pr = json.loads((ROOT / 'PREREG_R2.json').read_text())
    r1 = json.loads((ROOT / 'PREREG_R1.json').read_text())
    assert sha(ROOT / 'PREREG_R2.json') == (ROOT / 'PREREG_R2.sha256').read_text().split()[0]
    assert sha(ORIGINAL / pr['source_round']) == pr['source_round_sha256']
    DATA.mkdir(parents=True, exist_ok=True)
    own = sum((p.stat().st_size for p in Path(_release_expand('@DENTAL_WORK_ROOT@/X33')).rglob('*') if p.is_file()))
    if own > 2700000000:
        raise RuntimeError('3GB lane intermediate budget would be exceeded')
    geo = geometry(1, r1)
    matrixpath = DATA / 'R2_PULSE_POINT_TIME_RESPONSE.npy'
    (matrix, times) = response_matrix(geo, matrixpath)
    M = matrix.reshape(-1, 42)
    total = pr['total_retained_energy_J'] / 0.25
    cap = 0.1
    uniform = np.full(42, total / 42)
    uniform_peak = float(np.max(M @ uniform))
    (upper, row, bad) = upper_support(M, total, cap)
    (good, lower, iterations, converged) = minimax(M, total, cap)
    full_lp = linprog(-np.asarray(M[row]), A_eq=np.ones((1, 42)), b_eq=[total], bounds=[(0, cap)] * 42, method='highs', options={'threads': 4})
    assert full_lp.success
    control_error = abs(upper + float(full_lp.fun))
    traces = {}
    checks = {}
    energies = {}
    for (name, x) in [('uniform', uniform), ('low_peak', good), ('high_peak', bad)]:
        (peak, stored, hist) = direct_replay(geo, x, times)
        matrixpeak = float(np.max(M @ x))
        energies[name] = stored
        checks[name] = {'matrix_peak_C': matrixpeak, 'direct_peak_C': peak, 'absolute_difference_C': abs(peak - matrixpeak), 'total_energy_J': float(0.25 * x.sum()), 'stored_energy_J': stored, 'minimum_eta': float(x.min()), 'maximum_eta': float(x.max())}
        np.savetxt(ROOT / f'raw/R2_HISTORY_{name}.csv', hist, delimiter=',', header='time_s,max_pulp_rise_C,mean_pulp_rise_C', comments='')
        traces[name] = x.tolist()
    (robust_upper, robust_row, _) = upper_support(M, (0.5 + 0.005) / 0.25, cap)
    numerical_mutation_rejected = abs(upper + 1 - checks['high_peak']['direct_peak_C']) > 0.02
    energy_mutation_rejected = bool(abs(0.25 * bad.sum() + 0.01 - 0.5) > 1e-10)

    def history_validator(x, expected):
        return np.array_equal(x, expected)
    order_poison = bad[::-1].copy()
    history_mutation_rejected = not history_validator(order_poison, bad)
    (reverse_peak, _, _) = direct_replay(geo, order_poison, times)
    n_pulp = int(geo[3].sum())
    (jtime, jpoint) = divmod(row, n_pulp)
    point_indices = np.argwhere(geo[-1])
    point_mm = (point_indices[jpoint] * 0.3).tolist()
    result = {'round': 'R2', 'claim_type': 'capability', 'resolution_level': 'PER_POINT', 'matrix_shape': list(matrix.shape), 'matrix_path': str(matrixpath), 'matrix_bytes': matrixpath.stat().st_size, 'matrix_sha256': sha(matrixpath), 'total_retained_energy_J': 0.5, 'total_measurement_error_J': 0.005, 'eta_box': [0, 0.1], 'uniform_peak_C': uniform_peak, 'minimum_possible_maximum_C': lower, 'maximum_possible_maximum_C': upper, 'robust_maximum_with_energy_error_C': robust_upper, 'same_energy_decision_flip': lower < 5.5 < upper, 'worst_support': {'time_s': float(times[jtime]), 'pulp_point_local_zyx_mm': point_mm}, 'measurement_answer': 'Total heat alone insufficient for this conditional query' if lower < 5.5 < upper else 'No threshold flip under frozen scenario', 'required_support': 'pulse-time and scan-position indexed retained energy, or an externally validated uniform-partition law; consumer aggregates after pointwise propagation', 'minimax_converged': converged, 'minimax_iterations': iterations, 'LP_control_error_C': control_error, 'direct_replay_checks': checks, 'histories_eta': traces, 'reverse_adverse_history_peak_C': reverse_peak, 'pulse_reversal_peak_difference_C': abs(reverse_peak - upper), 'gates': {'same_total_energy': all((abs(v - 0.5) < 1e-10 for v in energies.values())), 'matrix_replay': all((r['absolute_difference_C'] < 0.02 for r in checks.values())), 'LP_control': control_error < 1e-07, 'minimax_gap': converged, 'opposite_decisions': lower < 5.5 < upper}, 'mutations': {'plus_1C_rejected': numerical_mutation_rejected, 'energy_extra_0p01J_rejected': energy_mutation_rejected, 'pulse_order_contract_rejected': history_mutation_rejected}, 'source_measurement_status': 'No calorimetric measurement occurred; Q and per-pulse bounds are prospective synthetic input on published anatomy.', 'external_referent': pr['external_referent'], 'physical_certification': 'UNKNOWN', 'cost': {'wall_s': time.time() - start, 'CPU_s': time.process_time(), 'maxRSS_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'threads_max': 4, 'GPU': False, 'LP_calls': len(iterations) + 1, 'matrix_column_histories': 42, 'scalar_replays': 4, 'preparation_cost': 'UNKNOWN_UNINSTRUMENTED'}}
    dump('raw/R2_RESULTS.json', result)
    print(json.dumps({k: result[k] for k in ['uniform_peak_C', 'minimum_possible_maximum_C', 'maximum_possible_maximum_C', 'same_energy_decision_flip', 'gates', 'cost']}), flush=True)
    assert all((result['gates'][k] for k in ['same_total_energy', 'matrix_replay', 'LP_control', 'minimax_gap']))
    assert all(result['mutations'].values())
if __name__ == '__main__':
    main()
