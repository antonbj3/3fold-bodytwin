from dental_release.paths import expand as _release_expand
from bench import *
from scipy.optimize import linprog
import time, resource

def lp(obj, M, b, bounds, Aub=None, bub=None):
    r = linprog(obj, A_eq=M, b_eq=b, A_ub=Aub, b_ub=bub, bounds=bounds, method='highs', options={'threads': 4})
    assert r.success, r.message
    return r

def product_response_bound(M, b, bounds, i, j, h, t):
    n = len(bounds)
    aug = np.column_stack((M, np.zeros(len(M))))
    (li, ui) = bounds[i]
    (lj, uj) = bounds[j]
    rows = []
    rhs = []

    def inequality(ci, cj, cz, r):
        v = np.zeros(n + 1)
        v[i] += ci
        v[j] += cj
        v[-1] = cz
        rows.append(v)
        rhs.append(r)
    if i == j:
        for x in np.linspace(li, ui, 17):
            inequality(2 * x, 0, -1, x * x)
        inequality(-(li + ui), 0, 1, -li * ui)
    else:
        inequality(lj, li, -1, li * lj)
        inequality(uj, ui, -1, ui * uj)
        inequality(-lj, -ui, 1, -ui * lj)
        inequality(-uj, -li, 1, -li * uj)
    obj = np.zeros(n + 1)
    obj[j] = h * (i == j)
    obj[-1] = -h * t
    bz = [*bounds, (li * lj, ui * uj)]
    lo = lp(obj, aug, b, bz, np.array(rows), np.array(rhs))
    hi = lp(-obj, aug, b, bz, np.array(rows), np.array(rhs))
    return [float(lo.fun) - 1e-09, float(-hi.fun) + 1e-09]

def box_bound(bounds, i, j, h, t):
    (l, u) = bounds[i]
    (ll, uu) = bounds[j]
    if i == j:
        vals = [h * (x - t * x * x) for x in (l, u)]
        if t > 0 and l <= 1 / (2 * t) <= u:
            vals.append(h / (4 * t))
    else:
        vals = [-h * t * l * ll, -h * t * u * uu]
    return [float(min(vals)) - 1e-09, float(max(vals)) + 1e-09]

def main():
    check_frozen('R2')
    start = time.perf_counter()
    cpu = time.process_time()
    panel = json.loads((ROOT / 'PREREG_R1.json').read_text())['source_panel']
    rows = []
    eps = 0.001
    h = 0.01
    for case in panel:
        (ids, A, M, src) = geometry(case)
        n = len(ids)
        truth = 750 * (1 + 0.2 * np.cos(2 * np.pi * np.arange(n) / n))
        (w, q0) = base_w(A, truth)
        K = A.T @ (truth[:, None] * A)
        b = M @ truth
        Ki = np.linalg.inv(K)
        probe = np.array([-Ki @ a * h for a in A])
        y = probe * truth[:, None]
        exact = np.sum(probe * y, axis=1) / np.sum(probe * probe, axis=1)
        rawbounds = []
        for site in range(n):
            (lo, hi) = (500.0, 1130.0)
            for j in range(3):
                if abs(probe[site, j]) > 1e-15:
                    ends = sorted([(y[site, j] - eps) / probe[site, j], (y[site, j] + eps) / probe[site, j]])
                    lo = max(lo, ends[0])
                    hi = min(hi, ends[1])
            rawbounds.append((lo, hi))
        bounds = []
        for i in range(n):
            obj = np.eye(n)[i]
            lo = lp(obj, M, b, rawbounds)
            hi = lp(-obj, M, b, rawbounds)
            bounds.append((float(lo.fun), float(-hi.fun)))
        control = np.empty((n, n, 2))
        cand = np.empty_like(control)
        true = np.empty((n, n))
        for i in range(n):
            hh = np.eye(n)[i] * h
            (q, f, _) = fixed(A, truth, w, hh)
            (_, f0, _) = fixed(A, truth, w, np.zeros(n))
            true[:, i] = f - f0
            for j in range(n):
                t = float(A[j] @ Ki @ A[i])
                cand[j, i] = box_bound(rawbounds, i, j, h, t)
                control[j, i] = product_response_bound(M, b, bounds, i, j, h, t)
        enclose_c = np.all((true >= cand[:, :, 0] - 1e-07) & (true <= cand[:, :, 1] + 1e-07))
        enclose_l = np.all((true >= control[:, :, 0] - 1e-07) & (true <= control[:, :, 1] + 1e-07))
        width_c = float(np.max(cand[:, :, 1] - cand[:, :, 0]))
        width_l = float(np.max(control[:, :, 1] - control[:, :, 0]))
        rows.append({'case': case, 'n': n, 'source': src, 'truth_k_N_per_mm': truth.tolist(), 'K': K.tolist(), 'pose_probe_delta_mm': y.tolist(), 'pose_noise_bound_mm': eps, 'max_probe_component_um': float(np.max(np.abs(y)) * 1000), 'noiseless_relative_support_error': float(np.max(np.abs(exact - truth) / truth)), 'individual_probe_bounds': rawbounds, 'joint_LP_bounds': bounds, 'truth_held_delta_f_N': true.tolist(), 'candidate_intervals_N': cand.tolist(), 'control_intervals_N': control.tolist(), 'max_candidate_width_N': width_c, 'max_control_width_N': width_l, 'truth_enclosed': bool(enclose_c and enclose_l), 'candidate_information_gate': 'PASS' if enclose_c and width_c <= 2 else 'FAIL', 'control_information_gate': 'PASS' if enclose_l and width_l <= 2 else 'FAIL'})
    out = {'round': 'R2', 'claim_type': 'information_link', 'prereg_sha256': sha(ROOT / 'PREREG_R2.json'), 'n_cases': len(rows), 'noiseless_gate': all((r['noiseless_relative_support_error'] <= 1e-08 for r in rows)), 'truth_containment_gate': all((r['truth_enclosed'] for r in rows)), 'candidate_information_gate': 'PASS' if all((r['candidate_information_gate'] == 'PASS' for r in rows)) else 'FAIL', 'strong_control_information_gate': 'PASS' if all((r['control_information_gate'] == 'PASS' for r in rows)) else 'FAIL', 'max_candidate_width_N': max((r['max_candidate_width_N'] for r in rows)), 'max_control_width_N': max((r['max_control_width_N'] for r in rows)), 'rows': rows, 'wall_s': time.perf_counter() - start, 'cpu_s': time.process_time() - cpu, 'peak_rss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'external_matched_facit': 'ABSENT', 'physical_measurement': 'NOT_RUN'}
    write(ROOT / 'raw/R2.json', out)
    print(json.dumps({k: v for (k, v) in out.items() if k != 'rows'}, indent=2))
    write(ROOT / 'CURRENT_WORK_STATE.json', {'lane': _release_expand('XBREAK_HUNT_1'), 'status': 'R2_COMPLETE', 'latest_gate': out['candidate_information_gate'], 'milestone': 'Local pose reveals support exactly, but declared noise sets force interval width', 'next_operation': 'Change measurement from global displacement to regional force response; avoid recovering tissue support', 'evidence': 'raw/R2.json', 'physical_support': 'UNKNOWN'})
if __name__ == '__main__':
    main()
