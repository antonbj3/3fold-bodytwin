"""Finite ramp acquisition. Hidden oracle supports never enter work_bounds()."""
import json, resource, time
from pathlib import Path
import numpy as np
from run_r1 import ROOT, Oracle, model, sha, write

def sample_ramp(oracle, g, T, n):
    t = time.perf_counter()
    ts = np.linspace(0, T, n + 1)
    ds = [float(np.min(g))]
    active = None
    coeff = None
    slow = 0
    fast = 0
    maxres = 0.0
    for load in ts[1:]:
        accepted = False
        if coeff is not None:
            pol = coeff[:, 0] + load * coeff[:, 1]
            f = np.zeros(len(g))
            f[active] = pol[:-1]
            d = pol[-1]
            grad = oracle.H @ f + g
            residual = max(float(np.max(np.abs(grad[active] - d))), float(max(0.0, d - grad.min())), float(max(0.0, -f.min())))
            if residual <= 1e-07 and abs(f.sum() - load) <= 1e-07:
                ds.append(float(d))
                fast += 1
                accepted = True
                maxres = max(maxres, residual)
        if not accepted:
            ans = oracle.solve(float(load), g)
            slow += 1
            ds.append(ans['d'])
            active = np.asarray(ans['active'], dtype=int)
            K = np.block([[oracle.H[np.ix_(active, active)], -np.ones((len(active), 1))], [np.ones((1, len(active))), np.zeros((1, 1))]])
            rhs = np.zeros((len(active) + 1, 2))
            rhs[:-1, 0] = -g[active]
            rhs[-1, 1] = 1
            coeff = np.linalg.lstsq(K, rhs, rcond=1e-13)[0]
    ds = np.array(ds)
    if np.min(np.diff(ds)) < -2e-07:
        raise ArithmeticError('Closure ramp is not monotone')
    return (ts, ds, {'readings': n + 1, 'slow_QP': slow, 'fast_KKT': fast, 'max_fast_residual': maxres, 'seconds': time.perf_counter() - t})

def work_bounds(ts, observed, eta, numerical_guard):
    dt = np.diff(ts)
    err = eta + numerical_guard
    return (float(dt @ (observed[:-1] - err)), float(dt @ (observed[1:] + err)))

def force_bounds(minus, zero, plus, h, T):
    return (max(0.0, (plus[0] - zero[1]) / h), min(T, (zero[1] - minus[0]) / h))

def main():
    start = time.perf_counter()
    cpu = time.process_time()
    p = ROOT / 'PREREG_R2.json'
    cfg = json.loads(p.read_text())
    frozen = sha(p)
    write('raw/PREREG_R2_FROZEN.json', {'sha256': frozen, 'file': str(p), 'physical_creation_before_run': True})
    worlds = json.loads((ROOT / 'raw/R1_WORLDS.json').read_text())
    records = []
    raw = []
    cost = []
    coverage_work = []
    maxn = max(cfg['ramp_intervals'])
    T = cfg['T_N']
    guard = cfg['numerical_displacement_guard_mm']
    for world in worlds:
        (oracle, teeth, _, _) = model(world['case'], world['scenario'])
        assert np.array_equal(oracle.g, world['g_mm']) and np.array_equal(oracle.c, world['c_mm_per_N'])
        upper = [t for t in teeth if t < 30]
        targets = sorted(upper, key=lambda t: (-world['native_force_N'][str(t)], t))[:2]
        ramps = {}
        (ts, ds, stats) = sample_ramp(oracle, oracle.g, T, maxn)
        cost.append(stats)
        ramps['base'] = (ts, ds, world['native_V_Nmm'])
        raw.append({'case': world['case'], 'scenario': world['scenario'], 'ramp': 'base', 'force_N': ts.tolist(), 'd_mm': ds.tolist()})
        for tooth in targets:
            direction = oracle.B[teeth.index(tooth)]
            for h in cfg['h_mm']:
                for sign in [-1, 1]:
                    g = oracle.g + sign * h * direction
                    (tt, dd, stats) = sample_ramp(oracle, g, T, maxn)
                    cost.append(stats)
                    ideal = oracle.solve(T, g)['V']
                    ramps[tooth, h, sign] = (tt, dd, ideal)
                    raw.append({'case': world['case'], 'scenario': world['scenario'], 'ramp': f'{tooth}_{h}_{sign}', 'force_N': tt.tolist(), 'd_mm': dd.tolist()})
                for n in cfg['ramp_intervals']:
                    ix = np.arange(0, maxn + 1, maxn // n)
                    for eta in cfg['displacement_error_mm']:
                        bounds = []
                        for (j, key) in enumerate([(tooth, h, -1), 'base', (tooth, h, 1)]):
                            (tt, dd, ideal) = ramps[key]
                            observed = dd[ix] + eta * (0.6 * np.sin(np.arange(len(ix)) * 0.73 + j) + 0.4 * np.cos(j))
                            wb = work_bounds(tt[ix], observed, eta, guard)
                            bounds.append(wb)
                            coverage_work.append(wb[0] - 1e-05 <= ideal <= wb[1] + 1e-05)
                        (lo, hi) = force_bounds(*bounds, h, T)
                        q = world['native_force_N'][str(tooth)]
                        records.append({'case': world['case'], 'scenario': world['scenario'], 'FDI': tooth, 'h_mm': h, 'intervals': n, 'error_mm': eta, 'facit_N': q, 'lower_N': lo, 'upper_N': hi, 'width_pp': 100 * (hi - lo) / T, 'covered': lo - 1e-05 <= q <= hi + 1e-05, 'work_bounds_Nmm': bounds})
        write('raw/R2_PROGRESS.json', {'worlds_completed': len(cost) // 13, 'records': records, 'cost': cost})
        print(world['case'], world['scenario'], 'ramps completed', flush=True)
    summary = []
    for eta in cfg['displacement_error_mm']:
        for n in cfg['ramp_intervals']:
            for h in cfg['h_mm']:
                rr = [r for r in records if r['error_mm'] == eta and r['intervals'] == n and (r['h_mm'] == h)]
                widths = [r['width_pp'] for r in rr]
                summary.append({'error_mm': eta, 'intervals': n, 'h_mm': h, 'n': len(rr), 'covered': sum((r['covered'] for r in rr)), 'fraction_le5pp': float(np.mean(np.array(widths) <= 5)), 'median_width_pp': float(np.median(widths)), 'max_width_pp': max(widths)})

    def gate(eta):
        return any((r['fraction_le5pp'] >= 0.8 for r in summary if r['error_mm'] == eta and r['intervals'] == 512))
    result = {'round': 'R2', 'claim_type': 'information_link', 'prereg_sha256': frozen, 'status': 'PENDING_INDEPENDENT_REVIEW', 'coverage_gate': all((r['covered'] for r in records)), 'work_recovery_gate': all(coverage_work), 'practical_gate_1um': gate(0.001), 'ideal_metrology_gate': gate(0), 'numerical_guard_status': 'Declared conservative 1e-7 mm guard; not a formal roundoff certificate', 'strong_equal_control': 'Same monotone work enclosure and secants; no algorithm advantage', 'cost': {'ramp_count': len(cost), 'maximum_force_readings': sum((c['readings'] for c in cost)), 'slow_QPs': sum((c['slow_QP'] for c in cost)), 'fast_KKT': sum((c['fast_KKT'] for c in cost)), 'sampling_seconds': sum((c['seconds'] for c in cost)), 'wall_seconds': time.perf_counter() - start, 'CPU_seconds': time.process_time() - cpu, 'max_RSS_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'physical_duration': 'UNKNOWN; NOT_RUN', 'fit': 'NONE', 'R1_upstream_wall_seconds': 4.156158004887402}, 'summaries': summary, 'records': records, 'physical_facit': 'UNKNOWN; exact mathematical native facit inherited from R1'}
    write('raw/R2_RAMPS.json', raw)
    write('rounds/R2/results.json', result)
    assert sha(p) == frozen
    print(json.dumps({k: result[k] for k in ['coverage_gate', 'work_recovery_gate', 'practical_gate_1um', 'ideal_metrology_gate', 'cost', 'summaries']}, indent=2))
if __name__ == '__main__':
    main()
