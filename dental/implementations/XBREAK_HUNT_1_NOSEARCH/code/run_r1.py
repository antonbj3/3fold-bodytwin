"""Frozen ideal-work experiment. All worlds are simulations, not measurements."""
from dental_release.paths import expand as _release_expand
import hashlib, json, resource, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
ROOT = Path(__file__).resolve().parents[1]
X21 = Path(_release_expand('@DENTAL_IMPLEMENTATIONS@/X21'))
FIELD = Path(_release_expand('@DENTAL_IMPLEMENTATIONS@/FALT_TANDLAST'))

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def write(p, x):
    p = ROOT / p
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(x, indent=2, allow_nan=False) + '\n')

class Oracle:
    """Hidden convex support model; consumer below receives no model state."""

    def __init__(self, B, c, g):
        self.B = B
        self.c = c
        self.H = B.T @ (c[:, None] * B)
        self.g = g
        self.calls = 0
        self.seconds = 0.0
        self.max_kkt = 0.0
        self.max_conservation = 0.0
        self.solver_statuses = []

    def solve(self, T, g=None):
        t = time.perf_counter()
        self.calls += 1
        g = self.g if g is None else np.asarray(g)
        if T == 0:
            return dict(V=0.0, d=float(g.min()), force=np.zeros(len(g)), q=np.zeros(len(self.c)), active=[])
        n = len(g)
        scale = max(float(np.ptp(g)), float(T * np.max(self.H)), 0.0001)
        A = T * self.H / scale
        b = g / scale
        res = minimize(lambda x: 0.5 * x @ A @ x + b @ x, np.ones(n) / n, jac=lambda x: A @ x + b, constraints={'type': 'eq', 'fun': lambda x: x.sum() - 1, 'jac': lambda x: np.ones(n)}, bounds=[(0, 1)] * n, method='SLSQP', options={'ftol': 1e-14, 'maxiter': 2500})
        self.solver_statuses.append({'success': bool(res.success), 'message': str(res.message)})
        force = res.x * T
        active = np.flatnonzero(force > 1e-07)
        K = np.block([[self.H[np.ix_(active, active)], -np.ones((len(active), 1))], [np.ones((1, len(active))), np.zeros((1, 1))]])
        rhs = np.r_[-g[active], T]
        pol = np.linalg.lstsq(K, rhs, rcond=1e-13)[0]
        trial = np.zeros(n)
        trial[active] = pol[:-1]
        grad = self.H @ trial + g
        d = float(pol[-1])
        if trial.min() >= -1e-07 and grad.min() >= d - 1e-07:
            force = trial
        else:
            grad = self.H @ force + g
            d = float(np.mean(grad[active]))
        residual = max(float(np.max(np.abs(grad[active] - d))), max(0.0, float(d - grad.min())), max(0.0, float(-force.min())))
        conserve = abs(float(force.sum() - T))
        self.max_kkt = max(self.max_kkt, residual)
        self.max_conservation = max(self.max_conservation, conserve)
        self.seconds += time.perf_counter() - t
        if residual > 1e-07 or conserve > 1e-07:
            raise ArithmeticError(f'KKT residual {residual}, conservation {conserve}')
        return dict(V=float(0.5 * force @ self.H @ force + g @ force), d=d, force=force, q=self.B @ force, active=active.tolist())

def consumer(vminus, vzero, vplus, h, T):
    """Only bounded observations, step, and total load cross this interface."""
    lo = max(0.0, (vplus - vzero) / h)
    hi = min(T, (vzero - vminus) / h)
    return (lo, hi)

def model(case, scenario):
    p = X21 / 'raw/refined_certificates' / f'{case}.json'
    rows = json.loads(p.read_text())['pairs']
    teeth = sorted({int(r[z]) for r in rows for z in ['upper_fdi', 'lower_fdi']})
    B = np.zeros((len(teeth), len(rows)))
    idx = {t: i for (i, t) in enumerate(teeth)}
    for (e, r) in enumerate(rows):
        B[idx[int(r['upper_fdi'])], e] = 1
        B[idx[int(r['lower_fdi'])], e] = 1
    g = np.array([sum(r['minimum_gap_enclosure_mm']) / 2 for r in rows])
    if scenario == 'constant_750':
        k = np.full(len(teeth), 750.0)
    elif scenario == 'seed_20261003_loguniform_150_5000':
        k = np.exp(np.random.default_rng(20261003).uniform(np.log(150), np.log(5000), len(teeth)))
    else:
        k = np.exp(np.random.default_rng(20261004).uniform(np.log(30), np.log(30000), len(teeth)))
    return (Oracle(B, 1 / k, g), teeth, rows, p)

def main():
    start = time.perf_counter()
    cpu = time.process_time()
    prereg = ROOT / 'PREREG_R1.json'
    cfg = json.loads(prereg.read_text())
    frozen = sha(prereg)
    write('raw/PREREG_R1_FROZEN.json', {'path': str(prereg), 'sha256': frozen})
    sys.path.insert(0, str(FIELD / 'code'))
    import model as field
    T = cfg['T_N']
    records = []
    worlds = []
    observations = []
    exact_checks = []
    targets = list(range(11, 19)) + list(range(21, 29))
    count = 0
    oracle_s = 0.0
    maxkkt = 0.0
    for case in cfg['cases']:
        for scenario in cfg['support_scenarios']:
            (oracle, teeth, edges, p) = model(case, scenario)
            base = oracle.solve(T)
            ee = [(int(r['upper_fdi']), int(r['lower_fdi']), field.Q(str(g)), field.Q(str(g))) for (r, g) in zip(edges, oracle.g)]
            exact = field.exact_point(ee, list(map(lambda x: field.Q(str(x)), oracle.g)), field.Q(T), {t: field.Q(str(c)) for (t, c) in zip(teeth, oracle.c)})
            qexact = {int(t): float(q) for (t, q) in exact['tooth_force_N'].items()}
            maxdiff = max((abs(base['q'][i] - qexact[t]) for (i, t) in enumerate(teeth)))
            exact_checks.append({'case': case, 'scenario': scenario, 'exact_KKT': exact['exact_KKT'], 'max_force_difference_N': maxdiff, 'closure_exact_mm': str(exact['closure_mm'])})
            if maxdiff > 1e-05:
                raise ArithmeticError('Independent native force disagreement')
            worlds.append({'case': case, 'scenario': scenario, 'teeth': teeth, 'g_mm': oracle.g.tolist(), 'c_mm_per_N': oracle.c.tolist(), 'native_force_N': {str(t): qexact[t] for t in teeth}, 'native_V_Nmm': base['V'], 'native_d_mm': base['d'], 'source_path': str(p), 'source_sha256': sha(p)})
            for tooth in targets:
                direction = oracle.B[teeth.index(tooth)] if tooth in teeth else np.zeros(len(oracle.g))
                q = qexact.get(tooth, 0.0)
                for h in cfg['gap_steps_mm']:
                    minus = oracle.solve(T, oracle.g - h * direction)
                    plus = oracle.solve(T, oracle.g + h * direction)
                    (lo, hi) = consumer(minus['V'], base['V'], plus['V'], h, T)
                    control = sorted([(plus['V'] - base['V']) / h, (base['V'] - minus['V']) / h])
                    control = (max(0.0, control[0]), min(T, control[1]))
                    records.append({'case': case, 'scenario': scenario, 'FDI': tooth, 'h_mm': h, 'facit_N': q, 'lower_N': lo, 'upper_N': hi, 'width_pp': 100 * (hi - lo) / T, 'covered': lo - 1e-05 <= q <= hi + 1e-05, 'same_information_control_difference_N': max(abs(lo - control[0]), abs(hi - control[1])), 'contact_switched': minus['active'] != base['active'] or plus['active'] != base['active']})
                    observations.append({'case': case, 'scenario': scenario, 'FDI': tooth, 'h_mm': h, 'T_N': T, 'minus_work_Nmm': minus['V'], 'baseline_work_Nmm': base['V'], 'plus_work_Nmm': plus['V']})
            count += oracle.calls
            oracle_s += oracle.seconds
            maxkkt = max(maxkkt, oracle.max_kkt)
            worlds[-1]['solver_failed_statuses_with_passing_KKT'] = sum((not s['success'] for s in oracle.solver_statuses))
            write('raw/R1_PROGRESS.json', {'worlds_completed': worlds, 'records': records, 'exact_checks': exact_checks})
            print(case, scenario, 'completed', flush=True)
    summary = []
    for h in cfg['gap_steps_mm']:
        rr = [r for r in records if r['h_mm'] == h]
        ww = np.array([r['width_pp'] for r in rr])
        summary.append({'h_mm': h, 'n': len(rr), 'coverage': sum((r['covered'] for r in rr)), 'width_le5pp': int(sum(ww <= 5)), 'fraction_le5pp': float(np.mean(ww <= 5)), 'median_width_pp': float(np.median(ww)), 'max_width_pp': float(ww.max()), 'contact_switch_queries': sum((r['contact_switched'] for r in rr))})
    forces = [[T * 1000 / 1100, T * 100 / 1100], [T * 100 / 1100, T * 1000 / 1100]]
    obstruction = {'two_disjoint_contacts': True, 'g_mm': [0.0, 0.0], 'effective_stiffnesses_N_per_mm': [[1000, 100], [100, 1000]], 'each_tooth_stiffness_N_per_mm': [[2000, 200], [200, 2000]], 'native_curve_in_both': 'd(T)=T/1100 mm for all T>=0', 'native_force_N': forces, 'target_difference_pp': 100 * abs(forces[0][0] - forces[1][0]) / T}
    out = {'round': 'R1', 'claim_type': cfg['claim_type'], 'status': 'PENDING_INDEPENDENT_REVIEW', 'prereg_sha256': frozen, 'coverage_gate': all((r['covered'] for r in records)), 'useful_width_gate': any((r['fraction_le5pp'] >= 0.8 for r in summary)), 'numerical_gate': maxkkt <= 1e-07, 'information_obstruction_gate': obstruction['target_difference_pp'] >= 20, 'summaries': summary, 'same_information_control': 'same interval; no algorithm advantage', 'physical_validation': 'UNKNOWN; ideal work oracle; synthetic supports and unvalidated projected geometry', 'cost': {'all_equilibrium_calls': count, 'oracle_seconds': oracle_s, 'wall_seconds': time.perf_counter() - start, 'CPU_seconds': time.process_time() - cpu, 'max_RSS_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'fit': 'NONE', 'physical_work_acquisition': 'NOT_RUN', 'reasoning_token_cost': 'UNKNOWN'}, 'max_KKT_residual_mm': maxkkt, 'independent_native_checks': exact_checks, 'obstruction': obstruction, 'records': records}
    write('raw/R1_OBSERVATIONS.json', observations)
    write('raw/R1_WORLDS.json', worlds)
    write('rounds/R1/results.json', out)
    assert sha(prereg) == frozen
    print(json.dumps({k: out[k] for k in ['coverage_gate', 'useful_width_gate', 'summaries', 'cost']}, indent=2))
if __name__ == '__main__':
    main()
