"""Rate-state adverse fixture, independent ODE facit, bounded reaction inversion."""
import json, time, resource, warnings
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import linprog, OptimizeWarning
from run_r1 import ROOT, sha, write

def dynamic(K, G, T0, T1, duration):
    rate = (T1 - T0) / duration
    tau = G / K
    t = np.linspace(0, duration, 1001)
    d = T0 / K + rate / K * (t - tau * (1 - np.exp(-t / tau)))
    velocity = rate / K * (1 - np.exp(-t / tau))
    return (t, d, velocity)

def invert_reaction(observations, h, rates, error, d, v, derr, verr):
    """Native local-force consumer. No hidden constitutive coefficients enter."""
    A = np.column_stack([np.full(2, h), np.array(rates)])
    vertices = []
    for s0 in [-1, 1]:
        for s1 in [-1, 1]:
            z = np.linalg.solve(A, np.array(observations) + error * np.array([s0, s1]))
            if np.min(z) < 0:
                raise ArithmeticError('Negative coefficient corner; a clipped polygon is needed')
            vertices.append(z)
    low = min((z @ np.array([d - derr, v - verr]) for z in vertices))
    high = max((z @ np.array([d + derr, v + verr]) for z in vertices))
    AA = np.r_[A, -A]
    bb = np.r_[np.array(observations) + error, -np.array(observations) + error]
    lpmin = linprog([d - derr, v - verr], A_ub=AA, b_ub=bb, bounds=[(0, None), (0, None)], method='highs', options={'threads': 4})
    lpmax = linprog([-d - derr, -v - verr], A_ub=AA, b_ub=bb, bounds=[(0, None), (0, None)], method='highs', options={'threads': 4})
    if not lpmin.success or not lpmax.success:
        raise ArithmeticError('Reaction inversion LP failed')
    return (float(low), float(high), max(abs(float(low) - lpmin.fun), abs(float(high) + lpmax.fun)), [z.tolist() for z in vertices])

def perturbation(v, beta=5000):
    return beta * v * (v - 0.001) * (v - 0.002) * (v - 0.003)

def main():
    start = time.perf_counter()
    cpu = time.process_time()
    p = ROOT / 'PREREG_R5.json'
    cfg = json.loads(p.read_text())
    frozen = sha(p)
    write('raw/PREREG_R5_FROZEN.json', {'sha256': frozen, 'file': str(p)})
    warnings.filterwarnings('ignore', category=OptimizeWarning, message='Unrecognized options detected.*')
    k = np.array(cfg['k_N_per_mm'], dtype=float)
    K = float(k.sum())
    ramp = cfg['load_ramp']
    T0 = ramp['initial_N']
    T1 = ramp['final_N']
    dur = ramp['duration_s']
    h = cfg['probe']['gap_height_mm']
    rates = cfg['probe']['rates_mm_per_s']
    valrate = cfg['probe']['heldout_rate_mm_per_s']
    derr = cfg['native_observation_bounds']['displacement_mm']
    verr = cfg['native_observation_bounds']['velocity_mm_per_s']
    worlds = []
    records = []
    raw = []
    max_solver = 0.0
    maxlp = 0.0
    old_errors = []
    nonlinear_checks = []
    for (wi, eta_raw) in enumerate(cfg['eta_worlds_Ns_per_mm']):
        eta = np.array(eta_raw, dtype=float)
        G = float(eta.sum())
        (t, d, v) = dynamic(K, G, T0, T1, dur)
        force = k[:, None] * d + eta[:, None] * v
        sol = solve_ivp(lambda ts, x: [(T0 + (T1 - T0) * ts / dur - K * x[0]) / G], (0, dur), [T0 / K], t_eval=t, rtol=1e-11, atol=1e-13)
        independent_v = (T0 + (T1 - T0) * t / dur - K * sol.y[0]) / G
        independent_force = k[:, None] * sol.y[0] + eta[:, None] * independent_v
        error = float(np.max(np.abs(independent_force - force)))
        max_solver = max(max_solver, error)
        if error > 1e-07:
            raise ArithmeticError('Independent ODE disagrees')
        if np.min(force) <= 0:
            raise ArithmeticError('Contact loss in preload fixture')
        conservation = float(np.max(np.abs(force.sum(axis=0) - (T0 + (T1 - T0) * t / dur))))
        if conservation > 1e-07:
            raise ArithmeticError('Force conservation failed')
        work_difference = [quad(lambda load: h * k[i] / K, T0, T1, epsabs=1e-12)[0] for i in range(2)]
        old_estimate = k * T0 / K + np.array(work_difference) / h
        old_error = 100 * np.abs(old_estimate - force[:, -1]) / T1
        old_errors.extend(old_error.tolist())
        polynomial = perturbation(v)
        nonlinear = force + np.vstack([polynomial, -polynomial])
        dissipation_lower = float(np.min(eta) - 5000 * float(v.max()) ** 3)
        if dissipation_lower < 0 or np.min(nonlinear) <= 0:
            raise ArithmeticError('Passive adversary invalid')
        worlds.append({'world': wi, 'k_N_per_mm': k.tolist(), 'eta_Ns_per_mm': eta.tolist(), 'native_d_mm': float(d[-1]), 'native_velocity_mm_per_s': float(v[-1]), 'native_force_N': force[:, -1].tolist(), 'old_work_estimate_N': old_estimate.tolist(), 'old_error_pp': old_error.tolist(), 'nonlinear_force_N': nonlinear[:, -1].tolist(), 'nonlinear_native_delta_N': float(polynomial[-1]), 'force_conservation_error_N': conservation, 'independent_ODE_max_error_N': error, 'passive_viscous_coefficient_lower_bound_Ns_per_mm': dissipation_lower})
        raw.append({'world': wi, 'time_s': t.tolist(), 'closure_mm': d.tolist(), 'velocity_mm_per_s': v.tolist(), 'native_force_N': force.tolist(), 'nonlinear_force_N': nonlinear.tolist()})
        for err in cfg['reaction_error_bounds_N']:
            for i in range(2):
                truth = k[i] * h + eta[i] * np.array(rates)
                measured = truth + err * np.array([0.4 * np.sin(wi + i + 1), 0.7 * np.cos(wi - i + 1)])
                (lo, hi, lpdiff, vertices) = invert_reaction(measured, h, rates, err, float(d[-1]), float(v[-1]), derr, verr)
                maxlp = max(maxlp, lpdiff)
                covered = lo - 1e-05 <= force[i, -1] <= hi + 1e-05
                nonlinear_calibration_defect = max((abs(float(perturbation(rate))) for rate in [0, *rates, valrate]))
                nonlin_covered = lo - 1e-05 <= nonlinear[i, -1] <= hi + 1e-05
                records.append({'world': wi, 'support': i, 'reaction_error_N': err, 'observed_reaction_N': measured.tolist(), 'true_reaction_N': truth.tolist(), 'lower_N': lo, 'upper_N': hi, 'facit_N': float(force[i, -1]), 'width_pp': 100 * (hi - lo) / T1, 'covered': bool(covered), 'LP_difference_N': float(lpdiff), 'coefficient_vertices': vertices, 'heldout_rate_response_N': float(k[i] * h + eta[i] * valrate), 'nonlinear_calibration_defect_N': nonlinear_calibration_defect, 'nonlinear_native_force_N': float(nonlinear[i, -1]), 'nonlinear_force_in_linear_interval': bool(nonlin_covered), 'unverified_constitutive_law_consumer_status': 'UNKNOWN'})
                if err == 0.005:
                    nonlinear_checks.append(not nonlin_covered)
    swapped_pp = 100 * abs(worlds[0]['native_force_N'][0] - worlds[1]['native_force_N'][0]) / T1
    summary = []
    for err in cfg['reaction_error_bounds_N']:
        rr = [r for r in records if r['reaction_error_N'] == err]
        summary.append({'reaction_error_N': err, 'n': len(rr), 'covered': sum((r['covered'] for r in rr)), 'n_width_le5pp': sum((r['width_pp'] <= 5 for r in rr)), 'max_width_pp': max((r['width_pp'] for r in rr)), 'nonlinear_force_exclusions': sum((not r['nonlinear_force_in_linear_interval'] for r in rr))})
    out = {'round': 'R5', 'claim_type': 'information_link', 'status': 'PENDING_INDEPENDENT_REVIEW', 'prereg_sha256': frozen, 'indistinguishability_gate': swapped_pp >= 20, 'swapped_native_force_difference_pp': swapped_pp, 'native_total_curve_identical': True, 'old_work_failure_gate': max(old_errors) > 5, 'rate_port_gate': all((r['covered'] and r['width_pp'] <= 5 for r in records if r['reaction_error_N'] == 0.005)), 'ODE_gate': max_solver <= 1e-07, 'nonlinear_adverse_gate': any(nonlinear_checks), 'max_same_information_LP_difference_N': maxlp, 'worlds': worlds, 'summaries': summary, 'records': records, 'cost': {'wall_seconds': time.perf_counter() - start, 'CPU_seconds': time.process_time() - cpu, 'max_RSS_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'LP_calls': 2 * len(records), 'hypothetical_rate_probe_travel_seconds_per_support': sum((h / rate for rate in rates)), 'physical_probe_setup_settling_metrology': 'UNKNOWN; NOT_RUN', 'fit': 'two local coefficients per support, all force-reading cost charged'}, 'physical_native_force_certificate': 'UNKNOWN; law regime and measurement calibration absent', 'diagnosis': 'Even calibration plus a held-out intermediate rate cannot license extrapolation to native velocity without a constitutive law. The passive nonlinear adversary matches all three rates exactly.'}
    write('raw/R5_RAMPS.json', raw)
    write('rounds/R5/results.json', out)
    assert sha(p) == frozen
    print(json.dumps({k: out[k] for k in ['indistinguishability_gate', 'swapped_native_force_difference_pp', 'old_work_failure_gate', 'rate_port_gate', 'ODE_gate', 'nonlinear_adverse_gate', 'summaries', 'cost']}, indent=2))
if __name__ == '__main__':
    main()
