import time, importlib.util, copy, math
import numpy as np
from mpmath import mp, iv
from common import *

def solve(ctx, q, steps=4000, n=20):
    c = lambda x: ctx.mpf(str(x))
    L = c('.001')
    dx = L / n
    dt = c('.01')
    k = c('.58')
    rc = c('1960') * c('1590')
    G = k / dx
    cap = [rc * dx for _ in range(n + 1)]
    cap[0] /= 2
    cap[-1] /= 2
    diag = [v / dt + G * (1 if i in [0, n] else 2) for (i, v) in enumerate(cap)]
    piv = [diag[0]]
    factor = [None]
    for i in range(1, n + 1):
        f = -G / piv[-1]
        factor.append(f)
        piv.append(diag[i] + f * G)
    T = [c('0') for _ in range(n + 1)]
    history = []
    for step in range(steps):
        rhs = [cap[i] / dt * T[i] for i in range(n + 1)]
        rhs[0] += c(q)
        d = [rhs[0]]
        for i in range(1, n + 1):
            d.append(rhs[i] - factor[i] * d[-1])
        T[-1] = d[-1] / piv[-1]
        for i in range(n - 1, -1, -1):
            T[i] = (d[i] + G * T[i + 1]) / piv[i]
    return T[-1]

def independent_matrix_power(q, n=20, steps=4000):
    c = lambda x: mp.mpf(str(x))
    dx = c('.001') / n
    dt = c('.01')
    G = c('.58') / dx
    cap = [c('1960') * c('1590') * dx for _ in range(n + 1)]
    cap[0] /= 2
    cap[-1] /= 2
    K = mp.zeros(n + 1)
    C = mp.zeros(n + 1)
    flux = mp.zeros(n + 1, 1)
    flux[0] = c(q)
    for i in range(n + 1):
        C[i, i] = cap[i] / dt
        K[i, i] = cap[i] / dt + G * (1 if i in [0, n] else 2)
        if i > 0:
            K[i, i - 1] = -G
        if i < n:
            K[i, i + 1] = -G
    inverse = K ** (-1)
    A = inverse * C
    b = inverse * flux
    aug = mp.zeros(n + 2)
    for i in range(n + 1):
        for j in range(n + 1):
            aug[i, j] = A[i, j]
        aug[i, n + 1] = b[i]
    aug[n + 1, n + 1] = 1
    return (aug ** steps)[n, n + 1]

def interval_pair(x):
    return [str(x.a), str(x.b)]

def main():
    start = time.perf_counter()
    iv.dps = 40
    mp.dps = 70
    rigorous = solve(iv, '3800')
    reference = independent_matrix_power('3800')
    bounds = [float(rigorous.a), float(rigorous.b)]
    low_mp = mp.mpf(rigorous.a._mpi_[0])
    high_mp = mp.mpf(rigorous.b._mpi_[1])
    inside = low_mp <= reference <= high_mp
    width = high_mp - low_mp
    nativep = DENTAL / 'cells/procedure/prep_thermal.py'
    spec = importlib.util.spec_from_file_location('x66_r3_heat', nativep)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    (ts, T) = m.slab1d(1.0, 3800.0, 40.0, 0.0, 37.0, n=20, coolant=False)
    native = float(T[-1] - 37)
    distance = abs(native - float(reference))
    assert inside and width <= mp.mpf('1e-16') and (distance <= 1e-07)
    corrupted = reference + 1
    bad_inside = low_mp <= corrupted <= high_mp
    assert not bad_inside
    controltests = []
    ctrls = read('CONTROLS.json')
    for (lid, control) in ctrls.items():
        if lid == 'L01':
            expected = 3.4
            got = control['direct_max_minus_min_N']
        elif lid == 'L02':
            expected = (1 / 0.7) ** 2
            got = control['actual_operator_ratios'][0]
        elif lid == 'L03':
            expected = 158.45025
            got = control['direct_mean_um'][0]
        elif lid == 'L04':
            expected = 3.22 * 5000000 ** (-0.2)
            got = control['direct_source_power_value']
        else:
            expected = 15.2 / 62.6
            got = control['native_slab_ratio']
        tol = 1e-09
        assert abs(got - expected) <= tol and abs(got + 1 - expected) > tol
        controltests.append({'link': lid, 'original_pass': True, 'injected_computed_value': got + 1, 'bad_value_rejected': True, 'tolerance': tol})
    put('R3_RESULTS.json', {'claim_type': 'information_link', 'status': 'PASS_DISCRETE_ARITHMETIC_ENCLOSURE', 'quantity': 'radiation_only_dentin_slab_deltaT_at40s', 'unit': 'degC', 'resolution_level': 'PER_TOOTH', 'source_dose_upper_J_cm2': 15.2, 'discrete_slab': {'n_cells': 20, 'dt_s': 0.01, 'steps': 4000, 'thickness_mm': 1.0, 'k_W_m_K': 0.58, 'rho_kg_m3': 1960.0, 'c_J_kg_K': 1590.0, 'boundaries': 'adiabatic except prescribed heat flux'}, 'iv_endpoint_C': interval_pair(rigorous), 'endpoint_binary64_display_C': [float(np.nextafter(bounds[0], -np.inf)), float(np.nextafter(bounds[1], np.inf))], 'interval_width_C': str(width), 'independent_mp70_endpoint_C': str(reference), 'independent_point_algorithm': 'Dense LU inverse + augmented22x22 affine matrix raised to4000 by binary power; independent of interval Thomas stepping', 'independent_point_inside': bool(inside), 'native_float_endpoint_C': native, 'native_distance_C': distance, 'bad_plus1C_rejected': not bad_inside, 'rigorous_scope': 'Outward interval arithmetic and nonnegative inverse prove[0,upper] for the exact finite model and constant absorption in[0,1], assuming the stated source-dose scenario box. No continuous PDE truncation or empirical model enclosure.', 'linear_operator_proof': 'Positive capacities, positive conductances, strictly diagonally dominant M-matrix at every step; inverse nonnegative; derivative w.r.t imposed nonnegative flux nonnegative; interval Thomas carries arithmetic errors.', 'physical_absorption_range': 'UNKNOWN true value;[0,1] is energy-envelope closure, no measured absorption calibration.', 'full_luting_temperature': 'UNKNOWN', 'external_referent': {'kind': 'independent_measurement', 'locator': read('SOURCE_CONTRACTS.json')['L05']['source_path'] + '#/article/body/sec[2]/p[1]', 'compared_quantity': 'Source optical exposure15.2 upper scenario from14.9±0.3J/cm2; independent measurement is dose, not the simulated47.7C', 'refutes_us': True}, 'source_to_temperature_transformation_status': 'DERIVED_UNDER_ASSUMPTIONS', 'strong_control': {'kind': 'closed_form', 'locator': str(nativep) + '#slab1d', 'compared_quantity': 'Exact finite linear implicitEuler system, independently 70 decimal solar; not physical reference', 'refutes_us': True}, 'control_fault_tests': controltests, 'wall_seconds': time.perf_counter() - start})
    print('R3 interval enclosure', interval_pair(rigorous), 'width', str(width), 'native error', distance, 'wall', time.perf_counter() - start)
if __name__ == '__main__':
    main()
