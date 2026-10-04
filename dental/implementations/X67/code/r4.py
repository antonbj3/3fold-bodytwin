"""Continuous joint feasible parameter enclosure, conditional on the two-volume interpolant model.

The rational dual residual correction certifies the LP, including solver stationarity
rounding. It does not certify physical sensor accuracy or model reduction.
"""
import numpy as np, time
from fractions import Fraction as F
from scipy.optimize import linprog
from common import *
from r3 import fixture, TRUE
ENDS = [1, 2, 3, 5, 8, 10, 12, 15, 20]
BOUNDS = [(0.15, 0.6), (0.6, 2.4), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0)]

def interval_rows(a, end, eps):
    i = int(round(end / 0.02))
    d = a[:i + 1]
    t = d[:, 0]
    dt = np.diff(t)
    x = d[:, 1:3] - d[:, 3, None]
    ints = np.sum((x[:-1] + x[1:]) * 0.5 * dt[:, None], axis=0)
    delta = x[-1] - x[0]
    diff = ints[0] - ints[1]
    X = np.array([[delta[0], 0.0, diff, ints[0], 0.0], [0.0, delta[1], -diff, 0.0, ints[1]]])
    B = np.array([[2 * eps, 0.0, 2 * eps * end, eps * end, 0.0], [0.0, 2 * eps, 2 * eps * end, 0.0, eps * end]])
    Q = float(np.sum(d[:-1, 4] * dt))
    W = float(np.sum(d[:-1, 5] * dt))
    return (X, B, Q, W)

def problem(power, eps):
    cal = fixture(TRUE, 'heater', 0.02)
    cal[:, 1:3] = 22 + power * (cal[:, 1:3] - 22)
    cal[:, 4] *= power
    drill = fixture(TRUE, 'drill', 0.02)
    A = []
    b = []
    for end in ENDS:
        for (phase, a) in [('heater', cal), ('drill', drill)]:
            (X, B, Q, W) = interval_rows(a, end, eps)
            for j in range(2):
                low = np.r_[X[j] - B[j], [0.0, 0.0]]
                high = np.r_[-X[j] - B[j], [0.0, 0.0]]
                if phase == 'heater':
                    q = Q if j == 0 else 0.0
                    blo = q * 1.01
                    bhi = -q * 0.99
                else:
                    low[5 + j] = -W * 1.01
                    high[5 + j] = W * 0.99
                    blo = bhi = 0.0
                A.extend([low, high])
                b.extend([blo, bhi])
    A.append([0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0])
    b.append(1.0)
    return (np.array(A), np.array(b), cal, drill)

def frac(x):
    return F(float(x))

def exact_lower(A, b, c, dual):
    lam = [min(frac(v), F(0)) for v in dual]
    r = [frac(c[j]) - sum((lam[i] * frac(A[i, j]) for i in range(len(b))), F(0)) for j in range(len(c))]
    v = sum((lam[i] * frac(b[i]) for i in range(len(b))), F(0))
    v += sum((rj * frac(BOUNDS[j][0] if rj >= 0 else BOUNDS[j][1]) for (j, rj) in enumerate(r)), F(0))
    return (v, r, lam)

def certify(A, b, c):
    res = linprog(c, A_ub=A, b_ub=b, bounds=BOUNDS, method='highs', options={'threads': 4, 'parallel': False})
    if not res.success:
        raise ValueError('LP feasibility failure: ' + res.message)
    (val, r, lam) = exact_lower(A, b, c, res.ineqlin.marginals)
    cert = dict(objective=c.tolist(), dual_nonpositive=[float(v) for v in lam], lower_bound_rational=str(val), stationarity_remainder=[str(v) for v in r], finite_box=[list(v) for v in BOUNDS], solver_objective=float(res.fun), feasible_solution=res.x.tolist(), solver_max_constraint_defect=float(np.max(A @ res.x - b)))
    return (val, cert)

def check(A, b, cert):
    (v, r, lam) = exact_lower(A, b, np.array(cert['objective']), cert['dual_nonpositive'])
    return F(cert['lower_bound_rational']) <= v

def main():
    start = time.perf_counter()
    verify('PREREG_R4.json')
    out = []
    certificates = []
    for power in [1, 4]:
        tmp = fixture(TRUE, 'heater', 0.02)
        maxT = float(np.max(22 + power * (tmp[:, 1:3] - 22)))
        for (label, eps) in [('published_spec', 0.8 + 0.002 * maxT), ('target_0.1C', 0.1), ('target_0.02C', 0.02)]:
            (A, b, cal, drill) = problem(power, eps)
            c = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0])
            (lo, lc) = certify(A, b, c)
            (neg_hi, hc) = certify(A, b, -c)
            assert check(A, b, lc) and check(A, b, hc)
            wrong = dict(lc, lower_bound_rational=str(F(lc['lower_bound_rational']) + F(1, 10)))
            assert not check(A, b, wrong)
            lower = max(0.0, float(np.nextafter(float(lo), -np.inf)))
            upper = min(1.0, float(np.nextafter(float(-neg_hi), np.inf)))
            width = upper - lower
            out.append(dict(heater_power_W=power, sensor_error_C=eps, scenario=label, temperature_max_C=maxT, eta_lower=lower, eta_upper=upper, eta_width=width, acquisition_width_gate_pass=width <= 0.2, certificates_verified=True, wrong_plus0p1_bound_rejected=True, resolution='PER_SURFACE_REGION', scope='Declared two-volume linear-interpolant model only; sensor-to-region and intersample physics absent'))
            certificates.append(dict(scenario=label, power_W=power, A_ub=A.tolist(), b_ub=b.tolist(), lower_certificate=lc, negative_upper_certificate=hc))
    dump('raw/R4_certificates.json', certificates)
    result = dict(round='R4', claim_type='information_link', external_referent={'kind': 'published_dataset', 'locator': 'PMC10299697 Methods: instrument accuracy ±0.8 C ±0.2%, raw/source_manifest.json', 'compared_quantity': 'published instrument specification used as an assumed pointwise error ceiling; not independently verified device accuracy', 'refutes_us': any((not x['acquisition_width_gate_pass'] for x in out if x['scenario'] == 'published_spec'))}, scenarios=out, certificate_form='For lambda<=0, c*x >= lambda*b + min over finite box of (c-A^T lambda)*x; all arithmetic exact rational on frozen binary float coefficients', certificates_sha256=sha(ROOT / 'raw/R4_certificates.json'), physical_remainder='UNKNOWN: actual time history between samples, two-volume reduction, point-to-region observer, heat partition changing during removal, bone injury law', independent_physical_measurements=0, protocol_by_human_bone_class='UNKNOWN', runtime_s=time.perf_counter() - start)
    dump('attempts/R4_results.json', result)
    (ROOT / 'attempts/HANDOFF_R4.md').write_text('R4: published sensor specification applied to global joint thermal/source parameter bounds. Rational dual certificates enclose continuous LP parameter set, conditional on two-volume linear-interpolant model. No real sensor calibration or physiological validity. See R4_results.json for width gates. Next: obtain independent capacities and actual two-sensor/known-power/heater-off/drilling-work traces with the specified temperature precision; validate point-to-region and moving-source transfer on held-out specimens.\n')
    state('R4_DECIDED', 'Sensor-error interval width measured, every LP certificate verified', 'Assemble demo, negative results, physical acquisition contract and graph feedback')
    print('R4', [(d['heater_power_W'], d['scenario'], d['eta_width'], d['acquisition_width_gate_pass']) for d in out])
if __name__ == '__main__':
    main()
