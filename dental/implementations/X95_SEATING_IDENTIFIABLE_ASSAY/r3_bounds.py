"""Exact rational resistance enclosures; no physical rheogram is manufactured."""
import json, time
import sympy as s
from mpmath import iv
import numpy as np
from common import P, check, dump, state
from assay import hull, interval

def run(output=None):
    check('R3')
    t0 = time.perf_counter()
    n = 4
    root = P if output is None else __import__('pathlib').Path(output)
    V = [(s.Rational(x, n), s.Rational(y, n)) for y in range(n + 1) for x in range(n + 1)]
    faces = []
    for y in range(n):
        for x in range(n):
            a = y * (n + 1) + x
            b = a + n + 1
            faces.extend([[a, a + 1, b + 1], [a, b + 1, b]])
    local = []
    loads = s.zeros(len(V), 1)
    for f in faces:
        M = s.Matrix([[1, *V[i]] for i in f])
        area = abs(M.det()) / 2
        coeff = M.inv()
        grad = coeff[1:, :].T
        L = area * (grad * grad.T) / 12
        assert all((L[i, j] <= 0 for i in range(3) for j in range(3) if i != j))
        for i in f:
            loads[i] += area / 3
        local.append((f, L, sum((V[i][0] for i in f)) / 3, area))
    free = [i for (i, (x, y)) in enumerate(V) if 0 < x < 1 and 0 < y < 1]
    b = loads.extract(free, [0])

    def resistance(z, delta):
        K = s.zeros(len(V))
        mean = 0
        for (f, L, x, area) in local:
            h = z + delta * (2 * x - 1)
            assert h > 0
            mean += area * h
            for i in range(3):
                for j in range(3):
                    K[f[i], f[j]] += h ** 3 * L[i, j]
        p = K.extract(free, free).inv() * b
        assert min(p) >= 0
        return ((b.T * p)[0], mean)
    (D0, _) = resistance(s.Rational(1), s.Rational(0))
    rows = []
    for z in [s.Rational(1), s.Rational(5, 4), s.Rational(2)]:
        for d in [s.Rational(0), s.Rational(1, 10), s.Rational(1, 5)]:
            (D, mean) = resistance(z, d)
            lo = D0 / (z + d) ** 3
            hi = D0 / (z - d) ** 3
            valid = lo <= D <= hi
            rows.append(dict(z=str(z), delta=str(d), resistance_exact=str(D), resistance=float(D), lower_exact=str(lo), upper_exact=str(hi), lower=float(lo), upper=float(hi), bound_valid=bool(valid), summary_mean_exact=str(mean), mean_identity_error_exact=str(mean - z)))
            assert valid and mean == z
            if d == 0:
                assert D == D0 / z ** 3
    base = rows[0]
    tilt = rows[2]
    diff = s.Rational(tilt['resistance_exact']) - s.Rational(base['resistance_exact'])
    faults = []
    for (badlo, badhi, name) in [(2 * D0, D0, '2x_lower_bound'), (D0, D0 / 2, '0.5x_upper_bound')]:
        rejected = not badlo <= D0 <= badhi
        faults.append(dict(name=name, rejected=bool(rejected)))
        assert rejected
    z0 = interval(1.25)
    d = interval(0.2)
    (jlo, jhi) = (0.05, 0.15)
    lower = hull(1 / iv.sqrt(1 / (z0 + d) ** 2 + 2 * interval(jhi)) - d)[0]
    upper = hull(1 / iv.sqrt(1 / (z0 - d) ** 2 + 2 * interval(jlo)) + d)[1]
    assert lower > 0.2
    result = dict(claim_type='capability', outcome='LOCAL_TILT_AND_FUTURE_DOSE_ENCLOSURE_DERIVED', grid=dict(cells_per_axis=n, free_unknowns=len(free), triangles=len(faces), positive_edge_projection='inactive exactly'), uniform_D0_exact=str(D0), states=rows, bound_violations=sum((not r['bound_valid'] for r in rows)), mean_sufficiency_witness=dict(summary_identity_error_exact='0', summary_identity_error_float=0.0, summary_quantity='Area-weighted film height', states=['z=1,delta=0', 'z=1,delta=1/5'], downstream_resistance_difference_exact=str(diff), downstream_relative_difference=float(diff / D0), smallest_extension='Locked local film deviation interval delta suffices for a resistance enclosure; full local field for an exact value'), held_pulse_bound=dict(initial_z=1.25, delta=0.2, force_fluidity_dose_over_D0_interval=[jlo, jhi], final_z_interval=[lower, upper], conditional_validity='Full Newtonian film; locked pose; positive forcing; known future dose interval; lower>delta', physical_future_dose='UNKNOWN_MISSING_INDEPENDENT_AGE_SHEAR_TEMPERATURE_BOUND_RHEOGRAM'), faults=faults, seconds=time.perf_counter() - t0, resolution='PER_SURFACE_REGION', timescale='HANDOVER', physical_status='UNKNOWN', external_referent=dict(kind='closed_form', locator='doi:10.1098/rstl.1886.0005; ASSAY_DERIVATION_R3.md explicit variational proof', compared_quantity='Positive Reynolds resistance and nonlinear comparison envelope; exact rational fixture is algebraic control only', refutes_us=False), limitations=['No empirical crown response or physical rheogram', 'No general obtuse-mesh clipping enclosure', 'No time-jitter enclosure', 'No changing tilt, contact, particles or gel force'], review_state='PENDING_INDEPENDENT_REVIEW')
    assert diff != 0
    dump(root / 'RESULTS_R3.json', result)
    state('R3_EXACT_BOUNDS_DECIDED', result['outcome'], 'Package bench CSV/acquisition/scorer, source comparison and independent replay', root)
    print(result['outcome'], 'exact violations=', result['bound_violations'], 'same mean difference=', float(diff / D0))
    return result
if __name__ == '__main__':
    run()
