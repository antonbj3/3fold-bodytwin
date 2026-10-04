"""Read-only reuse of the prior physical operator; no statistical field fabrication."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
import sys, importlib.util, hashlib, json, time
import numpy as np
from scipy.sparse.linalg import spsolve
P = Path(__file__).resolve().parent
SOURCE = Path(_release_expand('@DENTAL_INPUT_ROOT@/artifacts/LANE_NEXT_E_CEMENT_SQUEEZE/squeeze.py'))

def load_operator():
    sys.dont_write_bytecode = True
    spec = importlib.util.spec_from_file_location('bte1_readonly', SOURCE)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

def hydraulic_input_from_group_mean(row):
    return dict(status='UNKNOWN_MISSING_SPATIAL_FILM', row_id=row['row_id'], reported_mean_um=row['measured_mean_um'], observation_state=row['state'], missing=['same-object paired die and intaglio geometry', 'regional normals and film distribution/connectivity', 'dry versus replica versus cemented seating separation', 'cement viscosity as shear/time/temperature law', 'fill state, vent geometry and boundary conditions'], no_additive_seating_prediction=True)

def strip(m, n, low, high, kind):
    L = W = 0.006
    mu = 100.0
    V = np.array([[x, y, 0] for y in np.linspace(0, W, n + 1) for x in np.linspace(0, L, n + 1)])
    F = []
    for iy in range(n):
        for ix in range(n):
            a = iy * (n + 1) + ix
            b = a + n + 1
            F.extend([[a, a + 1, b + 1], [a, b + 1, b]])
    s = m.Surface(V, np.array(F))
    coordinate = s.center[:, 0 if kind == 'series' else 1]
    h = np.where(coordinate < L / 2, low, high)
    K = s.stiffness(h ** 3 / (12 * mu))
    left = np.flatnonzero(np.isclose(V[:, 0], 0, atol=1e-15))
    right = np.flatnonzero(np.isclose(V[:, 0], L, atol=1e-15))
    free = np.setdiff1d(np.arange(len(V)), np.r_[left, right])
    p = np.zeros(len(V))
    p[left] = 1.0
    p[free] = spsolve(K[free][:, free], -(K[free][:, left] @ np.ones(len(left))))
    reaction = K @ p
    Q = float(reaction[left].sum())
    out = -float(reaction[right].sum())
    exact = W / L / (12 * mu) * ((0.5 / low ** 3 + 0.5 / high ** 3) ** (-1) if kind == 'series' else 0.5 * low ** 3 + 0.5 * high ** 3)
    return dict(kind=kind, n=n, mean_film_um=float(np.average(h, weights=s.area) * 1000000.0), conductance_m3_Pa_s=Q, exact_m3_Pa_s=exact, relative_error=abs(Q / exact - 1), mass_relative_error=abs(Q - out) / Q, pressure_min=float(p.min()), pressure_max=float(p.max()), faces=len(F))

def run():
    from predict_r1 import check, sha, dump
    check('R4')
    started = time.perf_counter()
    rows = json.loads((P / 'measurements.json').read_text())
    q = [d for d in rows if d['study'] == 'PMC10333096' and d['region'] == 'axial']
    lo = min(q, key=lambda d: d['measured_mean_um'])
    hi = max(q, key=lambda d: d['measured_mean_um'])
    low = lo['measured_mean_um'] * 1e-06
    high = hi['measured_mean_um'] * 1e-06
    m = load_operator()
    results = [strip(m, n, low, high, k) for n in [16, 32, 64] for k in ['series', 'parallel']]
    mean = (low + high) / 2
    uniform = mean ** 3 / 1200
    alias = results[-1]['conductance_m3_Pa_s'] / results[-2]['conductance_m3_Pa_s']
    numerical = all((d['relative_error'] <= 0.02 and d['mass_relative_error'] <= 1e-08 for d in results))
    equalmean = max((abs(d['mean_film_um'] - mean * 1000000.0) / (mean * 1000000.0) for d in results)) <= 1e-12
    tests = []
    for d in results:
        bad = 2 * d['conductance_m3_Pa_s']
        error = abs(bad / d['exact_m3_Pa_s'] - 1)
        assert error > 0.02
        tests.append(dict(control='BTE1 numerical conductance', kind=d['kind'], n=d['n'], injection='2x conductance', rejected=True, error=error))
        bad = 2 * d['exact_m3_Pa_s']
        error = abs(d['conductance_m3_Pa_s'] / bad - 1)
        assert error > 0.02
        tests.append(dict(control='closed-form conductance', kind=d['kind'], n=d['n'], injection='2x analytic value', rejected=True, error=error))
        error = abs(2 * uniform / d['exact_m3_Pa_s'] - 1)
        assert error > 0.02
        tests.append(dict(control='uniform-mean information ablation', kind=d['kind'], n=d['n'], injection='2x mean-film conductance', rejected=True, error=error))
    result = dict(round='R4', outcome='MODEL_MEAN_INSUFFICIENCY_DEMONSTRATED' if numerical and equalmean and (alias >= 1.25) else 'FAIL', numerical_verification_pass=numerical, same_mean_pass=equalmean, aliasing_conductance_ratio=alias, low_um=low * 1000000.0, high_um=high * 1000000.0, source_rows=[lo['row_id'], hi['row_id']], constructed_fields=True, empirical_spatial_field_validation=False, uniform_mean_conductance_m3_Pa_s=uniform, runs=results, corruption_tests=tests, operator_path=str(SOURCE), operator_sha256=sha(SOURCE), seconds=time.perf_counter() - started, external_referent=json.loads((P / 'PREREG_R4.json').read_text())['external_referent'], review_state='PENDING_INDEPENDENT_REVIEW')
    dump(P / 'RESULTS_R4.json', result)
    dump(P / 'CEMENT_FIELD_PORT.json', dict(status='UNKNOWN_MISSING_SPATIAL_FILM', records=[hydraulic_input_from_group_mean(d) for d in rows if d['region'] in ['marginal', 'axial', 'occlusal']], source_operator_sha256=sha(SOURCE), warning='Do not add a simulated seating lift to a source film that already contains seating. Region means are not normal gap fields.'))
    print(result['outcome'], 'conductance ratio', alias, 'max relative error', max((d['relative_error'] for d in results)))
if __name__ == '__main__':
    run()
