"""Finite normal cone; field exact verifier is the trusted arithmetic component.

Float LPs only propose witnesses. Six affine charts cover all nonzero directions.
Quantization envelope certifies the mathematical mesh-normal model, not scanner truth.
"""
import sys, time, json, hashlib
from pathlib import Path
from fractions import Fraction as Q
import numpy as np
from scipy.optimize import linprog
import sympy as sp
FIELD = Path(__file__).resolve().parent
sys.path.insert(0, str(FIELD / 'vendor'))
sys.path.insert(0, str(FIELD / 'code'))
import highspy
from dual_witness_v1 import check_dual, check_primal
SCALE = 10 ** 6

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda : f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def highs(A, b, c=None):
    A = np.asarray(A, float)
    b = np.asarray(b, float)
    (m, n) = A.shape
    h = highspy.Highs()
    h.setOptionValue('output_flag', False)
    h.setOptionValue('threads', 4)
    h.setOptionValue('presolve', 'off')
    h.setOptionValue('primal_feasibility_tolerance', 1e-09)
    h.setOptionValue('dual_feasibility_tolerance', 1e-09)
    h.addVars(n, np.full(n, -highspy.kHighsInf), np.full(n, highspy.kHighsInf))
    if c is not None:
        h.changeColsCost(n, np.arange(n, dtype=np.int32), np.asarray(c, float))
    h.addRows(m, np.full(m, -highspy.kHighsInf), b, m * n, np.arange(0, (m + 1) * n, n, dtype=np.int32), np.tile(np.arange(n, dtype=np.int32), m), A.ravel())
    h.run()
    status = h.modelStatusToString(h.getModelStatus())
    ray = h.getDualRay()
    return {'status': status, 'x': h.getSolution().col_value, 'ray': -ray[2] if ray[1] else None}

def reconstruct(A, b, ray):
    """Exact support reconstruction; never accepts a small numerical residual."""
    support = np.flatnonzero(np.abs(ray) > 1e-09).tolist()
    if not support:
        return None
    M = sp.Matrix([[sp.Rational(str(A[i][j])) for i in support] for j in range(len(A[0]))])
    for v in M.nullspace():
        for sign in [1, -1]:
            y = [Q(str(sign * q)) for q in v]
            if min(y) < 0:
                continue
            As = [A[i] for i in support]
            bs = [b[i] for i in support]
            d = check_dual('inequality', As, bs, y)
            if d.status == 'NEJ':
                return {'rows': support, 'A': [[str(z) for z in r] for r in As], 'b': [str(z) for z in bs], 'y': [str(z) for z in y], 'separation': str(d.separation), 'verified': True}
    return None

def chart(N, j, s, relax=Q(0)):
    A = [[-Q(int(v), SCALE) for v in r] for r in N]
    b = [relax] * len(A)
    for k in range(3):
        for sign in [1, -1]:
            r = [Q(0)] * 3
            r[k] = Q(sign)
            A.append(r)
            b.append(Q(1))
    r = [Q(0)] * 3
    r[j] = Q(-s)
    A.append(r)
    b.append(Q(-1))
    return (A, b)

def classify(normals, rounding_envelope=True):
    t0 = time.perf_counter()
    unit = np.asarray(normals, float)
    unit = unit / np.linalg.norm(unit, axis=1)[:, None]
    N = np.rint(unit * SCALE).astype(np.int64)
    N = np.unique(N, axis=0)
    relax = Q(2, SCALE) if rounding_envelope else Q(0)
    A = np.column_stack([-N / SCALE, np.ones(len(N))])
    opt = linprog([0, 0, 0, -1], A_ub=A, b_ub=np.zeros(len(N)), bounds=[(-1, 1)] * 3 + [(None, None)], method='highs', options={'dual_feasibility_tolerance': 1e-09, 'primal_feasibility_tolerance': 1e-09})
    prep_s = time.perf_counter() - t0
    if opt.success and opt.x[3] > 4e-06:
        d = [Q(float(x)).limit_denominator(10 ** 9) for x in opt.x[:3]]
        Af = [[-Q(int(v), SCALE) for v in r] for r in N]
        bf = [-relax] * len(N)
        tick = time.perf_counter()
        valid = check_primal('inequality', Af, bf, d).status == 'JA'
        vs = time.perf_counter() - tick
        if valid:
            df = np.array([float(x) for x in d])
            margin = float(np.min(unit @ df) / np.linalg.norm(df))
            return {'status': 'YES_STRICT', 'direction': [str(x) for x in d], 'unit_direction': (df / np.linalg.norm(df)).tolist(), 'min_dot_input_normals': margin, 'guaranteed_margin_exported': float(np.min(N / SCALE @ df) - float(relax)) / np.linalg.norm(df), 'unique_normals': len(N), 'certificates': [], 'discovery_s': prep_s, 'exact_check_s': vs, 'total_s': time.perf_counter() - t0, 'charts': []}
    certs = []
    charts = []
    candidate = None
    solve_s = 0
    check_s = 0
    for j in range(3):
        for s in [1, -1]:
            tick = time.perf_counter()
            (Ar, br) = chart(N, j, s, relax)
            v = highs(Ar, br)
            ds = time.perf_counter() - tick
            solve_s += ds
            tick = time.perf_counter()
            cert = reconstruct(Ar, br, v['ray']) if v['ray'] is not None else None
            check_s += time.perf_counter() - tick
            charts.append({'axis': j, 'sign': s, 'solver': v['status'], 'exact_no': cert is not None, 'solve_s': ds})
            if cert:
                cert.update(axis=j, sign=s)
                certs.append(cert)
            elif v['status'] == 'Optimal':
                d = [Q(float(x)).limit_denominator(10 ** 9) for x in v['x']]
                if not rounding_envelope and check_primal('inequality', Ar, br, d).status == 'JA':
                    candidate = d
    status = 'NO_NONZERO_DIRECTION' if len(certs) == 6 else 'YES_BOUNDARY' if candidate else 'UNKNOWN_BOUNDARY_OR_RECONSTRUCTION'
    return {'status': status, 'direction': [str(x) for x in candidate] if candidate else None, 'unique_normals': len(N), 'certificates': certs, 'charts': charts, 'discovery_s': prep_s + solve_s, 'exact_check_s': check_s, 'total_s': time.perf_counter() - t0, 'min_dot_input_normals': None}

def recheck(record):
    return all((check_dual('inequality', c['A'], c['b'], c['y']).status == 'NEJ' for c in record['certificates']))
