"""Independent archived Field arithmetic, two dyadic support worlds.
A tiny local copy of one archived certificate kernel; no source tree copied.
"""
import hashlib, json
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from field_observable_snapshot import certify_observable, verify_observable
R = Path(__file__).resolve().parents[1]

def exact(C, g, T, h):
    for active in [(0, 1), (0,), (1,)]:
        E = [[C[i][j] for j in active] + [F(-1)] for i in active] + [[F(1)] * len(active) + [F(0)]]
        rhs = [h[i] - g[i] for i in active] + [T]
        readout = [[F(i == j) for j in range(len(active) + 1)] for i in range(len(active) + 1)]
        cert = certify_observable(E, rhs, readout)
        if cert['status'] != 'UNIQUE' or not verify_observable(E, rhs, readout, cert):
            continue
        x = cert['value']
        f = [F(0), F(0)]
        for (i, val) in zip(active, x[:-1]):
            f[i] = val
        slack = [sum((C[i][j] * f[j] for j in range(2))) + g[i] - h[i] - x[-1] for i in range(2)]
        if all((v >= 0 for v in f + slack)) and all((f[i] * slack[i] == 0 for i in range(2))):
            corrupt = dict(cert)
            corrupt['value'] = [v + 1 if i == 0 else v for (i, v) in enumerate(cert['value'])]
            return (f, not verify_observable(E, rhs, readout, corrupt))
    raise ArithmeticError('No exact contact branch')

def main():
    C0 = [[F(1, 64), F(0)], [F(0), F(1, 64)]]
    r = F(1, 128)
    C1 = [[v + r for v in row] for row in C0]
    gs = [[-sum((C[i][j] * F(16) for j in range(2))) for i in range(2)] for C in [C0, C1]]
    G = [C[0][0] + C[1][1] - 2 * C[0][1] for C in [C0, C1]]
    p = [(C[1][1] - C[0][1]) / gg for (C, gg) in zip([C0, C1], G)]
    summary = [np.array([float(gg), float(pp), 16.0, 16.0]) for (gg, pp) in zip(G, p)]
    rows = []
    for T in [F(24), F(32), F(40)]:
        for ht in [F(-1, 2), F(-1, 64), F(0), F(1, 64), F(1, 2)]:
            ff = []
            faults = []
            for (C, g) in zip([C0, C1], gs):
                (f, det) = exact(C, g, T, [ht, F(0)])
                ff.append(f)
                faults.append(det)
            expected = max(F(0), min(T, T * p[0] + ht / G[0]))
            assert ff[0][0] == expected and ff[1][0] == expected
            rows.append({'T_N': float(T), 'height_mm': float(ht), 'forces_N': [[float(x) for x in f] for f in ff], 'difference_N': float(abs(ff[0][0] - ff[1][0])), 'corrupt_certificate_detected': all(faults)})
    out = {'summary_identity_error': float(np.max(np.abs(summary[0] - summary[1]))), 'summary_bit_identical': summary[0].tobytes() == summary[1].tobytes(), 'downstream_max_difference_N': max((x['difference_N'] for x in rows)), 'questions': len(rows), 'closure': 'H/projected G + measured load response p sufficient in the declared linear two-support contact model', 'external_referent': {'kind': 'closed_form', 'locator': 'Exact two-support balance f1=clip(T*p+h/G,0,T), independently replayed by archived Field observable_cert.py', 'compared_quantity': 'Exact unilateral reactions after height/load changes', 'refutes_us': False}, 'field_kernel_sha256': hashlib.sha256((R / 'code/field_observable_snapshot.py').read_bytes()).hexdigest(), 'rows': rows, 'gate': all((x['difference_N'] == 0 and x['corrupt_certificate_detected'] for x in rows))}
    (R / 'raw/FIELD_EXACT_SUFFICIENCY.json').write_text(json.dumps(out, indent=2))
    print(json.dumps({k: out[k] for k in ['summary_identity_error', 'summary_bit_identical', 'downstream_max_difference_N', 'questions', 'gate']}))
if __name__ == '__main__':
    main()
