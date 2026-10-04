import json, time, os
from pathlib import Path
import numpy as np
from scipy.linalg import solve_banded
from scipy.special import erfc
ROOT = Path(os.environ.get('X12_RUN_ROOT', str(Path(__file__).resolve().parents[1])))

def solve(alpha):
    dx = 0.025
    dt = 0.001
    x = np.arange(0, 20 + dx / 2, dx)
    n = len(x) - 2
    r = alpha * dt / dx ** 2
    ab = np.zeros((3, n))
    ab[0, 1:] = -r
    ab[1, :] = 1 + 2 * r
    ab[2, :-1] = -r
    u = np.zeros(n)
    out = {}
    for step in range(1, 10001):
        rhs = u.copy()
        rhs[0] += r
        u = solve_banded((1, 1), ab, rhs, overwrite_b=True, check_finite=False)
        if step in [1000, 5000, 10000]:
            out[step * dt] = np.r_[1, u, 0]
    return (x, out)

def main():
    t0 = time.time()
    (x, y) = solve(0.2)
    (_, bad) = solve(0.6)
    checks = []
    for (t, v) in y.items():
        exact = erfc(x / (2 * np.sqrt(0.2 * t)))
        err = float(np.max(abs(v - exact)))
        poison = float(np.max(abs(bad[t] - exact)))
        checks.append({'time_s': t, 'max_abs_T_fraction_error': err, 'gate': err <= 0.002, 'poisoned_alpha_error': poison, 'poison_rejected': poison > 0.002})
    result = {'checks': checks, 'all_pass': all((r['gate'] for r in checks)), 'all_poison_rejected': all((r['poison_rejected'] for r in checks)), 'external_referent': {'kind': 'closed_form', 'locator': 'https://web.mit.edu/snively/www/2_005%20Coursenotes-1.pdf', 'compared_quantity': 'normalizedtemperature halfspace step, Chapter6 Eq6.206', 'refutes_us': False}, 'wall_s': time.time() - t0, 'source_amplitude': 'UNKNOWN', 'clinical_temperature': 'UNKNOWN'}
    (ROOT / 'raw/THERMAL_CONTROL.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))
    assert result['all_pass'] and result['all_poison_rejected']
if __name__ == '__main__':
    main()
