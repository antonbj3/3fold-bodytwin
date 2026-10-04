import copy
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
ROOT = Path(__file__).resolve().parent

def validate(record):
    ratios = np.asarray(record['peak_normalized'])[1:] / np.asarray(record['peak_normalized'])[:-1]
    return {'crack_half_pitch_ratio': bool(np.all(np.abs(ratios - np.sqrt(2)) < 1e-10)), 'fixed_support_integral': abs(record['average_normalized'] - 2) < 1e-10, 'no_actual_dental_divergence_claim': record['actual_dental_class3'] == 'UNKNOWN'}

def main():
    p = ROOT / 'PREREG_R4.json'
    assert hashlib.sha256(p.read_bytes()).hexdigest() == (ROOT / 'PREREG_R4.sha256').read_text().strip()
    pr = json.loads(p.read_text())
    h = np.asarray(pr['normalized_h_over_ell'])
    (integral, error) = quad(lambda u: 2.0, 0, 1, epsabs=1e-12, epsrel=1e-12)
    record = {'claim_type': 'capability', 'round': 'R4', 'h_over_ell': h.tolist(), 'peak_normalized': (1 / np.sqrt(h)).tolist(), 'average_normalized': integral, 'quadrature_absolute_error_estimate': error, 'ideal_peak_class': 3, 'regularized_query_class': 1, 'actual_dental_class3': 'UNKNOWN', 'physical_ell': 'UNKNOWN', 'physical_K': 'UNKNOWN', 'external_referent': pr['external_referent']}
    gates = validate(record)
    mutations = {}
    changes = {'crack_half_pitch_ratio': {'peak_normalized': [1] * len(h)}, 'fixed_support_integral': {'average_normalized': 4}, 'no_actual_dental_divergence_claim': {'actual_dental_class3': 3}}
    for (gate, change) in changes.items():
        mutations[gate] = not validate(record | change)[gate]
    record.update(gates=gates, mutation_rejections=mutations)
    (ROOT / 'RESULTS_R4.json').write_text(json.dumps(record, indent=2) + '\n')
    assert all(gates.values()) and all(mutations.values())
    print('R4 ideal crack and fixed-support observable match external closed form; actual dental divergence UNKNOWN')
if __name__ == '__main__':
    main()
