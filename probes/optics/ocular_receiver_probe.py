"""Check incident-power normalization against frozen packet ledgers."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from bodytwin.geometry.optics.ocular_receiver_v1 import normalize


def main():
    source = ROOT/'reports/ocular_transport_1000000.json'
    transport = json.loads(source.read_text())
    rows, arrays = [], {}
    with np.load(ROOT/'reports/ocular_transport_1000000_arrays.npz') as archive:
        for r in transport['rows']:
            key = f"{r['wavelength']}_{r['leg']}"
            density, area, overflow, total = normalize(
                archive[key+'_retinal_histogram'], transport['retinal_bin_edges'],
                r['photons']*(1 << 30), r['histogram_overflow_energy'])
            error = abs(float(np.sum(density*area))+overflow-total)
            arrays[key] = density
            rows.append(dict(wavelength=r['wavelength'], leg=r['leg'],
                             retinal_fraction=total, overflow_fraction=overflow,
                             integral_error=error, sha256=hashlib.sha256(density.tobytes()).hexdigest()))
    rejected = 0
    for args in (([[-1]], [0, 1], 1, 0), ([[0]], [1, 0], 1, 0),
                 ([[2]], [0, 1], 1, 0), ([[0]], [0, 1], 1, 2)):
        try:
            normalize(*args)
        except ValueError:
            rejected += 1
    zero, area, overflow, total = normalize(np.zeros((1, 1), np.int64), [0, 1], 1)
    gates = dict(integral=all(r['integral_error'] <= 1e-12 for r in rows) and
                 float(np.sum(zero*area))+overflow == total == 0,
                 repeat=all(arrays[f'{w}_0'].tobytes() == arrays[f'{w}_1'].tobytes() for w in (532, 650)),
                 invalid_rejected=rejected == 4)
    report = dict(rows=rows, gates=gates, rejected_controls=rejected,
                  input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  scope='Unit incident power per synthetic wavelength; explicit map overflow; no physical source power or clinical dose.')
    (ROOT/'reports/ocular_receiver.json').write_text(json.dumps(report, indent=2)+'\n')
    np.savez_compressed(ROOT/'reports/ocular_receiver_arrays.npz', **arrays)
    print(json.dumps(report))
    return 0 if all(gates.values()) else 2


if __name__ == '__main__':
    raise SystemExit(main())
