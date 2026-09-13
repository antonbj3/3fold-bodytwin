"""Observe the pinned reference spectral objective, retaining its zero singularity."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from pulse_reference_mechanism import PINS

SOURCE_PINS = dict(PINS, **{
    'generate_analytical_L.m': 'ce31724e4895a56bdeb7e998de2b0e4dc798c81c903df756298face86054ce1d',
    'generate_measured_L.m': '8f18326ec9fb94965d23498a3bf3b89f12dffe5a435c981d56cd04e45bdca98e'})


def observe(data, coeff):
    wavelengths = data['wavelengths_all'].ravel()
    indices = [int(np.flatnonzero(coeff['wavelength'].ravel() == w)[0]) for w in wavelengths]
    eo = coeff['ext_hbo'].ravel()[indices] * 2.3026 / 1000
    ed = coeff['ext_hb'].ravel()[indices] * 2.3026 / 1000
    grid = np.arange(101) / 100
    mixture = eo[:, None] * grid + ed[:, None] * (1-grid)
    mua = .05 * mixture
    musp = (260.6724 * wavelengths**(-.4668))[:, None]
    intensity = data['I_r_all']
    od = np.log(intensity[..., 0] / intensity[..., 1])
    rows, literal_all, limit_all = [], [], []
    maximum_difference = 0.
    for distance_index, distance in enumerate(data['r_all'].ravel()):
        r = distance / 10
        root = r * np.sqrt(3*mua*musp)
        length = np.sqrt(3*musp)*root/(2*np.sqrt(mua)*(root+1))*r
        analytical = length / length[:1]
        for saturation_index, saturation in enumerate(data['oxy_sim_all'].ravel()):
            relative_od = od[:, saturation_index, distance_index] / od[0, saturation_index, distance_index]
            with np.errstate(divide='ignore', invalid='ignore'):
                k = eo[:, None] + ed[:, None]*(1-grid)/grid
                literal = relative_od[:, None] / (k/k[:1])
            literal[0] = 1
            limiting = relative_od[:, None]*mixture[:1]/mixture
            maximum_difference = max(maximum_difference, float(np.max(np.abs(literal[:, 1:]-limiting[:, 1:]))))
            literal_objective = np.sum((analytical-literal)**2, axis=0)
            limit_objective = np.sum((analytical-limiting)**2, axis=0)
            finite_indices = np.flatnonzero(np.isfinite(literal_objective))
            finite_best = int(finite_indices[np.argmin(literal_objective[finite_indices])])
            limit_best = int(np.argmin(limit_objective))
            rows.append(dict(distance=float(distance), saturation=float(saturation),
                             literal_nonfinite=int(np.sum(~np.isfinite(literal_objective))),
                             literal_finite_index=finite_best, limiting_index=limit_best,
                             index_error=abs(limit_best-int(round(float(saturation)*100)))))
            literal_all.append(literal_objective)
            limit_all.append(limit_objective)
    return rows, np.asarray(literal_all), np.asarray(limit_all), maximum_difference


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source_dir', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    for name, expected in SOURCE_PINS.items():
        if hashlib.sha256((args.source_dir/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Source hash mismatch')
    data = loadmat(args.source_dir/'demo_mc_data.mat')
    coeff = loadmat(args.source_dir/'ext_coef.mat')
    legs = [observe(data, coeff) for _ in range(2)]
    rows, literal, limiting, difference = legs[0]
    gates = dict(full_repeat=rows == legs[1][0] and literal.tobytes() == legs[1][1].tobytes() and limiting.tobytes() == legs[1][2].tobytes(),
                 positive_grid_agreement=difference <= 1e-12,
                 literal_grid_finite=bool(np.isfinite(literal).all()),
                 demo_distance_index_error=max(row['index_error'] for row in rows if row['distance'] == 30) <= 1)
    arrays = {f'{name}_{leg}': legs[leg][index] for leg in range(2) for name, index in [('literal', 1), ('limiting', 2)]}
    np.savez(args.output.with_suffix('.npz'), **arrays)
    report = dict(rows=rows, gates=gates, maximum_positive_grid_difference=difference,
                  hashes={name: hashlib.sha256(array.tobytes()).hexdigest() for name, array in arrays.items()},
                  source_pins=SOURCE_PINS, source='https://doi.org/10.1184/R1/24530353.v1',
                  scope='Independent expression observer on simulated inputs; finite-only argmin is a diagnostic, not a claim about MATLAB NaN handling or clinical calibration.')
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(gates=gates, maximum_positive_grid_difference=difference, rows=rows)))
    return 0 if all(gates.values()) else 2


if __name__ == '__main__':
    raise SystemExit(main())
