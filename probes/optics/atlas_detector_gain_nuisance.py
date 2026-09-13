"""Profile saturation information with independent detector gains."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def residual(matrix):
    nuisance = matrix[:, 1:]
    q, _ = np.linalg.qr(nuisance, mode='reduced')
    return matrix[:, 0]-q@(q.T@matrix[:, 0])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid_arrays', type=Path)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    legs, arrays = [], {}
    with np.load(a.grid_arrays, allow_pickle=False) as source:
        for leg in (0, 1):
            rows, expanded, spectra, profiles = [], [], [], []
            for index, matrices in enumerate(source[f'jacobians_{leg}']):
                pair = []
                for start in (0, 3):
                    shared = matrices[start:start+3].reshape(6, 3)
                    independent = np.zeros((6, 5))
                    independent[:, :2] = shared[:, [0, 2]]
                    for detector in range(3):
                        independent[2*detector:2*detector+2, detector+2] = shared[2*detector:2*detector+2, 1]
                    values = np.linalg.svd(independent, compute_uv=False)
                    shared_residual = residual(shared)
                    independent_residual = residual(independent)
                    # A separate unknown gain for every channel spans the observation space.
                    channel = np.column_stack((shared[:, 0], np.diag(shared[:, 1])))
                    channel_residual = residual(channel)
                    pair.append(dict(shared_sensitivity=float(np.linalg.norm(shared_residual)),
                                     independent_sensitivity=float(np.linalg.norm(independent_residual)),
                                     rank=int(np.linalg.matrix_rank(independent)), singular_values=values.tolist(),
                                     channel_relative_residual=float(np.linalg.norm(channel_residual)/np.linalg.norm(shared[:, 0]))))
                    expanded.append(independent); spectra.append(values)
                    profiles.append(np.array([shared_residual, independent_residual, channel_residual]))
                rows.append(dict(grid_index=index, selected=pair[0], reference=pair[1],
                                 gain=pair[0]['independent_sensitivity']/pair[1]['independent_sensitivity']))
            legs.append(rows)
            for name, value in [('expanded', expanded), ('spectra', spectra), ('profiles', profiles)]:
                arrays[f'{name}_{leg}'] = np.asarray(value)
    gates = dict(full_repeat=legs[0]==legs[1] and all(arrays[f'{k}_0'].tobytes()==arrays[f'{k}_1'].tobytes() for k in ('expanded','spectra','profiles')),
                 full_rank=all(np.isfinite(v).all() for v in arrays.values()) and all(r[arm]['rank']==5 for rows in legs for r in rows for arm in ('selected','reference')),
                 profiled_gain=all(r['gain']>1 for rows in legs for r in rows),
                 nuisance_monotonic=all(r[arm]['independent_sensitivity']<=r[arm]['shared_sensitivity']*(1+1e-12) for rows in legs for r in rows for arm in ('selected','reference')),
                 channel_null=all(r[arm]['channel_relative_residual']<=1e-12 for rows in legs for r in rows for arm in ('selected','reference')))
    result = dict(legs=legs, gates={k:bool(v) for k,v in gates.items()},
                  source_sha256=hashlib.sha256(a.grid_arrays.read_bytes()).hexdigest(),
                  hashes={k:hashlib.sha256(v.tobytes()).hexdigest() for k,v in arrays.items()},
                  scope='Finite-grid local Poisson information; independent detector gains shared across wavelengths, plus an unidentifiable per-channel-gain control. No clinical precision certificate.')
    a.output.write_text(json.dumps(result, indent=2)+'\n')
    np.savez(a.output.with_suffix('.npz'), **arrays)
    print(json.dumps(dict(gates=result['gates'], gains=[r['gain'] for r in legs[0]])))
    return 0 if all(gates.values()) else 2


if __name__ == '__main__':
    raise SystemExit(main())
