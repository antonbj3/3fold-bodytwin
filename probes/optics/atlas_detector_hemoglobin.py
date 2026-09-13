"""Observe frozen detector information with unknown log hemoglobin concentration."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from atlas_detector_support import candidates
from atlas_oxygen_bounds import EXTINCTION, FRACTION, BACKGROUND


def information(paths, masks):
    scale = np.log(10.) * 15 / 64500.
    blood = (scale * (EXTINCTION[:, 0]*.8 + EXTINCTION[:, 1]*.2))[:, None]*FRACTION
    ds = (scale * (EXTINCTION[:, 0]-EXTINCTION[:, 1]))[:, None]*FRACTION
    depth = paths @ blood.T
    weights = np.exp(-(paths @ BACKGROUND)[:, None]-depth)
    saturation = -(paths @ ds.T)*weights
    hemoglobin = -depth*weights
    step = 1e-5
    perturbed = [np.exp(-(paths @ BACKGROUND)[:, None]-depth*np.exp(t)) for t in (-step, step)]
    matrices, differences = [], []
    for mask in masks.T:
        signal = weights[mask].sum(axis=0)
        if np.any(signal <= 0):
            raise ValueError('Positive detector support required')
        dh = hemoglobin[mask].sum(axis=0)
        fd = (perturbed[1][mask].sum(axis=0)-perturbed[0][mask].sum(axis=0))/(2*step)
        differences.append(float(np.max(np.abs(dh-fd)/np.maximum(np.abs(dh), np.finfo(float).tiny))))
        matrices.append(np.column_stack((saturation[mask].sum(axis=0), signal, dh))/np.sqrt(signal)[:, None])
    return np.asarray(matrices), max(differences)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original', required=True)
    p.add_argument('captures', nargs=2)
    p.add_argument('--output', required=True, type=Path)
    a = p.parse_args()
    with np.load(a.original, allow_pickle=False) as source:
        centers, _ = candidates(source['exits'], source['terminal'], source['source'])
    ids = [5, 7, 9, 0, 4, 8]
    rows, arrays = [], {}
    for leg, filename in enumerate(a.captures):
        with np.load(filename, allow_pickle=False) as capture:
            paths, exits, terminal = (capture[k] for k in ('paths', 'exits', 'terminal'))
            masks = (np.linalg.norm(exits[:, None, :3]-centers[None, ids], axis=2) <= 5) & (terminal[:, 1] > 0)[:, None]
            matrix, fd_error = information(paths, masks)
            selected = np.linalg.svd(matrix[:3].reshape(6, 3), compute_uv=False)
            reference = np.linalg.svd(matrix[3:].reshape(6, 3), compute_uv=False)
            rows.append(dict(selected_singular_values=selected.tolist(), reference_singular_values=reference.tolist(),
                             gain=float(selected[-1]/reference[-1]), finite_difference_relative_error=fd_error,
                             individual_nullity=[3-int(np.linalg.matrix_rank(m)) for m in matrix],
                             counts=masks.sum(axis=0).tolist(), caps=int(np.count_nonzero(terminal[:, 2])),
                             path_hash=hashlib.sha256(paths.tobytes()).hexdigest()))
            arrays.update({f'jacobian_{leg}': matrix, f'selected_svd_{leg}': selected, f'reference_svd_{leg}': reference})
    gates = dict(full_repeat=rows[0] == rows[1] and all(arrays[f'{k}_0'].tobytes() == arrays[f'{k}_1'].tobytes() for k in ('jacobian', 'selected_svd', 'reference_svd')),
                 derivative=all(r['finite_difference_relative_error'] <= 1e-7 for r in rows),
                 positive_rank=all(np.isfinite(arrays[f'jacobian_{leg}']).all() and all(np.linalg.matrix_rank(arrays[f'jacobian_{leg}'][start:start+3].reshape(6, 3)) == 3 for start in (0, 3)) for leg in (0, 1)),
                 selected_gain=all(r['gain'] > 1 for r in rows))
    result = dict(rows=rows, gates={k: bool(v) for k, v in gates.items()}, detector_ids=ids,
                  center_hash=hashlib.sha256(centers.tobytes()).hexdigest(),
                  hashes={k: hashlib.sha256(v.tobytes()).hexdigest() for k, v in arrays.items()},
                  scope='Fixed centers and unit parameter scales, conditional Poisson information with unknown log-Hb; no detector reselection, clinical calibration or capped-path acceptance.')
    np.savez(a.output.with_suffix('.npz'), **arrays)
    a.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))
    return 0 if all(gates.values()) else 2


if __name__ == '__main__':
    raise SystemExit(main())
