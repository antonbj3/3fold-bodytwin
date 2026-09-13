"""Fixed detector information over a preregistered optical parameter grid."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from atlas_detector_support import candidates
from atlas_oxygen_bounds import EXTINCTION, FRACTION, BACKGROUND


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original', required=True)
    p.add_argument('captures', nargs=2)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    with np.load(a.original, allow_pickle=False) as c:
        centers, _ = candidates(c['exits'], c['terminal'], c['source'])
    legs, arrays = [], {}
    for leg, filename in enumerate(a.captures):
        with np.load(filename, allow_pickle=False) as c:
            paths = c['paths']
            masks = (np.linalg.norm(c['exits'][:, None, :3]-centers[None, [5,7,9,0,4,8]], axis=2)<=5) & (c['terminal'][:, 1]>0)[:, None]
        rows, matrices, singular = [], [], []
        for sat in (.7, .8, .9, 1.):
            for hb in (12., 15., 18.):
                scale = np.log(10.)*hb/64500.
                blood = (scale*(EXTINCTION[:, 0]*sat+EXTINCTION[:, 1]*(1-sat)))[:, None]*FRACTION
                ds = (scale*(EXTINCTION[:, 0]-EXTINCTION[:, 1]))[:, None]*FRACTION
                depth = paths@blood.T
                w = np.exp(-(paths@BACKGROUND)[:, None]-depth)
                derivatives = [-(paths@ds.T)*w, w, -depth*w]
                matrix = np.array([np.column_stack([d[m].sum(axis=0) for d in derivatives])/np.sqrt(w[m].sum(axis=0))[:, None] for m in masks.T])
                spectra = np.array([np.linalg.svd(matrix[start:start+3].reshape(6,3), compute_uv=False) for start in (0,3)])
                rows.append(dict(saturation=sat, hemoglobin=hb, gain=float(spectra[0,-1]/spectra[1,-1]),
                                 full_rank=all(np.linalg.matrix_rank(matrix[start:start+3].reshape(6,3))==3 for start in (0,3)),
                                 selected_singular_values=spectra[0].tolist(), reference_singular_values=spectra[1].tolist()))
                matrices.append(matrix); singular.append(spectra)
        legs.append(rows)
        arrays[f'jacobians_{leg}'] = np.array(matrices)
        arrays[f'singular_values_{leg}'] = np.array(singular)
    gates = dict(full_repeat=legs[0]==legs[1] and all(arrays[f'{k}_0'].tobytes()==arrays[f'{k}_1'].tobytes() for k in ('jacobians','singular_values')),
                 finite_rank=all(np.isfinite(v).all() for v in arrays.values()) and all(r['full_rank'] for rows in legs for r in rows),
                 all_gains=all(r['gain']>1 for rows in legs for r in rows))
    result = dict(legs=legs, gates={k:bool(v) for k,v in gates.items()},
                  hashes={k:hashlib.sha256(v.tobytes()).hexdigest() for k,v in arrays.items()},
                  scope='Finite assumed-optics grid at frozen centers; no continuous robustness, physiological or clinical certificate.')
    np.savez(a.output.with_suffix('.npz'), **arrays)
    a.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))
    return 0 if all(gates.values()) else 2


if __name__ == '__main__':
    raise SystemExit(main())
