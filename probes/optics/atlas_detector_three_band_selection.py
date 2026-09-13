"""Select three-band detectors on training packets with explicit gain nuisance."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
from atlas_detector_support import candidates
from atlas_detector_third_wavelength import matrix as spectral_matrix
from atlas_detector_gain_nuisance import residual


def information(paths, masks):
    out = []
    for start in range(0, 12, 3):
        j = spectral_matrix(paths, masks, [0,1,2], tuple(range(start,start+3)))
        for i in range(3):
            out.append(j[i*3:(i+1)*3][:,[0,2+i,1]])
    return np.array(out), None


def profile(matrix, triple):
    shared = matrix[list(triple)].reshape(9,3)
    expanded = np.zeros((9,5))
    expanded[:,:2] = shared[:,[0,2]]
    for i in range(3):
        expanded[3*i:3*i+3,2+i] = shared[3*i:3*i+3,1]
    return float(np.linalg.norm(residual(expanded))), expanded


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original', required=True)
    p.add_argument('captures', nargs=2)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    rows, arrays = [], {}
    for leg, filename in enumerate(a.captures):
        with np.load(a.original, allow_pickle=False) as c:
            centers, masks = candidates(c['exits'], c['terminal'], c['source'])
            train, _ = information(c['paths'][::2], masks[::2])
        counts = masks[::2].sum(axis=0)
        distance = np.linalg.norm(centers[:, None]-centers[None, :], axis=2)
        triples = [t for t in itertools.combinations(np.flatnonzero(counts>=20).tolist(), 3) if all(distance[i,j]>10 for i,j in itertools.combinations(t,2))]
        scores = np.array([profile(train,t)[0] for t in triples])
        chosen = min(zip(triples,scores), key=lambda x:(-x[1],x[0]))[0]
        with np.load(filename, allow_pickle=False) as c:
            masks = (np.linalg.norm(c['exits'][:,None,:3]-centers[None],axis=2)<=5)&(c['terminal'][:,1]>0)[:,None]
            heldout, _ = information(c['paths'], masks)
        selected, sj = profile(heldout, chosen)
        reference, rj = profile(heldout, (0,4,8))
        rows.append(dict(selected_ids=list(chosen), feasible_triples=len(triples), selected_sensitivity=selected,
                         reference_sensitivity=reference, gain=selected/reference,
                         training_counts=counts[list(chosen)].tolist(), heldout_counts=masks[:,chosen].sum(axis=0).tolist(),
                         minimum_separation=min(float(distance[i,j]) for i,j in itertools.combinations(chosen,2)),
                         ranks=[int(np.linalg.matrix_rank(j)) for j in (sj,rj)]))
        arrays.update({f'scores_{leg}':scores, f'training_{leg}':train, f'heldout_{leg}':heldout,
                       f'selected_{leg}':sj, f'reference_{leg}':rj})
    gates = dict(full_repeat=rows[0]==rows[1] and all(arrays[f'{k}_0'].tobytes()==arrays[f'{k}_1'].tobytes() for k in ('scores','training','heldout','selected','reference')),
                 finite_rank=all(np.isfinite(v).all() for v in arrays.values()) and all(r['ranks']==[5,5] and min(r['selected_sensitivity'],r['reference_sensitivity'])>0 for r in rows),
                 constraints=all(min(r['training_counts'])>=20 and r['minimum_separation']>10 for r in rows),
                 heldout_gain=all(r['gain']>1 for r in rows))
    result = dict(rows=rows,gates={k:bool(v) for k,v in gates.items()},
                  hashes={k:hashlib.sha256(v.tobytes()).hexdigest() for k,v in arrays.items()},
                  training_triples=[list(t) for t in triples],
                  scope='Training-only equal-exposure three-band selection; reused independent-seed evaluation, wavelength-shared detector gains, no clinical or sampling certificate.')
    np.savez(a.output.with_suffix('.npz'),**arrays)
    a.output.write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(dict(rows=rows,gates=result['gates'])))
    return 0 if all(gates.values()) else 2


if __name__=='__main__':
    raise SystemExit(main())
