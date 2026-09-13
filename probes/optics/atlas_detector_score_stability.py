"""Measure frozen candidate score instability and independent contrast algebra."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from atlas_detector_support import candidates
from atlas_detector_hemoglobin import information
from atlas_detector_nuisance_selection import profile


def contrast_score(matrix, triple):
    local = matrix[list(triple)]
    gain = local[:, :, 1]
    contrast = np.column_stack((gain[:, 1], -gain[:, 0]))/np.linalg.norm(gain, axis=1)[:, None]
    s = np.sum(contrast*local[:, :, 0], axis=1)
    h = np.sum(contrast*local[:, :, 2], axis=1)
    remaining = s-h*(h@s)/(h@h)
    return float(np.linalg.norm(remaining)), float(np.linalg.norm(remaining)/np.linalg.norm(s))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original', required=True)
    p.add_argument('--selection', required=True, type=Path)
    p.add_argument('captures', nargs=2)
    p.add_argument('--output', required=True, type=Path)
    a = p.parse_args()
    specification = json.loads(a.selection.read_text())
    triples = [tuple(t) for t in specification['training_triples']]
    chosen, reference = triples.index((2,6,10)), triples.index((0,4,8))
    rows, arrays = [], {}
    for leg, filename in enumerate(a.captures):
        with np.load(a.original, allow_pickle=False) as c:
            centers, masks = candidates(c['exits'], c['terminal'], c['source'])
            train, _ = information(c['paths'][::2], masks[::2])
        with np.load(filename, allow_pickle=False) as c:
            paths = c['paths']
            masks = (np.linalg.norm(c['exits'][:, None, :3]-centers[None], axis=2)<=5)&(c['terminal'][:,1]>0)[:,None]
        matrices = [train, information(paths,masks)[0]]
        matrices += [information(paths[b::10],masks[b::10])[0] for b in range(10)]
        scores = np.array([[profile(m,t)[0] for t in triples] for m in matrices])
        contrasts = np.array([[contrast_score(m,t) for t in triples] for m in matrices])
        allowed = np.maximum(1e-12,1e-10*np.abs(scores))
        order = np.argsort(-scores,axis=1,kind='stable')
        ranks = np.argsort(order,axis=1,kind='stable')
        gains = scores[2:,chosen]/scores[2:,reference]
        rows.append(dict(block_gains=gains.tolist(),block_winners=[list(triples[i]) for i in order[2:,0]],
                         original_winner_rank_full=int(ranks[1,chosen])+1,
                         full_winner=list(triples[order[1,0]]),
                         training_full_rank_correlation=float(np.corrcoef(ranks[:2])[0,1]),
                         full_chosen_sin_angle=float(contrasts[1,chosen,1]),
                         full_reference_sin_angle=float(contrasts[1,reference,1]),
                         max_contrast_difference=float(np.max(np.abs(scores-contrasts[:,:,0]))),
                         contrast_agreement=bool(np.all(np.abs(scores-contrasts[:,:,0])<=allowed))))
        for name,value in [('matrices',matrices),('scores',scores),('contrasts',contrasts),('ranks',ranks)]:
            arrays[f'{name}_{leg}']=np.asarray(value)
    gates=dict(full_repeat=rows[0]==rows[1] and all(arrays[f'{k}_0'].tobytes()==arrays[f'{k}_1'].tobytes() for k in ('matrices','scores','contrasts','ranks')),
               finite_nonnegative=all(np.isfinite(v).all() for v in arrays.values()) and all(np.min(arrays[f'scores_{leg}'])>=0 for leg in (0,1)),
               independent_contrast=all(r['contrast_agreement'] for r in rows),
               fixed_block_gain=all(min(r['block_gains'])>1 for r in rows))
    result=dict(rows=rows,gates={k:bool(v) for k,v in gates.items()},triples=[list(t) for t in triples],
                hashes={k:hashlib.sha256(v.tobytes()).hexdigest() for k,v in arrays.items()},
                scope='Deterministic packet-index blocks and reused-data rank diagnostic; no detector promotion or confirmatory sampling claim.')
    np.savez(a.output.with_suffix('.npz'),**arrays)
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(rows=rows,gates=result['gates'])))
    return 0 if all(gates.values()) else 2


if __name__=='__main__':
    raise SystemExit(main())
