"""Evaluate previously frozen three-band placement on the preregistered new seed."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from atlas_detector_support import candidates
from atlas_detector_third_wavelength import matrix
from atlas_detector_gain_nuisance import residual

ATLAS = '8f6c71ced645db796de1736ee77eb7514cdf6e7c6315407ecb2bebf81ee01672'
CENTER = '10fa3ea6c18c36c78c2c13c0a13245c28ca2d042f98400dbd2930f21e7a857bd'
TRANSPORT = '62e720ec695c3e6e7272e12e7db82603a1e5663ceaf72d0c9bd5889b1e93b5bc'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original',required=True)
    p.add_argument('--receipt',type=Path,required=True)
    p.add_argument('captures',nargs=2)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    receipt=json.loads(a.receipt.read_text())
    if receipt['seed']!=20260916 or receipt['atlas_sha256']!=ATLAS or receipt['source_hashes']['transport']!=TRANSPORT:
        raise ValueError('Confirmation receipt does not match frozen experiment')
    with np.load(a.original,allow_pickle=False) as c:
        centers,training_masks=candidates(c['exits'],c['terminal'],c['source'])
    if hashlib.sha256(centers.tobytes()).hexdigest()!=CENTER:
        raise ValueError('Frozen detector center identity mismatch')
    triples=((4,7,9),(0,4,8))
    distances=np.linalg.norm(centers[:,None]-centers[None],axis=2)
    rows, arrays=[],{}
    for leg,filename in enumerate(a.captures):
        with np.load(filename,allow_pickle=False) as c:
            raw={k:c[k] for k in c.files if k!='source'}
        if set(raw)!={'paths','exits','terminal','counters','absorption'}:
            raise ValueError('Unexpected transport array inventory')
        hashes={k:hashlib.sha256(v.tobytes()).hexdigest() for k,v in raw.items()}
        if hashes!=receipt['rows'][leg]['hashes']:
            raise ValueError('Full transport hashes do not match receipt')
        terminal=raw['terminal'];paths=raw['paths']
        if len(paths)!=100000:raise ValueError('Unexpected packet count')
        masks=(np.linalg.norm(raw['exits'][:,None,:3]-centers[None],axis=2)<=5)&(terminal[:,1]>0)[:,None]
        scores=[];ranks=[]
        for name,triple in zip(('selected','reference'),triples):
            j=matrix(paths,masks,[0,1,2],triple)
            spectrum=np.linalg.svd(j,compute_uv=False)
            projected=residual(j)
            scores.append(float(np.linalg.norm(projected)));ranks.append(int(np.linalg.matrix_rank(j)))
            arrays[f'{name}_jacobian_{leg}']=j
            arrays[f'{name}_spectrum_{leg}']=spectrum
            arrays[f'{name}_profile_{leg}']=projected
        sums=terminal.sum(axis=0,dtype=np.int64)
        rows.append(dict(hashes=hashes,selected_sensitivity=scores[0],reference_sensitivity=scores[1],gain=scores[0]/scores[1],
                         ranks=ranks,counts=[masks[:,t].sum(axis=0).tolist() for t in triples],
                         absorbed=int(sums[0]),escaped=int(sums[1]),residual=int(sums[2]),leaks=int(sums[3]),caps=int(raw['counters'][2])))
    constraints=all(min(training_masks[::2,t].sum(axis=0))>=20 and min(distances[i,j] for i in t for j in t if i!=j)>10 for t in triples)
    gates=dict(full_repeat=rows[0]==rows[1] and all(v.tobytes()==arrays[k[:-1]+'1'].tobytes() for k,v in arrays.items() if k.endswith('_0')),
               energy_leaks=all(r['absorbed']+r['escaped']+r['residual']==100000*2**30 and r['leaks']==0 for r in rows),
               finite_rank=all(np.isfinite(v).all() for v in arrays.values()) and all(r['ranks']==[5,5] and min(r['selected_sensitivity'],r['reference_sensitivity'])>0 for r in rows),
               frozen_constraints=constraints,
               fresh_gain=all(r['gain']>1 for r in rows))
    result=dict(rows=rows,gates={k:bool(v) for k,v in gates.items()},seed=20260916,center_hash=CENTER,selected_ids=[4,7,9],reference_ids=[0,4,8],
                receipt_sha256=hashlib.sha256(a.receipt.read_bytes()).hexdigest(),
                hashes={k:hashlib.sha256(v.tobytes()).hexdigest() for k,v in arrays.items()},
                scope='Frozen three-band placement on preregistered fresh-seed paths; conditional optical/gain model, capped proposals retained, no clinical certificate.')
    np.savez(a.output.with_suffix('.npz'),**arrays)
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
    return 0 if all(gates.values()) else 2


if __name__=='__main__':
    raise SystemExit(main())
