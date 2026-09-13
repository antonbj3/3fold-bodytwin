"""Equal-exposure third-band information with fixed detectors and nuisance gains."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from atlas_detector_support import candidates
from atlas_oxygen_bounds import FRACTION, BACKGROUND
from atlas_detector_gain_nuisance import residual

COEFFICIENTS = np.array([[586.,1548.52],[816.,761.72],[1058.,691.32]])


def matrix(paths, masks, bands, triple):
    extinction = COEFFICIENTS[bands]
    scale = np.log(10.)*15/64500.
    blood = (scale*(extinction[:,0]*.8+extinction[:,1]*.2))[:,None]*FRACTION
    ds = (scale*(extinction[:,0]-extinction[:,1]))[:,None]*FRACTION
    depth = paths@blood.T
    weights = np.exp(-(paths@BACKGROUND)[:,None]-depth)*(2/len(bands))
    sat = -(paths@ds.T)*weights
    hb = -depth*weights
    out = np.zeros((3*len(bands),5))
    for detector, identity in enumerate(triple):
        mask = masks[:,identity]
        signal = weights[mask].sum(axis=0)
        if np.any(signal<=0): raise ValueError('Positive support required')
        block = slice(detector*len(bands),(detector+1)*len(bands))
        out[block,0] = sat[mask].sum(axis=0)/np.sqrt(signal)
        out[block,1] = hb[mask].sum(axis=0)/np.sqrt(signal)
        out[block,2+detector] = np.sqrt(signal)
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original',required=True)
    p.add_argument('captures',nargs=2)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    with np.load(a.original,allow_pickle=False) as c:
        centers,_=candidates(c['exits'],c['terminal'],c['source'])
    legs, arrays=[],{}
    for leg,filename in enumerate(a.captures):
        with np.load(filename,allow_pickle=False) as c:
            paths=c['paths']
            masks=(np.linalg.norm(c['exits'][:,None,:3]-centers[None],axis=2)<=5)&(c['terminal'][:,1]>0)[:,None]
        rows=[]
        for triple in ((5,7,9),(2,6,10),(0,4,8)):
            scores=[]; ranks=[]
            for name,bands in [('two',[0,2]),('three',[0,1,2]),('duplicate',[0,2,0,2])]:
                j=matrix(paths,masks,bands,triple)
                spectrum=np.linalg.svd(j,compute_uv=False)
                projected=residual(j)
                scores.append(float(np.linalg.norm(projected)));ranks.append(int(np.linalg.matrix_rank(j)))
                key=name+'_'+''.join(map(str,triple))
                arrays[f'{key}_jacobian_{leg}']=j
                arrays[f'{key}_spectrum_{leg}']=spectrum
                arrays[f'{key}_profile_{leg}']=projected
            rows.append(dict(ids=list(triple),two_sensitivity=scores[0],three_sensitivity=scores[1],
                             ratio=scores[1]/scores[0],duplicate_difference=abs(scores[2]-scores[0]),ranks=ranks))
        legs.append(rows)
    gates=dict(full_repeat=legs[0]==legs[1] and all(v.tobytes()==arrays[k[:-1]+'1'].tobytes() for k,v in arrays.items() if k.endswith('_0')),
               finite_rank=all(np.isfinite(v).all() for v in arrays.values()) and all(r['ranks']==[5,5,5] for rows in legs for r in rows),
               improved=all(r['ratio']>1 for rows in legs for r in rows),
               duplicate_control=all(r['duplicate_difference']<=max(1e-12,1e-10*r['two_sensitivity']) for rows in legs for r in rows))
    result=dict(legs=legs,gates={k:bool(v) for k,v in gates.items()},coefficients=COEFFICIENTS.tolist(),
                source='https://doi.org/10.1184/R1/24530353.v1',
                hashes={k:hashlib.sha256(v.tobytes()).hexdigest() for k,v in arrays.items()},
                scope='Equal total exposure, frozen reweighted paths and assumed wavelength-independent scattering; not an experimental spectral optimum.')
    np.savez(a.output.with_suffix('.npz'),**arrays)
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(legs=legs,gates=result['gates'])))
    return 0 if all(gates.values()) else 2


if __name__=='__main__':
    raise SystemExit(main())
