"""Frozen three-band placement bootstrap with algebraic nuisance projection."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from atlas_detector_support import candidates
from atlas_detector_third_wavelength import COEFFICIENTS, matrix
from atlas_oxygen_bounds import FRACTION, BACKGROUND
from atlas_detector_gain_nuisance import residual


def sensitivities(moments):
    signal=moments[...,0]
    if np.any(signal<=0):raise ValueError('Nonpositive detector signal')
    gain=np.sqrt(signal)
    s=moments[...,1]/gain;h=moments[...,2]/gain
    norm=np.sum(gain*gain,axis=-1,keepdims=True)
    s=s-gain*np.sum(gain*s,axis=-1,keepdims=True)/norm
    h=h-gain*np.sum(gain*h,axis=-1,keepdims=True)/norm
    out=[]
    for start in (0,3):
        sv=s[:,start:start+3].reshape(len(s),9)
        hv=h[:,start:start+3].reshape(len(h),9)
        hnorm=np.sum(hv*hv,axis=1,keepdims=True)
        if np.any(hnorm<=0):raise ValueError('Unresolved hemoglobin direction')
        projected=sv-hv*np.sum(sv*hv,axis=1,keepdims=True)/hnorm
        out.append(np.linalg.norm(projected,axis=1))
    return np.asarray(out)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original',required=True)
    p.add_argument('--capture',required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    with np.load(a.original,allow_pickle=False) as c:
        centers,_=candidates(c['exits'],c['terminal'],c['source'])
    with np.load(a.capture,allow_pickle=False) as c:
        paths=c['paths'];terminal=c['terminal']
        masks=(np.linalg.norm(c['exits'][:,None,:3]-centers[None],axis=2)<=5)&(terminal[:,1]>0)[:,None]
    ids=[4,7,9,0,4,8]
    # Shared detector4 is intentionally retained in both triples; resampling is paired.
    chosen_masks=masks[:,ids]
    scale=np.log(10.)*15/64500.
    blood=(scale*(COEFFICIENTS[:,0]*.8+COEFFICIENTS[:,1]*.2))[:,None]*FRACTION
    ds=(scale*(COEFFICIENTS[:,0]-COEFFICIENTS[:,1]))[:,None]*FRACTION
    depth=paths@blood.T
    weights=np.exp(-(paths@BACKGROUND)[:,None]-depth)*(2/3)
    d_s=-(paths@ds.T)*weights;d_h=-depth*weights
    features=np.stack([chosen_masks[:,:,None]*v[:,None,:] for v in (weights,d_s,d_h)],axis=-1)
    active=chosen_masks.any(axis=1);n=len(paths);k=int(active.sum())
    probability=np.r_[np.full(k,1/n),(n-k)/n]
    direct=sensitivities(features.sum(axis=0)[None])[...,0]
    qr=np.array([np.linalg.norm(residual(matrix(paths,masks,[0,1,2],t))) for t in ((4,7,9),(0,4,8))])
    difference=np.abs(direct-qr)
    algebraic=bool(np.all(difference<=np.maximum(1e-12,1e-10*np.abs(qr))))
    rows=[];arrays={}
    for leg in (0,1):
        counts=np.random.default_rng(20260914).multinomial(n,probability,size=2000)
        moments=(counts[:,:k]@features[active].reshape(k,-1)).reshape(2000,6,3,3)
        selected,reference=sensitivities(moments)
        gain=selected/reference
        for name,value in [('selected',selected),('reference',reference),('gain',gain)]:arrays[f'{name}_{leg}']=value
        rows.append(dict(percentiles=np.percentile(gain,[2.5,50,97.5]).tolist(),count_not_above_one=int(np.count_nonzero(gain<=1)),
                         counts_hash=hashlib.sha256(counts.tobytes()).hexdigest()))
    gates=dict(full_repeat=rows[0]==rows[1] and all(v.tobytes()==arrays[key[:-1]+'1'].tobytes() for key,v in arrays.items() if key.endswith('_0')),
               finite_positive=all(np.isfinite(v).all() and np.min(v)>0 for v in arrays.values()),
               algebraic=algebraic,lower_percentile=all(r['percentiles'][0]>1 for r in rows))
    result=dict(rows=rows,gates={key:bool(v) for key,v in gates.items()},draws=2000,seed=20260914,
                histories=n,active_histories=k,caps=int(np.count_nonzero(terminal[:,2])),
                input_paths_hash=hashlib.sha256(paths.tobytes()).hexdigest(),algebraic_difference=difference.tolist(),
                hashes={key:hashlib.sha256(v.tobytes()).hexdigest() for key,v in arrays.items()},
                scope='Paired empirical bootstrap of frozen three-band placement, conditional optical/gain model; no certified coverage or clinical precision.')
    np.savez(a.output.with_suffix('.npz'),**arrays)
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
    return 0 if all(gates.values()) else 2


if __name__=='__main__':raise SystemExit(main())
