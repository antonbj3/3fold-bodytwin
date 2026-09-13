"""Paired packet bootstrap for fixed detector information gain."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from atlas_detector_support import candidates
from atlas_oxygen_bounds import EXTINCTION,FRACTION,BACKGROUND


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original',required=True);p.add_argument('--capture',required=True)
    p.add_argument('--output',required=True,type=Path);a=p.parse_args()
    with np.load(a.original,allow_pickle=False) as c:centers,_=candidates(c['exits'],c['terminal'],c['source'])
    ids=[5,7,9,0,4,8]
    with np.load(a.capture,allow_pickle=False) as c:
        paths=c['paths'];escaped=c['terminal'][:,1]>0;exits=c['exits']
    masks=(np.linalg.norm(exits[:,None,:3]-centers[None,ids,:],axis=2)<=5)&escaped[:,None]
    scale=np.log(10.)*15/64500.
    mu=BACKGROUND[None,:]+(scale*(EXTINCTION[:,0]*.8+EXTINCTION[:,1]*.2))[:,None]*FRACTION
    dmu=(scale*(EXTINCTION[:,0]-EXTINCTION[:,1]))[:,None]*FRACTION
    weights=np.exp(-(paths@mu.T));derivative=-(paths@dmu.T)*weights
    features=np.stack((masks[:,:,None]*weights[:,None,:],masks[:,:,None]*derivative[:,None,:]),axis=-1)
    active=masks.any(axis=1);n=len(paths);k=int(active.sum())
    # Exactly aggregate identical zero-contribution histories into one category.
    probability=np.r_[np.full(k,1/n),(n-k)/n]
    rows=[];arrays={}
    for leg in (0,1):
        counts=np.random.default_rng(20260914).multinomial(n,probability,size=2000)
        moments=(counts[:,:k]@features[active].reshape(k,-1)).reshape(2000,6,2,2)
        signal=moments[:,:,:,0];gradient=moments[:,:,:,1]
        if np.any(signal<=0):raise ValueError('Bootstrap produced zero detector signal')
        J=np.stack((gradient/np.sqrt(signal),np.sqrt(signal)),axis=-1)
        selected=np.linalg.svd(J[:,:3].reshape(2000,6,2),compute_uv=False)[:,-1]
        reference=np.linalg.svd(J[:,3:].reshape(2000,6,2),compute_uv=False)[:,-1]
        gain=selected/reference
        arrays[f'gain_{leg}']=gain;arrays[f'selected_{leg}']=selected;arrays[f'reference_{leg}']=reference
        rows.append(dict(percentiles=np.percentile(gain,[2.5,50,97.5]).tolist(),fraction_not_above_one=float(np.mean(gain<=1)),
                         count_not_above_one=int(np.count_nonzero(gain<=1)),gain_hash=hashlib.sha256(gain.tobytes()).hexdigest(),
                         counts_hash=hashlib.sha256(counts.tobytes()).hexdigest(),finite_positive=bool(np.isfinite(gain).all() and np.min(gain)>0)))
    gates=dict(full_repeat=rows[0]==rows[1] and all(arrays[k].tobytes()==arrays[k[:-1]+'1'].tobytes() for k in arrays if k.endswith('_0')),
               finite_positive=all(r['finite_positive'] for r in rows),lower_percentile_gain=all(r['percentiles'][0]>1 for r in rows))
    r=dict(rows=rows,gates=gates,draws=2000,seed=20260914,packet_histories=n,active_histories=k,
           scope='Paired empirical packet bootstrap; no certified coverage, anatomy or physiological uncertainty.')
    a.output.write_text(json.dumps(r,indent=2)+'\n');np.savez_compressed(a.output.with_suffix('.npz'),**arrays)
    print(json.dumps(r));return 0 if all(gates.values()) else 2


if __name__=='__main__':raise SystemExit(main())
