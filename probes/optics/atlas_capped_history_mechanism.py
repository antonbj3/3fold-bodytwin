"""Measure conservative signal and derivative envelopes for capped histories."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from atlas_detector_support import candidates
from atlas_detector_third_wavelength import COEFFICIENTS
from atlas_oxygen_bounds import FRACTION,BACKGROUND


def envelope(depth):
    return np.where(depth<=1,1/np.e,depth*np.exp(-depth))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original',required=True)
    p.add_argument('captures',nargs=2)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    with np.load(a.original,allow_pickle=False) as source:
        centers,_=candidates(source['exits'],source['terminal'],source['source'])
    scale=np.log(10.)*15/64500.
    blood=(scale*(COEFFICIENTS[:,0]*.8+COEFFICIENTS[:,1]*.2))[:,None]*FRACTION
    derivative=(scale*(COEFFICIENTS[:,0]-COEFFICIENTS[:,1]))[:,None]*FRACTION
    mu=blood+BACKGROUND
    k_s=np.max(np.divide(np.abs(derivative),mu,out=np.zeros_like(mu),where=mu>0),axis=1)
    k_h=np.max(np.divide(blood,mu,out=np.zeros_like(mu),where=mu>0),axis=1)
    legs=[]
    for filename in a.captures:
        with np.load(filename,allow_pickle=False) as source:
            paths=source['paths'];terminal=source['terminal'];exits=source['exits']
        capped=terminal[:,2]>0;depth=paths@mu.T;w=np.exp(-depth)*(2/3)
        moments=(w,-(paths@derivative.T)*w,-(paths@blood.T)*w)
        masks=(np.linalg.norm(exits[:,None,:3]-centers[None,[4,7,9,0,4,8]],axis=2)<=5)&(terminal[:,1]>0)[:,None]
        observed=np.array([[v[m].sum(axis=0) for v in moments] for m in masks.T])
        remainder=np.array([np.exp(-depth[capped]).sum(axis=0),k_s*envelope(depth[capped]).sum(axis=0),k_h*envelope(depth[capped]).sum(axis=0)])*(2/3)
        relative=remainder[None]/np.abs(observed)
        legs.append(dict(histories=len(paths),caps=int(capped.sum()),path_hash=hashlib.sha256(paths.tobytes()).hexdigest(),
                         minimum_capped_depth=depth[capped].min(axis=0).tolist(),remainder=remainder.tolist(),
                         observed=observed.tolist(),relative=relative.tolist(),maximum_relative=float(relative.max())))
    floor=np.linspace(0,40,10001);test=np.max([(floor+d)*np.exp(-(floor+d))-envelope(floor) for d in (0,.01,1,10)])
    gates=dict(full_repeat=legs[0]==legs[1],input_binding=all(r['histories']==100000 and r['caps']==1901 and r['path_hash']=='f808f4ddba539f8f2071333593f33cb54ff4d5d7da8db8286db85aa43910afd6' for r in legs),
               finite_observations=all(np.isfinite(r['observed']).all() and np.all(np.abs(r['observed'])>0) for r in legs),
               envelope_control=bool(test<=1e-14),relative_remainder=all(r['maximum_relative']<=1e-6 for r in legs))
    report=dict(legs=legs,gates={k:bool(v) for k,v in gates.items()},envelope_control_maximum=float(test),
                scope='Fixed optical/scattering model, worst per-detector capped-history envelope; floating-point evaluation, not interval or clinical certification.')
    a.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(gates=report['gates'],maximum_relative=legs[0]['maximum_relative'],remainder=legs[0]['remainder'])))
    return 0 if all(gates.values()) else 2


if __name__=='__main__':raise SystemExit(main())
