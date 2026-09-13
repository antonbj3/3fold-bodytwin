"""Propagate capped-history moment envelopes through frozen nuisance profiling."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import numpy as np
from bodytwin.geometry.optics.profiled_information_bound_v1 import profile_bound,selftest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--tails',type=Path,required=True)
    p.add_argument('--jacobians',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();tails=json.loads(a.tails.read_text());rows=[];arrays={}
    if not all(v is True for v in tails['gates'].values()):raise ValueError('Tail mechanism gates failed')
    with np.load(a.jacobians,allow_pickle=False) as archive:
        for leg in (0,1):
            observed=np.array(tails['legs'][leg]['observed']);delta=np.array(tails['legs'][leg]['remainder'])
            pair=[]
            for start,name in ((0,'selected'),(3,'reference')):
                nominal=np.zeros((9,7));radius=np.zeros_like(nominal)
                for detector in range(3):
                    signal,ds,dh=observed[start+detector];di,dds,ddh=delta
                    if np.any(signal<=0) or np.any(delta<0):raise ValueError('Positive signal and nonnegative envelopes required')
                    root=np.sqrt(signal);upper=np.sqrt(signal+di)
                    dr=di/(root+upper)
                    inverse_difference=di/(root*upper*(root+upper))
                    sl=slice(detector*3,detector*3+3)
                    nominal[sl,0]=ds/root;nominal[sl,1]=dh/root;nominal[sl,detector+2]=root
                    radius[sl,0]=dds/root+np.abs(ds)*inverse_difference
                    radius[sl,1]=ddh/root+np.abs(dh)*inverse_difference;radius[sl,detector+2]=dr
                    nominal[detector*3+1,5]=root[1];nominal[detector*3+2,6]=root[2]
                    radius[detector*3+1,5]=dr[1];radius[detector*3+2,6]=dr[2]
                frozen=archive[f'{name}_extended_{leg}'];alignment=np.abs(nominal-frozen)
                radius+=alignment
                bound=profile_bound(frozen[:,:1],frozen[:,1:],radius[:,:1],radius[:,1:])
                pair.append(dict(name=name,bound=asdict(bound),maximum_alignment=float(alignment.max())))
                arrays[f'{name}_radius_{leg}']=radius;arrays[f'{name}_nominal_{leg}']=nominal
            rows.append(dict(pair=pair,gain_lower=pair[0]['bound']['lower']/pair[1]['bound']['upper']))
    gates=dict(selftest=all(selftest().values()),full_repeat=rows[0]==rows[1] and all(v.tobytes()==arrays[k[:-1]+'1'].tobytes() for k,v in arrays.items() if k.endswith('_0')),
               alignment=all(r['maximum_alignment']<=1e-10 for leg in rows for r in leg['pair']),gain_lower=all(r['gain_lower']>1 for r in rows))
    report=dict(rows=rows,gates=gates,input_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (a.tails,a.jacobians)},
                hashes={k:hashlib.sha256(v.tobytes()).hexdigest() for k,v in arrays.items()},
                scope='Evaluated capped-history perturbation envelope at fixed optical point; no outward-rounded interval, physiological, geometric or calibration certificate.')
    np.savez(a.output.with_suffix('.npz'),**arrays);a.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(rows=rows,gates=gates)))
    return 0 if all(gates.values()) else 2


if __name__=='__main__':raise SystemExit(main())
