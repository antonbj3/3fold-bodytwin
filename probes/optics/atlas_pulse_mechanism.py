"""Measure saturation bias from fixing a pulse pathlength ratio."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from atlas_oxygen_bounds import EXTINCTION,FRACTION,BACKGROUND


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('capture');p.add_argument('--output',required=True,type=Path);a=p.parse_args()
    with np.load(a.capture,allow_pickle=False) as c:
        radius=np.linalg.norm(c['exits'][:,:3]-c['source'],axis=1)
        mask=(c['terminal'][:,1]>0)&(radius>=5)&(radius<40)
        paths=c['paths'][mask];unfinished=c['paths'][c['terminal'][:,2]>0];n=len(c['paths'])
    scale=np.log(10.)*15/64500.;saturations=[.5,.6,.7,.8,.9,1.]
    legs=[];arrays={}
    for leg in (0,1):
        rows=[];all_weights=[]
        for sat in saturations:
            extinction=EXTINCTION[:,0]*sat+EXTINCTION[:,1]*(1-sat)
            blood=scale*extinction
            base_mu=BACKGROUND[None,:]+blood[:,None]*FRACTION
            pulse_mu=BACKGROUND[None,:]+1.01*blood[:,None]*FRACTION
            w0=np.exp(-(paths@base_mu.T));w1=np.exp(-(paths@pulse_mu.T))
            i0=w0.sum(axis=0)/n;i1=w1.sum(axis=0)/n
            od=np.log(i0/i1);effective=od/(.01*blood)
            tail=np.exp(-(unfinished@base_mu.T)).sum(axis=0)/n
            rows.append(dict(saturation=sat,pulse_od=od.tolist(),od_ratio=float(od[0]/od[1]),
                             effective_weighted_lengths=effective.tolist(),path_ratio=float(effective[1]/effective[0]),
                             unfinished_intensity_upper=tail.tolist()))
            all_weights.append(np.array([w0,w1]))
        fixed=rows[3]['path_ratio']
        for row in rows:
            ratio=row['od_ratio']*fixed
            value=(EXTINCTION[0,1]-ratio*EXTINCTION[1,1])/(ratio*(EXTINCTION[1,0]-EXTINCTION[1,1])-(EXTINCTION[0,0]-EXTINCTION[0,1]))
            row['constant_ratio_saturation']=float(value);row['absolute_error']=float(abs(value-row['saturation']))
        arrays[f'weights_{leg}']=np.array(all_weights);legs.append(rows)
    gates=dict(full_repeat=legs[0]==legs[1] and arrays['weights_0'].tobytes()==arrays['weights_1'].tobytes(),
               positive_pulse=all(np.isfinite(r['pulse_od']).all() and min(r['pulse_od'])>0 for leg in legs for r in leg),
               anchor=all(leg[3]['absolute_error']<=1e-12 for leg in legs),
               constant_ratio_bias=all(r['absolute_error']<=.01 for leg in legs for r in leg))
    report=dict(legs=legs,gates=gates,selected_histories=len(paths),relative_blood_fraction_pulse=.01,
                weight_hashes=[hashlib.sha256(arrays[f'weights_{i}'].tobytes()).hexdigest() for i in (0,1)],
                source='https://doi.org/10.1117/1.JBO.28.11.115002',
                scope='Illustrative fixed-geometry pulse/path approximation; no matching clinical calibration or source-experiment reproduction.')
    a.output.write_text(json.dumps(report,indent=2)+'\n');np.savez_compressed(a.output.with_suffix('.npz'),**arrays)
    print(json.dumps(report));return 0 if all(gates.values()) else 2


if __name__=='__main__':raise SystemExit(main())
