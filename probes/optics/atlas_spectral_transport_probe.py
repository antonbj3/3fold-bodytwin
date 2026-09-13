"""Compare direct atlas spectral transport with fixed proposal path weights."""
import hashlib
import numpy as np
from atlas_oxygen_bounds import EXTINCTION,FRACTION,BACKGROUND
PACKET=2**30


def absorption(saturation):
    blood=np.log(10.)*15/64500.*(EXTINCTION[:,0]*saturation+EXTINCTION[:,1]*(1-saturation))
    return (BACKGROUND[None,:]+blood[:,None]*FRACTION).astype(np.float32)


def assess(out,proposal,source,mua):
    radius=np.linalg.norm(out['exits'][:,:3]-source,axis=1)
    collected=(out['terminal'][:,1]>0)&(radius>=5)&(radius<40)
    pr=np.linalg.norm(proposal['exits'][:,:3]-source,axis=1)
    chosen=(proposal['terminal'][:,1]>0)&(pr>=5)&(pr<40)
    unfinished=proposal['terminal'][:,2]>0
    weights=np.exp(-(proposal['paths']@np.asarray(mua,dtype=np.float64)))
    n=len(out['terminal'])
    actual=float(out['terminal'][collected,1].sum(dtype=np.int64))/(n*PACKET)
    predicted=float(weights[chosen].sum())/n;tail=float(weights[unfinished].sum())/n
    t=out['terminal'];energy=t[:,:3].sum(axis=0,dtype=np.int64)
    return dict(hashes={k:hashlib.sha256(v.tobytes()).hexdigest() for k,v in out.items()},
                photons=n,absorbed=int(energy[0]),escaped=int(energy[1]),residual=int(energy[2]),
                leaks=int(t[:,3].sum()),caps=int(out['counters'][2]),detected=int(collected.sum()),
                intensity=actual,predicted_intensity=predicted,unfinished_upper=tail,
                intensity_difference=abs(actual-predicted),bound=1/PACKET+tail,mua=np.asarray(mua).tolist())


def evaluate(rows):
    groups=[[r for r in rows if r['saturation']==sat and r['wavelength']==w] for sat in (.7,.8,.9,1.) for w in (760,850)]
    if any(len(g)!=2 for g in groups):raise ValueError('Two complete captures for every spectral arm required')
    ratios=[]
    for i,sat in enumerate((.7,.8,.9,1.)):
        ratios.append(dict(saturation=sat,ratio=groups[2*i][0]['intensity']/groups[2*i+1][0]['intensity']))
    gates=dict(full_repeat=all(g[0]==g[1] for g in groups),energy=all(r['absorbed']+r['escaped']+r['residual']==r['photons']*PACKET for r in rows),
               zero_failures=all(r['leaks']==r['caps']==r['residual']==0 for r in rows),
               direct_path_agreement=all(r['intensity_difference']<=r['bound'] for r in rows),
               monotonic=all(a['ratio']<b['ratio'] for a,b in zip(ratios,ratios[1:])),matching_calibration=False)
    return dict(rows=rows,ratios=ratios,gates=gates,matching_calibration_anchors=0,
                scope='Direct transport under assumed spectral properties; represented-coefficient path comparison, no clinical or anatomical acceptance.')


if __name__=='__main__':
    import argparse,json
    from pathlib import Path
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('capture');p.add_argument('proposal');p.add_argument('--mua',nargs=6,type=float,required=True);p.add_argument('--output',required=True,type=Path)
    a=p.parse_args()
    with np.load(a.capture,allow_pickle=False) as c,np.load(a.proposal,allow_pickle=False) as q:
        row=assess({k:c[k] for k in c.files if k!='source'},{k:q[k] for k in q.files if k!='source'},c['source'],np.array(a.mua,dtype=np.float32))
    a.output.write_text(json.dumps(row,indent=2)+'\n')
    raise SystemExit(0 if row['intensity_difference']<=row['bound'] and row['leaks']==row['caps']==row['residual']==0 else 2)
