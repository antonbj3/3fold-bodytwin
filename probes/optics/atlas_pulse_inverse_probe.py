"""Compare frozen-path-ratio and path-dependent pulse inverses on held-out paths."""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from bodytwin.geometry.optics.pulse_path_inverse_v1 import PulsePathInverse
from atlas_oxygen_bounds import EXTINCTION,FRACTION,BACKGROUND


def selected(filename):
    with np.load(filename,allow_pickle=False) as c:
        radius=np.linalg.norm(c['exits'][:,:3]-c['source'],axis=1)
        return c['paths'][(c['terminal'][:,1]>0)&(radius>=5)&(radius<40)]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--training',required=True);p.add_argument('--independent',required=True);p.add_argument('--output',required=True,type=Path);a=p.parse_args()
    train=selected(a.training);independent=selected(a.independent)
    model=PulsePathInverse(train,BACKGROUND,FRACTION,EXTINCTION)
    target=PulsePathInverse(independent,BACKGROUND,FRACTION,EXTINCTION)
    e=EXTINCTION[:,0]*.8+EXTINCTION[:,1]*.2
    fixed=e[0]/e[1]/model.ratio(.8)
    legs=[]
    for _ in (0,1):
        rows=[]
        for sat in (.55,.65,.75,.85,.95):
            observed=target.ratio(sat);prediction=model.solve(observed)
            corrected=observed*fixed
            old=(EXTINCTION[0,1]-corrected*EXTINCTION[1,1])/(corrected*(EXTINCTION[1,0]-EXTINCTION[1,1])-(EXTINCTION[0,0]-EXTINCTION[0,1]))
            rows.append(dict(saturation=sat,observed_ratio=observed,prediction=prediction,error=abs(prediction-sat),frozen_prediction=float(old),frozen_error=float(abs(old-sat)),self_error=abs(model.solve(model.ratio(sat))-sat)))
        legs.append(rows)
    rejected=0
    for bad in (-1.,np.nan,100.):
        try:model.solve(bad)
        except ValueError:rejected+=1
    gates=dict(full_repeat=legs[0]==legs[1],self_consistency=all(r['self_error']<=1e-10 for r in legs[0]),independent_error=all(r['error']<=.01 for r in legs[0]),refusals=rejected==3)
    report=dict(legs=legs,gates=gates,refused_controls=rejected,training_count=len(train),independent_count=len(independent),path_hashes=[hashlib.sha256(x.tobytes()).hexdigest() for x in (train,independent)],scope='Held-out path-library optical-model inversion; no clinical or physiological acceptance.')
    a.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));return 0 if all(gates.values()) else 2


if __name__=='__main__':raise SystemExit(main())
