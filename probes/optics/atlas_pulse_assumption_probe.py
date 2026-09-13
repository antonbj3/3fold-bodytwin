"""Expose optical-assumption sensitivity of the frozen pulse path inverse."""
import argparse,json
from pathlib import Path
from atlas_pulse_inverse_probe import selected,PulsePathInverse,BACKGROUND,FRACTION,EXTINCTION


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--training',required=True);p.add_argument('--independent',required=True);p.add_argument('--output',required=True,type=Path);a=p.parse_args()
    train=selected(a.training);test=selected(a.independent)
    inverse=PulsePathInverse(train,BACKGROUND,FRACTION,EXTINCTION)
    conditions=[('matched',15.,.01,1.),('hb_low',12.,.01,1.),('hb_high',18.,.01,1.),('pulse_low',15.,.005,1.),('pulse_high',15.,.02,1.),('background_low',15.,.01,.8),('background_high',15.,.01,1.2)]
    legs=[]
    for _ in (0,1):
        rows=[]
        for name,hb,pulse,background in conditions:
            target=PulsePathInverse(test,BACKGROUND*background,FRACTION,EXTINCTION,hb_g_dl=hb,pulse_fraction=pulse)
            for sat in (.55,.65,.75,.85,.95):
                observation=target.ratio(sat)
                try:prediction=inverse.solve(observation);error=abs(prediction-sat);refused=False
                except ValueError:prediction=None;error=None;refused=True
                rows.append(dict(condition=name,saturation=sat,hb=hb,pulse=pulse,background_scale=background,observed_ratio=observation,prediction=prediction,error=error,refused=refused))
        legs.append(rows)
    gates=dict(full_repeat=legs[0]==legs[1],matched=all(not r['refused'] and r['error']<=.01 for r in legs[0] if r['condition']=='matched'),all_assumptions=all(not r['refused'] and r['error']<=.01 for r in legs[0]))
    summary={name:dict(refusals=sum(r['refused'] for r in legs[0] if r['condition']==name),maximum_error=max((r['error'] for r in legs[0] if r['condition']==name and r['error'] is not None),default=None)) for name,*_ in conditions}
    result=dict(legs=legs,gates=gates,summary=summary,scope='One-factor optical-assumption stress; numerical ranges only, no clinical calibration.')
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(gates=gates,summary=summary)));return 0 if all(gates.values()) else 2


if __name__=='__main__':raise SystemExit(main())
