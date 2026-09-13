"""Measure declared synthetic input uncertainty before a composed NIRS chain."""
import itertools
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from bodytwin.geometry.optics.nirs_forward_v1 import signal,sao2_hill,BLOOD_FRACTION
NOMINAL=np.array([60.,15.,1.]);BANDS=np.array([2.,.75,.1])


def evaluate(paths,values):
    po2,hb,fraction=values
    return signal(paths,float(sao2_hill(po2)),hb,BLOOD_FRACTION*fraction)[0]['ratio_760_850']


def measure(paths):
    nominal=evaluate(paths,NOMINAL);isolated=[];derivatives=[]
    for axis in range(3):
        values=[]
        for side in (-1,0,1):
            x=NOMINAL.copy();x[axis]+=side*BANDS[axis];values.append(evaluate(paths,x))
        isolated.append({'axis':axis,'ratios':values,'span':max(values)-min(values)})
        delta=NOMINAL[axis]*1e-5;x=NOMINAL.copy();x[axis]+=delta;plus=evaluate(paths,x)
        x[axis]-=2*delta;minus=evaluate(paths,x);derivatives.append((plus-minus)/(2*delta))
    rows=[]
    for factors in itertools.product((-1,0,1),repeat=3):
        values=NOMINAL+BANDS*np.array(factors);rows.append({'inputs':values.tolist(),'ratio':evaluate(paths,values)})
    joint=max(r['ratio'] for r in rows)-min(r['ratio'] for r in rows);linear=sum(r['span'] for r in isolated)
    degeneracy=abs(derivatives[1]*NOMINAL[1]-derivatives[2]*NOMINAL[2])/abs(derivatives[2]*NOMINAL[2])
    return {'nominal_ratio':nominal,'isolated':isolated,'grid':rows,'joint_span':joint,'sum_isolated_spans':linear,
        'composition_pass':joint<=linear,'composition_excess':joint-linear,'derivatives':derivatives,
        'hb_fraction_sensitivity_relative_difference':float(degeneracy)}


def main():
    paths=np.load(ROOT/'reports/photon_curved_detector_paths_l4.npz')['paths_0'];legs=[measure(paths),measure(paths)]
    gates={'tables_repeat_exact':legs[0]==legs[1],
        'finite_positive':all(np.isfinite(r['ratio']) and r['ratio']>0 for leg in legs for r in leg['grid']),
        'hb_fraction_degeneracy':all(r['hb_fraction_sensitivity_relative_difference']<=1e-7 for r in legs)}
    report={'legs':legs,'observation_gates':gates,'nominal':NOMINAL.tolist(),'half_bands':BANDS.tolist(),
        'scope':'declared synthetic input sensitivity; joint-span composition is measured separately, not a clinical confidence bound'}
    (ROOT/'reports/nirs_curved_uncertainty.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:legs[0][k] for k in ('nominal_ratio','joint_span','sum_isolated_spans','composition_pass','composition_excess','derivatives','hb_fraction_sensitivity_relative_difference')}))
    return 0 if all(gates.values()) else 1


if __name__=='__main__':raise SystemExit(main())
