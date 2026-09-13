"""Measure float32 boundary offsets and integer per-segment attenuation."""
import json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]


def measure():
    offsets=[]
    for p in (0.,60.,100.):
        for component in (-1.,-.1,-.01,-.001):
            before=np.float32(p);after=np.float32(before+np.float32(component)*np.float32(1e-5))
            offsets.append({'position':p,'direction_component':component,'nominal_offset':float(np.float32(component)*np.float32(1e-5)),
                            'represented_offset':float(after)-float(before),'unchanged':bool(before==after)})
    absorption=[]
    for weight in (1,32,64,128,512,1024,1073741824):
        mua=.01;length=.1
        raw=weight*(-np.expm1(-mua*length));integer=int(np.floor(raw+.5))
        minimum_length=-np.log1p(-.5/weight)/mua
        absorption.append({'weight':weight,'mua':mua,'length':length,'continuous_loss':float(raw),
                           'rounded_loss':integer,'minimum_length_for_nonzero_loss':float(minimum_length)})
    return {'offsets':offsets,'absorption':absorption}


def main():
    legs=[measure(),measure()]
    gates={'two_tables_exact':legs[0]==legs[1],
           'finite':all(np.isfinite(r['minimum_length_for_nonzero_loss']) for r in legs[0]['absorption'])}
    report={'legs':legs,'gates':gates,'scope':'scalar arithmetic diagnostic; no transport or failure attribution without canonical observer evidence'}
    (ROOT/'reports/photon_arithmetic_mechanism.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True);return 0 if all(gates.values()) else 1


if __name__=='__main__':raise SystemExit(main())
