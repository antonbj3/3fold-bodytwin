"""Homogeneous log-signal control for two-band concentration ambiguity."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from atlas_detector_third_wavelength import COEFFICIENTS


def design(bands):
    out=np.zeros((3*len(bands),5))
    for detector,length in enumerate((10.,20.,30.)):
        block=slice(detector*len(bands),(detector+1)*len(bands))
        out[block,:2]=-length*np.log(10.)/64500.*COEFFICIENTS[bands]
        out[block,detector+2]=1
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    legs,arrays=[],{}
    for leg in (0,1):
        two=design([0,2]);three=design([0,1,2])
        delta=np.zeros(5);delta[0]=.1
        delta[1]=-(COEFFICIENTS[0,0]-COEFFICIENTS[2,0])/(COEFFICIENTS[0,1]-COEFFICIENTS[2,1])*delta[0]
        delta[2:]=-two[::2,:2]@delta[:2]
        rows=[];signals=[];solutions=[]
        for sat in (.7,.8,.9):
            for hb in (12.,15.,18.):
                truth=np.array([sat*hb,(1-sat)*hb,.1,-.2,.3])
                # Generate forward signals separately from the inverse design multiplication.
                log_signal=[]
                for detector,length in enumerate((10.,20.,30.)):
                    for eo,ed in COEFFICIENTS:
                        log_signal.append(truth[detector+2]-length*(.003+np.log(10.)/64500.*hb*(sat*eo+(1-sat)*ed)))
                signal=np.array(log_signal)+np.repeat(np.array([10.,20.,30.])*.003,3)
                recovered=np.linalg.lstsq(three,signal,rcond=None)[0]
                recovered_hb=float(recovered[:2].sum());recovered_sat=float(recovered[0]/recovered_hb)
                rows.append(dict(saturation=sat,hemoglobin=hb,saturation_error=abs(recovered_sat-sat),hemoglobin_error=abs(recovered_hb-hb)))
                signals.append(signal);solutions.append(recovered)
        legs.append(dict(rows=rows,two_rank=int(np.linalg.matrix_rank(two)),three_rank=int(np.linalg.matrix_rank(three)),
                         two_null_difference=float(np.max(np.abs(two@delta))),three_null_difference=float(np.max(np.abs(three@delta)))))
        for name,value in [('two',two),('three',three),('delta',delta),('signals',signals),('solutions',solutions)]:
            arrays[f'{name}_{leg}']=np.asarray(value)
    gates=dict(full_repeat=legs[0]==legs[1] and all(v.tobytes()==arrays[k[:-1]+'1'].tobytes() for k,v in arrays.items() if k.endswith('_0')),
               ranks=all(r['two_rank']==4 and r['three_rank']==5 for r in legs),
               recovery=all(row['saturation_error']<=1e-10 and row['hemoglobin_error']<=1e-8 for r in legs for row in r['rows']),
               null_control=all(r['two_null_difference']<=1e-12 and r['three_null_difference']>1e-6 for r in legs))
    result=dict(legs=legs,gates=gates,hashes={k:hashlib.sha256(v.tobytes()).hexdigest() for k,v in arrays.items()},
                scope='Homogeneous analytic signals with wavelength-shared unknown gains; no clinical or noisy-data recovery claim.')
    np.savez(a.output.with_suffix('.npz'),**arrays)
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
    return 0 if all(gates.values()) else 2


if __name__=='__main__':
    raise SystemExit(main())
