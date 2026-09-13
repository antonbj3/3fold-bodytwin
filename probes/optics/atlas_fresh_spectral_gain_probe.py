"""Measure unknown common wavelength response on the frozen fresh placement."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from bodytwin.geometry.optics.profiled_information_v1 import profile_information


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('arrays',type=Path)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();rows=[];arrays={}
    with np.load(a.arrays,allow_pickle=False) as source:
        for leg in (0,1):
            row=[]
            for identity in ('selected','reference'):
                original=source[f'{identity}_jacobian_{leg}']
                extended=np.zeros((9,7));extended[:,:5]=original
                root_signal=original[:,2:].sum(axis=1)
                extended[1::3,5]=root_signal[1::3]
                extended[2::3,6]=root_signal[2::3]
                baseline=profile_information(original[:,:1],original[:,1:])
                result=profile_information(extended[:,:1],extended[:,1:])
                before=float(np.linalg.norm(baseline.projected));after=float(np.linalg.norm(result.projected))
                row.append(dict(identity=identity,before=before,after=after,retained_fraction=after/before,
                                rank=int(np.linalg.matrix_rank(extended)),nuisance_rank=result.nuisance_rank))
                arrays[f'{identity}_extended_{leg}']=extended
                arrays[f'{identity}_projected_{leg}']=result.projected
                arrays[f'{identity}_spectrum_{leg}']=np.linalg.svd(extended,compute_uv=False)
            rows.append(dict(triples=row,placement_gain=row[0]['after']/row[1]['after']))
    gates=dict(full_repeat=rows[0]==rows[1] and all(v.tobytes()==arrays[k[:-1]+'1'].tobytes() for k,v in arrays.items() if k.endswith('_0')),
               finite_rank=all(np.isfinite(v).all() for v in arrays.values()) and all(t['rank']==7 for r in rows for t in r['triples']),
               monotonic=all(t['after']<=t['before']*(1+1e-12) for r in rows for t in r['triples']),
               placement_gain=all(r['placement_gain']>1 for r in rows))
    report=dict(rows=rows,gates={k:bool(v) for k,v in gates.items()},input_sha256=hashlib.sha256(a.arrays.read_bytes()).hexdigest(),
                hashes={k:hashlib.sha256(v.tobytes()).hexdigest() for k,v in arrays.items()},
                scope='Common relative spectral response nuisance on fixed assumed-optics matrices; not new calibration or a clinical certificate.')
    np.savez(a.output.with_suffix('.npz'),**arrays)
    a.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(rows=rows,gates=report['gates'])))
    return 0 if all(gates.values()) else 2


if __name__=='__main__':
    raise SystemExit(main())
