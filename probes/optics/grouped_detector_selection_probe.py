"""Replay frozen atlas candidate scores through the reusable group selector."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from bodytwin.geometry.optics.grouped_detector_selection_v1 import select_groups, selftest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--selection',type=Path,required=True)
    p.add_argument('--training',type=Path,required=True)
    p.add_argument('--confirmation',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    specification=json.loads(a.selection.read_text())
    rows=[];arrays={}
    with np.load(a.training,allow_pickle=False) as training, np.load(a.confirmation,allow_pickle=False) as confirmation:
        for leg in (0,1):
            j=training[f'training_{leg}']
            selected=select_groups(j[:,:,[0]],j[:,:,[2]],j[:,:,[1]],specification['training_triples'])
            difference=np.abs(selected.scores-training[f'scores_{leg}'])
            fresh=[];fresh_errors=[]
            for name in ('selected','reference'):
                matrix=confirmation[f'{name}_jacobian_{leg}']
                interest=matrix[:,0].reshape(3,3,1)
                shared=matrix[:,1].reshape(3,3,1)
                local=np.array([matrix[3*i:3*i+3,2+i] for i in range(3)])[:,:,None]
                result=select_groups(interest,shared,local,[(0,1,2)])
                expected=np.linalg.norm(confirmation[f'{name}_profile_{leg}'])
                fresh.append(float(result.scores[0]));fresh_errors.append(float(abs(result.scores[0]-expected)))
            rows.append(dict(selected_ids=list(selected.indices),selftest=selftest(),maximum_score_difference=float(difference.max()),
                             training_agreement=bool(np.all(difference<=np.maximum(1e-12,1e-10*np.abs(training[f'scores_{leg}'])))),
                             fresh_sensitivities=fresh,fresh_errors=fresh_errors,
                             fresh_agreement=all(error<=max(1e-12,1e-10*score) for error,score in zip(fresh_errors,fresh))))
            arrays[f'scores_{leg}']=selected.scores
            arrays[f'ranks_{leg}']=selected.ranks
            arrays[f'fresh_{leg}']=np.array(fresh)
    gates=dict(selftest=all(all(r['selftest'].values()) for r in rows),
               frozen_selection=all(r['selected_ids']==[4,7,9] for r in rows),
               numerical_agreement=all(r['training_agreement'] and r['fresh_agreement'] for r in rows),
               full_repeat=rows[0]==rows[1] and all(v.tobytes()==arrays[k[:-1]+'1'].tobytes() for k,v in arrays.items() if k.endswith('_0')))
    report=dict(rows=rows,gates=gates,input_hashes={name:hashlib.sha256(getattr(a,name).read_bytes()).hexdigest() for name in ('selection','training','confirmation')},
                hashes={k:hashlib.sha256(v.tobytes()).hexdigest() for k,v in arrays.items()},
                scope='Public archived Jacobian replay of the reusable group selector, not a fresh transport run or new physiological claim.')
    np.savez(a.output.with_suffix('.npz'),**arrays)
    a.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
    return 0 if all(gates.values()) else 2


if __name__=='__main__':raise SystemExit(main())
