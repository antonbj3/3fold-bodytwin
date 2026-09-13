"""Compare reusable nuisance profiling with archived atlas QR observations."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from bodytwin.geometry.optics.profiled_information_v1 import profile_information, selftest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('arrays',type=Path)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    rows, output=[],{}
    with np.load(a.arrays,allow_pickle=False) as source:
        for leg in (0,1):
            results=[];differences=[];agreement=[]
            for index,matrix in enumerate(source[f'expanded_{leg}']):
                result=profile_information(matrix[:,:1],matrix[:,1:])
                expected=np.linalg.norm(source[f'profiles_{leg}'][index,1])
                difference=abs(np.linalg.norm(result.projected)-expected)
                differences.append(float(difference));agreement.append(difference<=max(1e-12,1e-10*expected))
                results.append(np.concatenate((result.projected.ravel(),result.gram.ravel(),result.singular_values,
                                               [result.nuisance_rank,result.interest_rank,result.nuisance_cutoff,result.interest_cutoff])))
            output[f'results_{leg}']=np.array(results)
            rows.append(dict(selftest=selftest(),maximum_difference=max(differences),agreement=bool(all(agreement)),cases=len(results)))
    gates=dict(selftest=all(all(r['selftest'].values()) for r in rows),
               atlas_agreement=all(r['agreement'] for r in rows),
               full_repeat=rows[0]==rows[1] and output['results_0'].tobytes()==output['results_1'].tobytes())
    report=dict(rows=rows,gates=gates,input_sha256=hashlib.sha256(a.arrays.read_bytes()).hexdigest(),
                hashes={k:hashlib.sha256(v.tobytes()).hexdigest() for k,v in output.items()},
                scope='Numerical nuisance projection and analytic controls; no clinical or interval identifiability certificate.')
    np.savez(a.output.with_suffix('.npz'),**output)
    a.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
    return 0 if all(gates.values()) else 2


if __name__=='__main__':
    raise SystemExit(main())
