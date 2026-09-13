"""Check neighbor traversal against the frozen exhaustive interval oracle."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(Path(__file__).resolve().parent))
from bodytwin.geometry.tetra_ray_walk_v1 import TetraRayWalk
from tetra_partition_mechanism import fixtures
from tetra_ray_mechanism import intervals,observe


def main():
    _,v,t,labels=fixtures()[1];_,reference=observe();origin=reference['origin'];directions=reference['directions']
    rows=[];arrays={}
    for name,cells in [('original',t),('inverted',t[:,::-1])]:
        walker=TetraRayWalk(v,cells,labels)
        for leg in (0,1):
            records=[];max_error=0.;sequence=True
            for direction in directions:
                got=walker.trace(origin,direction);expected,_=intervals(v,t,origin,direction)
                sequence &= got.shape==expected.shape and np.array_equal(got[:,0],expected[:,0])
                if got.shape==expected.shape:max_error=max(max_error,float(np.max(np.abs(got[:,1:]-expected[:,1:]))))
                else:max_error=float('inf')
                records.extend(got.tolist())
            a=np.asarray(records);arrays[f'{name}_{leg}']=a;rows.append(dict(case=name,leg=leg,sequence=bool(sequence),max_parameter_error=max_error,sha256=hashlib.sha256(a.tobytes()).hexdigest()))
    controls={}
    walker=TetraRayWalk(v,t,labels)
    cases={'edge_source':lambda:walker.trace([0,0,0],[1,0,0]),'edge_exit':lambda:walker.trace(origin,[-1,-1,0]),'zero_direction':lambda:walker.trace(origin,[0,0,0]),'nan_direction':lambda:walker.trace(origin,[np.nan,0,1]),'duplicate_cells':lambda:TetraRayWalk(v,np.repeat(t[:1],2,axis=0),[1,1])}
    for name,fn in cases.items():
        try:fn();controls[name]=False
        except ValueError:controls[name]=True
    gates=dict(oracle=all(r['sequence'] and r['max_parameter_error']<=1e-12 for r in rows),repeat=all(arrays[f'{name}_0'].tobytes()==arrays[f'{name}_1'].tobytes() for name in ('original','inverted')),inversion=np.array_equal(arrays['original_0'][:,0],arrays['inverted_0'][:,0]) and float(np.max(np.abs(arrays['original_0'][:,1:]-arrays['inverted_0'][:,1:])))<=1e-12,singular_refused=controls['edge_source'] and controls['edge_exit'],invalid_refused=all(controls[k] for k in ('zero_direction','nan_direction','duplicate_cells')))
    report=dict(rows=rows,gates=gates,controls=controls,scope='Synthetic conforming straight-ray cell traversal only; no invalid-atlas or photon transport promotion.')
    (ROOT/'reports/tetra_ray_walk.json').write_text(json.dumps(report,indent=2)+'\n');np.savez_compressed(ROOT/'reports/tetra_ray_walk_arrays.npz',**arrays);print(json.dumps(report));return 0 if all(gates.values()) else 2

if __name__=='__main__':raise SystemExit(main())
