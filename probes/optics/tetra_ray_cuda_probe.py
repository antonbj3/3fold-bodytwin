"""Full CUDA traversal comparison plus singular/capacity controls."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(Path(__file__).resolve().parent))
from bodytwin.geometry.tetra_ray_cuda_v1 import trace
from bodytwin.geometry.tetra_ray_walk_v1 import TetraRayWalk
from tetra_ray_grid_probe import fixture
from tetra_partition_mechanism import fixtures


def main():
    v,t,labels=fixture();prepared=TetraRayWalk(v,t,labels);origin=np.array([.17,.23,.31]);directions=np.random.default_rng(20260913).normal(size=(512,3));directions/=np.linalg.norm(directions,axis=1)[:,None]
    expected=[prepared.trace(origin,d) for d in directions];rows=[];arrays={}
    for leg in (0,1):
        out=trace(prepared,origin,directions);sequence=True;error=0.
        for ray,ref in enumerate(expected):
            n=int(out['counts'][ray]);sequence &= n==len(ref) and np.array_equal(out['cells'][ray,:n],ref[:,0])
            if n==len(ref):error=max(error,float(np.max(np.abs(out['intervals'][ray,:n]-ref[:,1:]))))
            else:error=float('inf')
        rows.append(dict(leg=leg,sequence=bool(sequence),max_parameter_error=error,status_histogram={str(int(x)):int(n) for x,n in zip(*np.unique(out['status'],return_counts=True))},hashes={k:hashlib.sha256(a.tobytes()).hexdigest() for k,a in out.items()}))
        for k,a in out.items():arrays[f'{k}_{leg}']=a
    controls=[]
    for leg in (0,1):
        capped=trace(prepared,origin,directions,capacity=2);expected_cap=np.array([len(x)>2 for x in expected]);capacity_exact=bool(np.array_equal(capped['status']==4,expected_cap) and np.all(capped['status'][~expected_cap]==0))
        _,sv,st,sl=fixtures()[1];small=TetraRayWalk(sv,st,sl);singular=trace(small,sv[st[0]].mean(axis=0),np.array([[-1.,-1.,0.]]))
        controls.append(dict(capacity_exact=capacity_exact,capacity_failures=int(expected_cap.sum()),singular_status=int(singular['status'][0]),hashes={k:hashlib.sha256(a.tobytes()).hexdigest() for k,a in capped.items()}))
        for k,a in capped.items():arrays[f'capped_{k}_{leg}']=a
    gates=dict(reference=all(r['sequence'] and r['max_parameter_error']<=1e-12 for r in rows),full_repeat=rows[0]['hashes']==rows[1]['hashes'] and controls[0]==controls[1],normal_status=all(r['status_histogram']=={'0':512} for r in rows),explicit_controls=all(r['capacity_exact'] and r['singular_status']==2 for r in controls))
    report=dict(rows=rows,controls=controls,gates=gates,scope='Synthetic CUDA straight-ray cell walk only; no photon/refraction/anatomy or speed claim.')
    (ROOT/'reports/tetra_ray_cuda.json').write_text(json.dumps(report,indent=2)+'\n');np.savez_compressed(ROOT/'reports/tetra_ray_cuda_arrays.npz',**arrays);print(json.dumps(report));return 0 if all(gates.values()) else 2

if __name__=='__main__':raise SystemExit(main())
