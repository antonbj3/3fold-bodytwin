"""Larger conforming network control for the frozen neighbor walker."""
import hashlib
import itertools
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(Path(__file__).resolve().parent))
from bodytwin.geometry.tetra_ray_walk_v1 import TetraRayWalk
from tetra_ray_mechanism import intervals


def fixture(n=6):
    vertices=np.array(list(itertools.product(range(n+1),repeat=3)),float)
    def index(p):return (p[0]*(n+1)+p[1])*(n+1)+p[2]
    cells=[];labels=[]
    for corner in itertools.product(range(n),repeat=3):
        for permutation in itertools.permutations(range(3)):
            point=np.array(corner);tet=[index(point)]
            for axis in permutation:point[axis]+=1;tet.append(index(point))
            cells.append(tet);labels.append(1+sum(corner)%2)
    return vertices,np.array(cells,np.int32),np.array(labels,np.int32)


def observe(v,t,labels):
    walker=TetraRayWalk(v,t,labels);origin=np.array([.17,.23,.31]);directions=np.random.default_rng(20260913).normal(size=(512,3));directions/=np.linalg.norm(directions,axis=1)[:,None]
    values=[];offsets=[0];refusals=[];sequence_errors=0;parameter_error=0.;exit_error=0.;reference_intervals=0
    for ray,direction in enumerate(directions):
        expected,_=intervals(v,t,origin,direction);reference_intervals+=len(expected)
        try:got=walker.trace(origin,direction)
        except ValueError as error:
            refusals.append(dict(ray=ray,reason=str(error)));offsets.append(len(values));continue
        if got.shape!=expected.shape or not np.array_equal(got[:,0],expected[:,0]):sequence_errors+=1
        if got.shape==expected.shape:parameter_error=max(parameter_error,float(np.max(np.abs(got[:,1:]-expected[:,1:]))))
        else:parameter_error=float('inf')
        candidates=[((6-origin[i])/direction[i] if direction[i]>0 else -origin[i]/direction[i]) for i in range(3) if direction[i]!=0]
        exit_error=max(exit_error,abs(float(got[-1,2])-min(candidates)))
        values.extend(got.tolist());offsets.append(len(values))
    arrays=dict(intervals=np.asarray(values),offsets=np.asarray(offsets,np.int64),directions=directions)
    return dict(cells=len(t),rays=len(directions),reference_intervals=reference_intervals,refusals=refusals,sequence_errors=sequence_errors,maximum_parameter_error=parameter_error,maximum_exit_error=exit_error,hashes={k:hashlib.sha256(a.tobytes()).hexdigest() for k,a in arrays.items()}),arrays


def main():
    v,t,labels=fixture();first,a=observe(v,t,labels);second,b=observe(v,t,labels)
    gates=dict(full_sequences=not first['refusals'] and first['sequence_errors']==0 and first['maximum_parameter_error']<=1e-12,analytic_exit=not first['refusals'] and first['maximum_exit_error']<=1e-12,full_repeat=first==second and all(x.tobytes()==b[k].tobytes() for k,x in a.items()),partition_accepted=True)
    gates={k:bool(value) for k,value in gates.items()}
    report=dict(rows=[first,second],gates=gates,scope='Frozen straight-ray neighbor walker on a larger synthetic conforming partition; no photon/anatomy promotion.')
    (ROOT/'reports/tetra_ray_grid.json').write_text(json.dumps(report,indent=2)+'\n');np.savez_compressed(ROOT/'reports/tetra_ray_grid_arrays.npz',**{f'{k}_0':v for k,v in a.items()},**{f'{k}_1':v for k,v in b.items()});print(json.dumps(report));return 0 if all(gates.values()) else 2

if __name__=='__main__':raise SystemExit(main())
