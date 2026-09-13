"""Exhaustive cell interval oracle before neighbor traversal design."""
import hashlib
import json
from pathlib import Path
import numpy as np
from tetra_partition_mechanism import fixtures
ROOT=Path(__file__).resolve().parents[2]


def intervals(v,t,origin,direction):
    p=v[t];matrix=np.stack((p[:,1]-p[:,0],p[:,2]-p[:,0],p[:,3]-p[:,0]),axis=2)
    inverse=np.linalg.inv(matrix);local=np.einsum('nij,nj->ni',inverse,origin-p[:,0]);slope=np.einsum('nij,j->ni',inverse,direction)
    a=np.column_stack((1-local.sum(axis=1),local));b=np.column_stack((-slope.sum(axis=1),slope))
    result=[]
    for cell,(offset,delta) in enumerate(zip(a,b)):
        lo=0.;hi=float('inf');valid=True
        for value,derivative in zip(offset,delta):
            if derivative>0:lo=max(lo,-value/derivative)
            elif derivative<0:hi=min(hi,-value/derivative)
            elif value<0:valid=False
        if valid and hi>lo:result.append((cell,lo,hi))
    return np.asarray(sorted(result,key=lambda r:(r[1],r[0])),dtype=np.float64),a


def observe():
    _,v,t,labels=fixtures()[1];origin=v[t[0]].mean(axis=0)
    directions=np.random.default_rng(20260913).normal(size=(512,3));directions/=np.linalg.norm(directions,axis=1)[:,None]
    normals=np.array([[x,y,z] for x in (-1,1) for y in (-1,1) for z in (-1,1)])
    rows=[];trace=[];offsets=[0]
    for direction in directions:
        values,a=intervals(v,t,origin,direction)
        projection=normals@direction;exits=(1-normals@origin)[projection>0]/projection[projection>0];analytic=float(exits.min())
        gap=float(np.max(np.abs(values[1:,1]-values[:-1,2]))) if len(values)>1 else 0.
        rows.append(dict(segments=len(values),min_length=float(np.min(values[:,2]-values[:,1])),gap_or_overlap=gap,exit_error=abs(float(values[-1,2])-analytic),strict_source_cells=int(np.count_nonzero(np.min(a,axis=1)>0))))
        trace.extend(values.tolist());offsets.append(len(trace))
    _,boundary=intervals(v,t,np.zeros(3),np.array([1.,0.,0.]));ambiguous=int(np.count_nonzero(np.min(boundary,axis=1)>=0))
    arrays=dict(directions=directions,intervals=np.asarray(trace),offsets=np.asarray(offsets,np.int64),origin=origin)
    report=dict(rays=len(rows),mean_segments=float(np.mean([r['segments'] for r in rows])),maximum_segments=max(r['segments'] for r in rows),minimum_length=min(r['min_length'] for r in rows),maximum_gap_or_overlap=max(r['gap_or_overlap'] for r in rows),maximum_exit_error=max(r['exit_error'] for r in rows),strict_source_unique=all(r['strict_source_cells']==1 for r in rows),edge_source_closed_cells=ambiguous,hashes={k:hashlib.sha256(a.tobytes()).hexdigest() for k,a in arrays.items()})
    return report,arrays


def main():
    first,arrays=observe();second,again=observe()
    gates=dict(repeat=first==second and all(a.tobytes()==again[k].tobytes() for k,a in arrays.items()),coverage=first['strict_source_unique'] and first['maximum_gap_or_overlap']<=1e-12 and first['maximum_exit_error']<=1e-12,ambiguous_refused=first['edge_source_closed_cells']>1)
    report=dict(measurement=first,gates=gates,scope='Exhaustive synthetic straight-ray interval oracle; ambiguous source identified, no transport acceptance.')
    (ROOT/'reports/tetra_ray_mechanism.json').write_text(json.dumps(report,indent=2)+'\n');np.savez_compressed(ROOT/'reports/tetra_ray_mechanism_arrays.npz',**arrays);print(json.dumps(report));return 0 if all(gates.values()) else 2

if __name__=='__main__':raise SystemExit(main())
