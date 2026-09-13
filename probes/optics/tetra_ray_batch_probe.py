"""Paired full-record comparison of shared-origin and frozen scalar traversal."""
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(Path(__file__).resolve().parent))
from bodytwin.geometry.tetra_ray_walk_v1 import TetraRayWalk
from bodytwin.geometry.tetra_ray_batch_v1 import TetraRayBatch
from tetra_ray_grid_probe import fixture
from tetra_partition_mechanism import fixtures
from tetra_ray_mechanism import intervals,observe


def packed(records):
    return dict(records=np.concatenate(records),offsets=np.r_[0,np.cumsum([len(a) for a in records])].astype(np.int64))


def digest(arrays):
    return {k:hashlib.sha256(a.tobytes()).hexdigest() for k,a in arrays.items()}


def main():
    v,t,labels=fixture();origin=np.array([.17,.23,.31])
    directions=np.random.default_rng(20260913).normal(size=(512,3))
    directions/=np.linalg.norm(directions,axis=1)[:,None]
    before=[a.tobytes() for a in (v,t,labels,origin,directions)]
    base=TetraRayWalk(v,t,labels);candidate=TetraRayBatch(v,t,labels)
    baseline=lambda:packed([base.trace(origin,d) for d in directions])
    batch=lambda:packed(candidate.trace_batch(origin,directions))
    baseline();batch();rows=[];saved={}
    for leg,order in enumerate((('baseline','batch'),('batch','baseline'))):
        timings={};outputs={}
        for name in order:
            start=time.perf_counter();outputs[name]=(baseline if name=='baseline' else batch)();timings[name]=time.perf_counter()-start
        rows.append(dict(leg=leg,seconds=timings,hashes={k:digest(a) for k,a in outputs.items()},speedup=timings['baseline']/timings['batch']))
        for name,arrays in outputs.items():
            for key,a in arrays.items():saved[f'{name}_{key}_{leg}']=a
    _,sv,st,sl=fixtures()[1];_,ref=observe();controls=[];oracle_error=0.;oracle_sequence=True
    for leg in range(2):
        for name,cells in (('original',st),('inverted',st[:,::-1])):
            b=TetraRayWalk(sv,cells,sl);c=TetraRayBatch(sv,cells,sl)
            for launch in (ref['origin'],sv[cells[0]].mean(axis=0)):
                # The alternate launch also checks that location is recomputed between calls.
                ds=ref['directions'] if np.array_equal(launch,ref['origin']) else directions[:32]
                got=c.trace_batch(launch,ds);expected=[b.trace(launch,d) for d in ds]
                controls.append(all(x.tobytes()==y.tobytes() for x,y in zip(got,expected)))
                for d,a in zip(ds,got):
                    oracle,_=intervals(sv,cells,launch,d)
                    oracle_sequence &= a.shape==oracle.shape and np.array_equal(a[:,0],oracle[:,0])
                    if a.shape==oracle.shape:oracle_error=max(oracle_error,float(np.max(np.abs(a[:,1:]-oracle[:,1:]))))
                    else:oracle_error=float('inf')
            cases=[([0,0,0],[1,0,0]),(ref['origin'],[-1,-1,0]),(ref['origin'],[0,0,0]),(ref['origin'],[np.nan,0,1]),([100,100,100],[1,0,0])]
            for o,d in cases:
                refused=[]
                for fn in (lambda:b.trace(o,d),lambda:c.trace_batch(o,[d])):
                    try:fn();refused.append(False)
                    except ValueError:refused.append(True)
                controls.append(all(refused))
            controls.append(c.trace_batch(ref['origin'],np.empty((0,3)))==[])
    gates=dict(full_baseline_equality=all(r['hashes']['baseline']==r['hashes']['batch'] for r in rows),
        full_repeat=rows[0]['hashes']==rows[1]['hashes'],controls=all(controls),
        oracle=bool(oracle_sequence and oracle_error<=1e-12),
        inputs_unchanged=before==[a.tobytes() for a in (v,t,labels,origin,directions)],
        paired_speedup=all(r['speedup']>1 for r in rows))
    report=dict(rows=rows,cells=len(t),rays=len(directions),controls_passed=sum(controls),controls_total=len(controls),maximum_oracle_error=oracle_error,gates=gates,
        source_sha256={n:hashlib.sha256((ROOT/'src/bodytwin/geometry'/n).read_bytes()).hexdigest() for n in ('tetra_ray_walk_v1.py','tetra_ray_batch_v1.py')},
        scope='Synthetic CPU shared-origin straight-ray traversal only; no CUDA, photon, anatomy or clinical claim.')
    (ROOT/'reports/tetra_ray_batch.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez_compressed(ROOT/'reports/tetra_ray_batch_arrays.npz',**saved)
    print(json.dumps(report));return 0 if all(gates.values()) else 2

if __name__=='__main__':raise SystemExit(main())
