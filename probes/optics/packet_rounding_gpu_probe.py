"""Strict integer CPU/GPU comparison at frozen represented rounding boundaries."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from bodytwin.geometry.tetra_ray_walk_v1 import TetraRayWalk
from bodytwin.geometry.optics.tetra_attenuation_v1 import attenuate,PACKET


def main():
    with np.load(ROOT/'reports/packet_rounding_mechanism_arrays.npz') as a:
        lengths=a['lengths'];expected=a['terminal']
    prepared=TetraRayWalk(np.array([[0,0,0],[32,0,0],[0,32,0],[0,0,32]],float),np.array([[0,1,2,3]],np.int32),np.array([1],np.int32))
    n=len(lengths);intervals=np.zeros((n,1,2));intervals[:,0,1]=lengths
    traversal=dict(cells=np.zeros((n,1),np.int32),intervals=intervals,counts=np.ones(n,np.int32),status=np.zeros(n,np.int32))
    rows=[];arrays={}
    for leg in (0,1):
        out=attenuate(prepared,traversal,[0,1]);delta=out['terminal']-expected
        rows.append(dict(leg=leg,differing_rays=int(np.count_nonzero(np.any(delta!=0,axis=1))),maximum_integer_difference=int(np.max(np.abs(delta))),absorbed=int(out['terminal'][:,0].sum()),escaped=int(out['terminal'][:,1].sum()),residual=int(out['terminal'][:,2].sum()),dose_sum=int(out['dose'].sum()),hashes={k:hashlib.sha256(v.tobytes()).hexdigest() for k,v in out.items()}))
        for k,v in out.items():arrays[f'{k}_{leg}']=v
    gates=dict(strict_reference=all(r['differing_rays']==0 and r['dose_sum']==int(expected[:,0].sum()) for r in rows),repeat=rows[0]['hashes']==rows[1]['hashes'],energy=all(r['absorbed']==r['dose_sum'] and r['absorbed']+r['escaped']+r['residual']==n*PACKET and r['residual']==0 for r in rows))
    report=dict(rows=rows,gates=gates,intervals=n,input_sha256=hashlib.sha256(lengths.tobytes()).hexdigest(),scope='Isolated scoring at represented rounding boundaries; no broader CPU/GPU integer identity or whole-geometry claim.')
    (ROOT/'reports/packet_rounding_gpu.json').write_text(json.dumps(report,indent=2)+'\n');np.savez_compressed(ROOT/'reports/packet_rounding_gpu_arrays.npz',**arrays);print(json.dumps(report));return 0 if all(gates.values()) else 2

if __name__=='__main__':raise SystemExit(main())
