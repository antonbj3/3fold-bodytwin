"""Measured integer-energy seam on frozen cell traversal arrays."""
import hashlib
import json
from pathlib import Path
import numpy as np
from tetra_ray_grid_probe import fixture
ROOT=Path(__file__).resolve().parents[2];PACKET=1<<30


def reference(cells,distances,counts,status,labels,mu):
    dose=np.zeros(len(labels),np.int64);terminal=np.zeros((len(counts),3),np.int64);depths=[]
    for ray,count in enumerate(counts):
        if status[ray]!=0:terminal[ray,2]=PACKET;depths.append(0.);continue
        continuous=float(PACKET);weight=PACKET;depth=0.
        for slot in range(int(count)):
            cell=int(cells[ray,slot]);step=float(mu[labels[cell]])*float(distances[ray,slot,1]-distances[ray,slot,0]);depth+=step
            continuous*=np.exp(-step);after=int(np.floor(continuous+.5));dose[cell]+=weight-after;weight=after
        terminal[ray,0]=PACKET-weight;terminal[ray,1]=weight;depths.append(depth)
    return dict(dose=dose,terminal=terminal),np.asarray(depths)


def main():
    _,_,labels=fixture();rows=[];saved={}
    with np.load(ROOT/'reports/tetra_ray_cuda_arrays.npz') as a:traversal={k:a[f'{k}_0'] for k in ('cells','intervals','counts','status')}
    for name,mu in [('transparent',[0,0,0]),('heterogeneous',[0,.01,.04]),('strong',[0,10,20])]:
        for leg in (0,1):
            out,depth=reference(traversal['cells'],traversal['intervals'],traversal['counts'],traversal['status'],labels,np.asarray(mu))
            rows.append(dict(case=name,leg=leg,depth_min=float(depth.min()),depth_max=float(depth.max()),zero_escape=int(np.count_nonzero(out['terminal'][:,1]==0)),absorbed=int(out['terminal'][:,0].sum()),escaped=int(out['terminal'][:,1].sum()),nonnegative=bool(all(np.min(x)>=0 for x in out.values())),hashes={k:hashlib.sha256(v.tobytes()).hexdigest() for k,v in out.items()}))
            for k,v in out.items():saved[f'{name}_{leg}_{k}']=v
    gates=dict(energy=all(r['absorbed']+r['escaped']==512*PACKET for r in rows),nonnegative=all(r['nonnegative'] for r in rows),repeat=all(rows[i]['hashes']==rows[i+1]['hashes'] for i in (0,2,4)))
    report=dict(rows=rows,gates=gates,scope='Synthetic straight-path attenuation seam; no scattering/refraction or physical dose calibration.')
    (ROOT/'reports/tetra_attenuation_mechanism.json').write_text(json.dumps(report,indent=2)+'\n');np.savez_compressed(ROOT/'reports/tetra_attenuation_mechanism_arrays.npz',**saved);print(json.dumps(report));return 0 if all(gates.values()) else 2

if __name__=='__main__':raise SystemExit(main())
