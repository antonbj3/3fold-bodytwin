"""Complete GPU cell-energy comparison and incomplete-path accounting."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(Path(__file__).resolve().parent))
from bodytwin.geometry.optics.tetra_attenuation_v1 import attenuate,PACKET
from bodytwin.geometry.tetra_ray_walk_v1 import TetraRayWalk
from tetra_ray_grid_probe import fixture
from tetra_attenuation_mechanism import reference


def main():
    v,t,labels=fixture();prepared=TetraRayWalk(v,t,labels);rows=[];arrays={}
    with np.load(ROOT/'reports/tetra_ray_cuda_arrays.npz') as a:
        full={k:a[f'{k}_0'] for k in ('cells','intervals','counts','status')};capped={k:a[f'capped_{k}_0'] for k in full}
    for name,mu,path in [('transparent',[0,0,0],full),('heterogeneous',[0,.01,.04],full),('strong',[0,10,20],full),('capped',[0,.01,.04],capped)]:
        expected,_=reference(path['cells'],path['intervals'],path['counts'],path['status'],labels,np.asarray(mu))
        for leg in (0,1):
            out=attenuate(prepared,path,mu);total=out['terminal'].sum(axis=0,dtype=np.int64)
            rows.append(dict(case=name,leg=leg,absorbed=int(total[0]),escaped=int(total[1]),residual=int(total[2]),dose_sum=int(out['dose'].sum()),reference_exact=all(out[k].tobytes()==expected[k].tobytes() for k in out),zero_escape=int(np.count_nonzero(out['terminal'][:,1]==0)),hashes={k:hashlib.sha256(a.tobytes()).hexdigest() for k,a in out.items()}))
            for k,a in out.items():arrays[f'{name}_{leg}_{k}']=a
    failed=capped['status']!=0
    gates=dict(reference=all(r['reference_exact'] for r in rows),repeat=all(rows[i]['hashes']==rows[i+1]['hashes'] for i in (0,2,4,6)),energy=all(r['absorbed']==r['dose_sum'] and r['absorbed']+r['escaped']+r['residual']==512*PACKET for r in rows),transparent=all(r['absorbed']==r['residual']==0 and r['escaped']==512*PACKET for r in rows if r['case']=='transparent'),incomplete=all(np.all(arrays[f'capped_{leg}_terminal'][failed,:2]==0) and np.all(arrays[f'capped_{leg}_terminal'][failed,2]==PACKET) for leg in (0,1)))
    report=dict(rows=rows,gates=gates,failed_trace_count=int(failed.sum()),scope='Synthetic Beer-Lambert scoring of known straight cell paths, no scattering/refraction or clinical dose.')
    (ROOT/'reports/tetra_attenuation.json').write_text(json.dumps(report,indent=2)+'\n');np.savez_compressed(ROOT/'reports/tetra_attenuation_arrays.npz',**arrays);print(json.dumps(report));return 0 if all(gates.values()) else 2

if __name__=='__main__':raise SystemExit(main())
