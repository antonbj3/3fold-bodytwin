"""Frozen v6 exact transport controls on a separately measured curved fixture."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import trimesh
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from bodytwin.geometry.layered_cylinder_v1 import layered_cylinder


def main():
    import warp as wp
    from bodytwin.geometry.optics.tissue_photon_mc_v6 import simulate,PACKET
    mesh=layered_cylinder();other=layered_cylinder()
    geometry_repeat=all(mesh[k].tobytes()==other[k].tobytes() for k in mesh)
    volumes=[];closed=True
    for region in (1,2):
        faces=np.concatenate((mesh['faces'][mesh['back']==region],mesh['faces'][mesh['front']==region][:,::-1]))
        surface=trimesh.Trimesh(mesh['vertices'],faces,process=False)
        closed &= bool(surface.is_watertight and surface.is_winding_consistent and surface.volume>0)
        volumes.append(float(surface.volume))
    ideal=[np.pi*(30**2*120-27**2*114),np.pi*27**2*114]
    errors=[abs(a-b)/b for a,b in zip(volumes,ideal)]
    wp.init();rows=[];outputs=[]
    properties=np.array([[0.,0.,1.,1.],[.005,1.,.01,1.37],[.01,2.,.8,1.4]])
    for _ in range(2):
        out=simulate(mesh['vertices'],mesh['faces'],properties,mesh['front'],mesh['back'],
                     source=[32.,32.,2.0001],direction=[0,0,1],initial_region=1,
                     origin=[0,0,0],shape=[64,64,124],reflect=True)
        outputs.append(out);sums=out['terminal'][:,:3].sum(axis=0,dtype=np.int64)
        rows.append(dict(hashes={k:hashlib.sha256(v.tobytes()).hexdigest() for k,v in out.items()},
                         absorbed=int(sums[0]),escaped=int(sums[1]),residual=int(sums[2]),
                         leaks=int(out['terminal'][:,3].sum()),caps=int(out['counters'][2]),
                         finite_nonnegative=bool(np.isfinite(out['absorption']).all() and np.min(out['absorption'])>=0)))
    gates=dict(geometry_repeat=geometry_repeat,closed_regions=closed,volume_bound=max(errors)<=.002,
               exact_energy=all(r['absorbed']+r['escaped']==1000000*PACKET for r in rows),
               zero_failed_packets=all(r['residual']==r['leaks']==r['caps']==0 for r in rows),
               full_repeat=all(outputs[0][k].tobytes()==outputs[1][k].tobytes() for k in outputs[0]),
               finite_nonnegative=all(r['finite_nonnegative'] for r in rows))
    report=dict(rows=rows,gates=gates,region_volumes=volumes,relative_volume_errors=errors,
                backend_sha256=hashlib.sha256((ROOT/'src/bodytwin/geometry/optics/tissue_photon_mc_v6.py').read_bytes()).hexdigest(),
                scope='Synthetic curved transport control, not anatomical or detector calibration.')
    (ROOT/'reports/photon_curved_layer.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report));return 0 if all(gates.values()) else 2


if __name__=='__main__':raise SystemExit(main())
