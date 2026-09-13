"""Capture frozen scattering/reflection arrays; no statistical-reference promotion."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import trimesh
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from bodytwin.geometry.optics.tissue_photon_mc_v6 import simulate,PACKET


def main():
    mesh=trimesh.creation.box(extents=[60]*3);mesh.apply_translation([30]*3)
    rows=[];arrays={};photons=10000
    for reflect in (False,True):
        for leg in range(2):
            out=simulate(mesh.vertices,mesh.faces,[[0,0,1,1],[.005,1,.01,1.37]],np.zeros(12,dtype=int),np.ones(12,dtype=int),source=[29,29,.0001],direction=[0,0,1],initial_region=1,origin=[0,0,0],shape=[60]*3,photons=photons,reflect=reflect)
            sums=out['terminal'][:,:3].sum(axis=0,dtype=np.int64)
            rows.append(dict(reflect=reflect,leg=leg,absorbed=int(sums[0]),escaped=int(sums[1]),residual=int(sums[2]),leaks=int(out['terminal'][:,3].sum()),caps=int(out['counters'][2]),events=int(out['counters'][0]),boundary_hits=int(out['counters'][1]),finite_nonnegative=bool(all(np.isfinite(a).all() and np.min(a)>=0 for a in out.values())),hashes={k:hashlib.sha256(a.tobytes()).hexdigest() for k,a in out.items()}))
            for k,a in out.items():arrays[f'{int(reflect)}_{leg}_{k}']=a
    gates=dict(energy=all(r['absorbed']+r['escaped']==photons*PACKET for r in rows),no_failed=all(r['residual']==r['leaks']==r['caps']==0 for r in rows),repeat=all(rows[i]['hashes']==rows[i+1]['hashes'] for i in (0,2)),finite=all(r['finite_nonnegative'] for r in rows))
    report=dict(rows=rows,gates=gates,photons=photons,source_sha256=hashlib.sha256((ROOT/'src/bodytwin/geometry/optics/tissue_photon_mc_v6.py').read_bytes()).hexdigest(),scope='Frozen scattering cube, reflection off/on, bounded determinism observer; no new MCX or anatomy validation.')
    (ROOT/'reports/photon_cross_capture.json').write_text(json.dumps(report,indent=2)+'\n');np.savez_compressed(ROOT/'reports/photon_cross_capture_arrays.npz',**arrays);print(json.dumps(report));return 0 if all(gates.values()) else 2

if __name__=='__main__':raise SystemExit(main())
