"""Bounded hundred-million-packet observer over the unchanged tissue transport."""
import gc
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import trimesh
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from bodytwin.geometry.optics.tissue_photon_mc_v6 import simulate,PACKET


def coarse(a):
    return a.reshape(12,5,12,5,12,5).sum(axis=(1,3,5))


def main():
    photons=100000000;mesh=trimesh.creation.box(extents=[60]*3);mesh.apply_translation([30]*3)
    reference_file=ROOT/'reports/mcx_boundary_reference_fields_l4.npz'
    with np.load(reference_file) as source:references=[source[f'cube60b_{leg}'].astype(float) for leg in (0,1)]
    rows=[];saved={}
    for leg in range(2):
        out=simulate(mesh.vertices,mesh.faces,[[0,0,1,1],[.005,1,.01,1.37]],np.zeros(12,dtype=int),np.ones(12,dtype=int),source=[29,29,.0001],direction=[0,0,1],initial_region=1,origin=[0,0,0],shape=[60]*3,photons=photons,reflect=True)
        sums=out['terminal'][:,:3].sum(axis=0,dtype=np.int64);fluence=out['absorption'].astype(float)/(photons*PACKET*.005)
        row=dict(leg=leg,absorbed=int(sums[0]),escaped=int(sums[1]),residual=int(sums[2]),leaks=int(out['terminal'][:,3].sum()),caps=int(out['counters'][2]),events=int(out['counters'][0]),boundary_hits=int(out['counters'][1]),finite_nonnegative=bool(all(np.isfinite(a).all() and np.min(a)>=0 for a in out.values())),hashes={k:hashlib.sha256(memoryview(a).cast('B')).hexdigest() for k,a in out.items()},reference_absorption_errors=[float(abs(fluence.sum()-a.sum())/a.sum()) for a in references],reference_binned_errors=[float(np.abs(coarse(fluence)-coarse(a)).sum()/a.sum()) for a in references])
        rows.append(row);saved[f'absorption_{leg}']=out['absorption'];saved[f'counters_{leg}']=out['counters']
        # Save each completed leg before releasing its large terminal ledger.
        (ROOT/'reports/photon_scale_partial.json').write_text(json.dumps(dict(photons=photons,rows=rows),indent=2)+'\n')
        print(json.dumps(row),flush=True);del out;gc.collect()
    gates=dict(exact_energy=all(r['absorbed']+r['escaped']==photons*PACKET for r in rows),no_failed=all(r['residual']==r['leaks']==r['caps']==0 for r in rows),full_repeat=rows[0]['hashes']==rows[1]['hashes'],finite=all(r['finite_nonnegative'] for r in rows),reference_absorption=all(max(r['reference_absorption_errors'])<=.02 for r in rows),reference_binned=all(max(r['reference_binned_errors'])<=.05 for r in rows))
    report=dict(photons_per_leg=photons,rows=rows,gates=gates,source_sha256=hashlib.sha256((ROOT/'src/bodytwin/geometry/optics/tissue_photon_mc_v6.py').read_bytes()).hexdigest(),reference_sha256=hashlib.sha256(reference_file.read_bytes()).hexdigest(),scope='Synthetic reflective cube at100000000 independent packet indices per leg, fixed tissue transport and MCX references; no anatomy or throughput claim.')
    (ROOT/'reports/photon_scale.json').write_text(json.dumps(report,indent=2)+'\n');np.savez_compressed(ROOT/'reports/photon_scale_arrays.npz',**saved);print(json.dumps(gates));return 0 if all(gates.values()) else 2

if __name__=='__main__':raise SystemExit(main())
