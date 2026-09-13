"""Full boundary/slab gates and byte-preserving backend package regression."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import trimesh
import warp as wp
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(Path(__file__).resolve().parent))
from bodytwin.geometry.optics.tissue_photon_mc_v3 import simulate,PACKET
from mcx_boundary_references import slab_profile,diffusion_profile


def hashes(result):return {k:hashlib.sha256(np.ascontiguousarray(v).tobytes()).hexdigest() for k,v in result.items()}


def coarse(a):
    n=a.shape[0]//5
    return a.reshape(n,5,n,5,n,5).sum(axis=(1,3,5))


def main():
    wp.init();references=np.load(ROOT/'reports/mcx_boundary_reference_fields_l4.npz')
    original=json.loads((ROOT/'reports/photon_cube_v3_l4.json').read_text())['legs'][0]['hashes']
    analytic=diffusion_profile();cases=[];all_gates=[]
    for case in ('cube60','cube60b','slab'):
        side=100 if case=='slab' else 60;mua=.01 if case=='slab' else .005
        medium=[.01,10,.9,1] if case=='slab' else [.005,1,.01,1.37]
        xy=50 if case=='slab' else 29
        mesh=trimesh.creation.box(extents=[side]*3);mesh.apply_translation([side/2]*3)
        outputs=[];rows=[]
        for _ in range(2):
            out=simulate(mesh.vertices,mesh.faces,[[0,0,1,1],medium],np.zeros(12,dtype=int),np.ones(12,dtype=int),
                source=[xy,xy,.0001],direction=[0,0,1],initial_region=1,origin=[0,0,0],shape=[side]*3,reflect=case=='cube60b')
            outputs.append(out);sums=out['terminal'][:,:3].sum(axis=0,dtype=np.int64)
            fluence=out['absorption'].astype(float)/(1000000*PACKET*mua)
            row={'hashes':hashes(out),'absorbed_integer':int(sums[0]),'escaped_integer':int(sums[1]),'residual_integer':int(sums[2]),
                 'leaked_photons':int(out['terminal'][:,3].sum()),'event_cap_photons':int(out['counters'][2]),
                 'events':int(out['counters'][0]),'boundary_hits':int(out['counters'][1]),
                 'finite_nonnegative':bool(np.isfinite(fluence).all() and fluence.min()>=0),
                 'window_absorbed_fraction':float(fluence.sum()*mua)}
            if case!='cube60':
                refs=[references[f'{case}_{leg}'].astype(float) for leg in (0,1)]
                row['absorption_relative_errors']=[float(abs(fluence.sum()-a.sum())/a.sum()) for a in refs]
                row['binned_fluence_l1_errors']=[float(np.abs(coarse(fluence)-coarse(a)).sum()/a.sum()) for a in refs]
            if case=='slab':
                profile=slab_profile(fluence)
                row['profile']=profile.tolist();row['diffusion_relative_errors']=(np.abs(profile-analytic)/analytic).tolist()
                row['reference_profile_errors']=[float(np.max(np.abs(profile-slab_profile(a))/slab_profile(a))) for a in refs]
            rows.append(row);print(json.dumps({'case':case,**row}),flush=True)
        gates={'full_arrays_bit_identical':all(outputs[0][k].tobytes()==outputs[1][k].tobytes() for k in outputs[0]),
               'exact_energy':all(r['absorbed_integer']+r['escaped_integer']==1000000*PACKET for r in rows),
               'no_failed_packets':all(r['residual_integer']==0 and r['leaked_photons']==0 and r['event_cap_photons']==0 for r in rows),
               'finite_nonnegative':all(r['finite_nonnegative'] for r in rows)}
        if case=='cube60':gates['pre_move_hashes_exact']=all(r['hashes']==original for r in rows)
        else:
            gates['reference_absorption']=all(max(r['absorption_relative_errors'])<=.02 for r in rows)
            gates['reference_binned_fluence']=all(max(r['binned_fluence_l1_errors'])<=.05 for r in rows)
        if case=='slab':
            gates['diffusion_profile']=all(max(r['diffusion_relative_errors'])<=.05 for r in rows)
            gates['reference_profile']=all(max(r['reference_profile_errors'])<=.05 for r in rows)
        cases.append({'case':case,'legs':rows,'gates':gates});all_gates.extend(gates.values())
    report={'cases':cases,'all_pass':all(all_gates),'scope':'fixed homogeneous/refractive/slab fixtures, not arbitrary anatomy or clinical calibration'}
    (ROOT/'reports/photon_boundary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True);return 0 if report['all_pass'] else 1


if __name__=='__main__':raise SystemExit(main())
