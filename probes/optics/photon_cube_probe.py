"""Two full integer mesh-transport runs against both frozen cube60 references."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import trimesh
import warp as wp
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from bodytwin.geometry.optics.tissue_photon_mc_v1 import simulate,PACKET


def hashes(result):
    return {k:hashlib.sha256(np.ascontiguousarray(v).tobytes()).hexdigest() for k,v in result.items()}


def coarse(a):
    return a.reshape(12,5,12,5,12,5).sum(axis=(1,3,5))


def main():
    wp.init();mesh=trimesh.creation.box(extents=[60,60,60]);mesh.apply_translation([30,30,30])
    reference=np.load(ROOT/'reports/mcx_reference_fields_l4.npz')
    refs=[reference[k][...,0].astype(float)*5e-9*.005 for k in ('first','second')]
    outputs=[];rows=[]
    for _ in range(2):
        result=simulate(mesh.vertices,mesh.faces,[[0,0,1,1],[.005,1,.01,1.37]],
                        np.zeros(len(mesh.faces),dtype=int),np.ones(len(mesh.faces),dtype=int),
                        source=[29,29,.0001],direction=[0,0,1],initial_region=1,origin=[0,0,0],shape=[60,60,60])
        outputs.append(result)
        hist=result['absorption'].astype(float)/(1000000*PACKET)
        sums=result['terminal'][:,:3].sum(axis=0,dtype=np.int64)
        row={'photons':1000000,'absorbed_integer':int(sums[0]),'escaped_integer':int(sums[1]),'residual_integer':int(sums[2]),
             'launched_integer':1000000*PACKET,'leaked_photons':int(result['terminal'][:,3].sum()),
             'event_cap_photons':int(result['counters'][2]),'events':int(result['counters'][0]),'boundary_hits':int(result['counters'][1]),
             'window_absorbed_fraction':float(hist.sum()),'all_time_absorbed_fraction':float(sums[0]/(1000000*PACKET)),
             'relative_absorption_errors':[float(abs(hist.sum()-a.sum())/a.sum()) for a in refs],
             'coarse_fluence_l1_errors':[float(np.abs(coarse(hist)-coarse(a)).sum()/a.sum()) for a in refs],
             'finite_nonnegative':bool(np.isfinite(hist).all() and hist.min()>=0),'hashes':hashes(result)}
        rows.append(row);print(json.dumps(row),flush=True)
    gates={'full_arrays_bit_identical':all(outputs[0][k].tobytes()==outputs[1][k].tobytes() for k in outputs[0]),
           'absorbed_plus_escaped_equals_launched':all(r['absorbed_integer']+r['escaped_integer']==r['launched_integer'] for r in rows),
           'no_residual_or_leak':all(r['residual_integer']==0 and r['leaked_photons']==0 and r['event_cap_photons']==0 for r in rows),
           'finite_nonnegative':all(r['finite_nonnegative'] for r in rows),
           'reference_absorption':all(max(r['relative_absorption_errors'])<=.02 for r in rows),
           'reference_binned_fluence':all(max(r['coarse_fluence_l1_errors'])<=.05 for r in rows)}
    report={'legs':rows,'gates':gates,'scope':'cube reference only; engineering comparison bands, no slab/clinical calibration claim',
            'source_inward_offset_mm':.0001,'interface_offset_mm':.00001,'packet_scale':PACKET,'max_events':50000}
    (ROOT/'reports/photon_cube.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez_compressed(ROOT/'reports/photon_cube_fields.npz',first=outputs[0]['absorption'],second=outputs[1]['absorption'])
    print(json.dumps(report),flush=True)
    return 0 if all(gates.values()) else 1


if __name__=='__main__':raise SystemExit(main())
