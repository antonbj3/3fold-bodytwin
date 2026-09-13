"""Certify unperturbed paths and a fixed detector against external path statistics."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import warp as wp
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from bodytwin.geometry.layered_cylinder_v1 import layered_cylinder
from bodytwin.geometry.optics.nirs_forward_v1 import signal,sao2_hill
from bodytwin.geometry.optics.tissue_photon_paths_v1 import simulate,PACKET
PROPERTIES=np.array([[0.,0.,1.,1.],[.005,1.,.01,1.37],[.01,2.,.8,1.4]])
N=1000000


def statistics(weights,paths):
    total=weights.sum();fraction=total/N;mean=np.sum(weights[:,None]*paths,axis=0)/total
    return {'count':len(weights),'fraction':float(fraction),'fraction_variance':float((np.sum(weights**2)/N-fraction**2)/(N-1)),
        'layer_path_mean':mean.tolist(),'layer_path_variance':(np.sum(weights[:,None]**2*(paths-mean)**2,axis=0)/total**2).tolist()}


def main():
    wp.init();m=layered_cylinder();base=json.loads((ROOT/'reports/photon_curved_layer_l4.json').read_text())['rows'][0]['hashes']
    rows=[];saved={}
    for leg in range(2):
        out=simulate(m['vertices'],m['faces'],PROPERTIES,m['front'],m['back'],source=[32,32,2.0001],direction=[0,0,1],
            initial_region=1,origin=[0,0,0],shape=[64,64,124],reflect=True)
        hashes={k:hashlib.sha256(v.tobytes()).hexdigest() for k,v in out.items()}
        canonical={k:hashes[k] for k in base};escaped=out['terminal'][:,1]>0
        reconstructed=np.floor(PACKET*np.exp(-(out['paths']*PROPERTIES.astype(np.float32).astype(float)[1:,0]).sum(axis=1))+.5).astype(np.int64)
        energy_error=int(np.max(np.abs(reconstructed[escaped]-out['terminal'][escaped,1])))
        distance2=np.sum((out['exits'][:,:3]-[42,32,2])**2,axis=1)
        detected=escaped&(out['exits'][:,6]<=5.)&(distance2<=4.)
        p=out['paths'][detected];w=out['terminal'][detected,1].astype(float)/PACKET
        stat=statistics(w,p) if len(p) else {'count':0}
        nirs=[signal(p,float(sao2_hill(po2)))[0] for po2 in (40.,60.,80.,100.)] if len(p) else []
        row={'hashes':hashes,'canonical_identical':canonical==base,'reconstructed_energy_max_integer_error':energy_error,
             'detector':stat,'nirs':nirs}
        rows.append(row);saved[f'paths_{leg}']=p;saved[f'weights_{leg}']=w;saved[f'exits_{leg}']=out['exits'][detected]
    gates={'canonical_transport_exact':all(r['canonical_identical'] for r in rows),'all_full_arrays_repeat_exact':rows[0]['hashes']==rows[1]['hashes'],
        'nonempty_detector':all(r['detector']['count']>0 for r in rows),
        'path_energy_reconstruction':all(r['reconstructed_energy_max_integer_error']<=1 for r in rows),
        'nirs_direction':all(len(r['nirs'])==4 and all(np.isfinite(x['ratio_760_850']) and x['ratio_760_850']>0 for x in r['nirs']) and all(a['ratio_760_850']<b['ratio_760_850'] for a,b in zip(r['nirs'],r['nirs'][1:])) for r in rows)}
    report={'legs':rows,'gates':gates,'scope':'synthetic curved detector and NIRS direction only; clinical calibration uncertified'}
    (ROOT/'reports/photon_curved_detector.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez_compressed(ROOT/'reports/photon_curved_detector_paths.npz',**saved)
    print(json.dumps(report));return 0 if all(gates.values()) else 1


if __name__=='__main__':raise SystemExit(main())
